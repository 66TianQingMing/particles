import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import image as mpimg
from matplotlib.animation import FFMpegWriter
from matplotlib.patches import Ellipse


FLUID = 1
BOUNDARY = 2
PARTICLE = 3


def read_flow_field(file_name):
    particles = []
    field_header = False
    grid_size = None
    with open(file_name, "r", encoding="utf-8") as file:
        for line in file:
            if not line.startswith("#"):
                break
            words = line.split()
            if words[1:] == ["X", "Y", "Type", "Density", "Vx", "Vy"]:
                field_header = True
            if len(words) == 4 and words[1] == "GridSize":
                grid_size = (int(words[2]), int(words[3]))
            if len(words) == 13 and words[1] == "Particle":
                particles.append({
                    "id": int(words[2]),
                    "sx": float(words[4]),
                    "sy": float(words[6]),
                    "a": float(words[8]),
                    "b": float(words[10]),
                    "fai": float(words[12]),
                })

    if not field_header:
        raise ValueError("not a flow-field snapshot")

    data = np.atleast_2d(np.loadtxt(file_name, comments="#"))
    if grid_size is not None and len(data) != grid_size[0]*grid_size[1]:
        raise ValueError(
            "incomplete file: %d rows, but GridSize %d x %d requires %d rows"
            % (len(data), grid_size[0], grid_size[1], grid_size[0]*grid_size[1])
        )
    return data, particles


def reshape_field(data):
    x_value = np.unique(data[:, 0])
    y_value = np.unique(data[:, 1])
    nx = len(x_value)
    ny = len(y_value)

    if nx*ny != len(data):
        raise ValueError(
            "incomplete rectangular grid: %d rows, but %d x %d = %d rows are required"
            % (len(data), nx, ny, nx*ny)
        )

    x = data[:, 0].reshape(nx, ny).T
    y = data[:, 1].reshape(nx, ny).T
    node_type = data[:, 2].reshape(nx, ny).T
    density = data[:, 3].reshape(nx, ny).T
    vx = data[:, 4].reshape(nx, ny).T
    vy = data[:, 5].reshape(nx, ny).T
    return x, y, node_type, density, vx, vy


def get_field_value(file_name, show_density):
    data, particles = read_flow_field(file_name)
    x, y, node_type, density, vx, vy = reshape_field(data)
    value = density if show_density else np.sqrt(vx*vx+vy*vy)
    return particles, x, y, node_type, density, vx, vy, value


def find_color_range(file_names, show_density):
    valid_files = []
    value_min = None
    value_max = None

    for file_name in file_names:
        try:
            particles, x, y, node_type, density, vx, vy, value = get_field_value(file_name, show_density)
            fluid_value = value[node_type == FLUID]
            fluid_value = fluid_value[np.isfinite(fluid_value)]
            if len(fluid_value) == 0:
                raise ValueError("no finite fluid value")
            value_min = fluid_value.min() if value_min is None else min(value_min, fluid_value.min())
            value_max = fluid_value.max() if value_max is None else max(value_max, fluid_value.max())
            valid_files.append(file_name)
        except ValueError as error:
            print("Skipped %s: %s" % (file_name, error))

    if value_min is not None and value_min == value_max:
        margin = max(abs(value_min)*1E-6, 1E-12)
        value_min -= margin
        value_max += margin
    return valid_files, value_min, value_max


def plot_flow_field(file_name, output_name, arrow_step, show_density, value_min, value_max):
    particles, x, y, node_type, density, vx, vy, value = get_field_value(file_name, show_density)
    value = np.ma.masked_where(node_type != FLUID, value)

    figure, axis = plt.subplots(figsize=(14, 4.5))
    image = axis.pcolormesh(x, y, value, shading="nearest", cmap="turbo", vmin=value_min, vmax=value_max)
    color_name = "Density" if show_density else "Velocity magnitude"
    figure.colorbar(image, ax=axis, label=color_name)

    wall = np.ma.masked_where(node_type != BOUNDARY, node_type)
    solid = np.ma.masked_where(node_type != PARTICLE, node_type)
    axis.pcolormesh(x, y, wall, shading="nearest", cmap="Greys", vmin=0, vmax=3)
    axis.pcolormesh(x, y, solid, shading="nearest", cmap="Greys", vmin=0, vmax=4)

    step = max(1, arrow_step)
    arrow_x = x[::step, ::step]
    arrow_y = y[::step, ::step]
    arrow_vx = np.ma.masked_where(node_type[::step, ::step] != FLUID, vx[::step, ::step])
    arrow_vy = np.ma.masked_where(node_type[::step, ::step] != FLUID, vy[::step, ::step])
    axis.quiver(arrow_x, arrow_y, arrow_vx, arrow_vy, color="black", alpha=0.65)

    for particle in particles:
        ellipse = Ellipse(
            (particle["sx"], particle["sy"]),
            2*particle["a"],
            2*particle["b"],
            angle=math.degrees(particle["fai"]),
            fill=False,
            edgecolor="red",
            linewidth=1.5,
        )
        axis.add_patch(ellipse)
        axis.text(particle["sx"], particle["sy"], "P%d" % particle["id"], color="red")

    axis.set_aspect("equal")
    axis.set_xlabel("X")
    axis.set_ylabel("Y")
    axis.set_title(Path(file_name).stem)
    figure.tight_layout()
    figure.savefig(output_name, dpi=200)
    plt.close(figure)
    print("Saved:", output_name)


def create_video(image_names, output_name, fps):
    if len(image_names) == 0:
        return False
    if not FFMpegWriter.isAvailable():
        print("Video skipped: FFmpeg is not installed or is not available in PATH.")
        return False

    first_frame = mpimg.imread(str(image_names[0]))
    height, width = first_frame.shape[0], first_frame.shape[1]
    dpi = 100
    figure = plt.figure(figsize=(width/float(dpi), height/float(dpi)), dpi=dpi)
    axis = figure.add_axes([0, 0, 1, 1])
    axis.set_axis_off()
    frame_image = axis.imshow(first_frame)

    writer = FFMpegWriter(fps=fps, codec="h264", extra_args=["-pix_fmt", "yuv420p"])
    with writer.saving(figure, str(output_name), dpi=dpi):
        for image_name in image_names:
            frame_image.set_data(mpimg.imread(str(image_name)))
            writer.grab_frame()
    plt.close(figure)
    print("Saved:", output_name)
    return True


def main():
    parser = argparse.ArgumentParser(description="Plot all particle and flow-field snapshots.")
    parser.add_argument("--input-dir", default="result", help="Directory containing result data files")
    parser.add_argument("--output-dir", default="result/figure", help="Directory for generated images")
    parser.add_argument("--arrow-step", type=int, default=10, help="Grid interval between velocity arrows")
    parser.add_argument("--fps", type=int, default=10, help="Video frames per second")
    parser.add_argument("--density", action="store_true", help="Plot density instead of velocity magnitude")
    parser.add_argument("--no-video", action="store_true", help="Generate PNG images without an MP4 video")
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for image_name in output_dir.glob("*.png"):
        image_name.unlink()
    video_name = output_dir / "simulation.mp4"
    if video_name.exists():
        video_name.unlink()

    file_names = sorted(input_dir.glob("*.dat"))
    valid_files, value_min, value_max = find_color_range(file_names, args.density)
    if len(valid_files) == 0:
        print("Finished: no complete flow-field snapshot was found in %s" % input_dir)
        return

    print("Color range: %.12g to %.12g" % (value_min, value_max))
    image_names = []
    for file_name in valid_files:
        output_name = output_dir / (file_name.stem+".png")
        plot_flow_field(file_name, output_name, args.arrow_step, args.density, value_min, value_max)
        image_names.append(output_name)

    if not args.no_video:
        create_video(image_names, video_name, max(1, args.fps))
    print("Finished: %d image(s) generated in %s" % (len(image_names), output_dir))


if __name__ == "__main__":
    main()
