# Project coding conventions

- Preserve the existing style of `Migration.cpp`, `Particle.cpp`, and `common.h` when modifying code.
- Prefer procedural, C-style code and plain functions, arrays, structs, and loops.
- Keep the current naming, indentation, brace placement, and file organization unless a change requires otherwise.
- Avoid introducing classes, templates, STL containers, lambdas, smart pointers, or other modern C++ abstractions unless they are clearly necessary.
- Keep compatibility with the project's existing Visual Studio 2012 (`v110`) toolset.
- Make focused changes and avoid unrelated formatting rewrites.
- After each completed change, check the diff, create a Git commit, and push `main` to the configured GitHub remote.
