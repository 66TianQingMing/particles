#include <fstream> 
#include <iomanip>
#include <iostream> 
#include <windows.h>
#include "common.h"
using namespace std;

//*************************************************************************************************

extern double	(*Dist)[DY][DQ];
extern double	(*Temp)[DY][DQ];
extern GridData (*Grid)[DY]; 
extern ParticleData Particle[ParticleNum];
extern LinkData (*Link);

extern int NowStep, LinkNum;

//*************************************************************************************************
bool InsideParticle( int p, double x, double y )
{
	ParticleData & t = Particle[p];
	double CosF = cos(t.ParFai), SinF = sin(t.ParFai);
	double X =  CosF*(x-t.Sx) + SinF*(y-t.Sy);
	double Y = -SinF*(x-t.Sx) + CosF*(y-t.Sy);
	return Sq(X/t.ParA)+Sq(Y/t.ParB) <= 1;
}

//*************************************************************************************************
double GetDelta( int p, int i, int j, int f )
{
	ParticleData & t = Particle[p];
	double CosF = cos(t.ParFai), SinF = sin(t.ParFai);
	double X0 =  CosF*(D(i)-t.Sx) + SinF*(D(j)-t.Sy);
	double Y0 = -SinF*(D(i)-t.Sx) + CosF*(D(j)-t.Sy);
	double Xe =  CosF*Ex[f] + SinF*Ey[f];
	double Ye = -SinF*Ex[f] + CosF*Ey[f];
	double aa = Sq(Xe/t.ParA)+Sq(Ye/t.ParB);
	double bb = 2*(X0*Xe/Sq(t.ParA)+Y0*Ye/Sq(t.ParB));
	double cc = Sq(X0/t.ParA)+Sq(Y0/t.ParB)-1;
	double dd = bb*bb-4*aa*cc;  	{ if( dd < 0 ) dd=0; } 	
	double Delta = (-bb-sqrt(dd))/aa/2; 

	if( Delta<0 && Delta>-1E-10 )       Delta = 0;
	if( Delta>1 && Delta<1+1E-10 )	   Delta = 1;
	if( Delta<0 || Delta>1 )
	{
		cout<<NowStep<<"  delta is error:  "<<setprecision(16)<<Delta<<"   "<<p<<"   "<<i<<"   "<<j<<"   "<<f<<endl;
	}

	return Delta;
}

//*************************************************************************************************
void ParticleCollide()
{
	#pragma omp parallel for
	for( int n=0; n<LinkNum; ++n )
	{
		int & p = Link[n].p;
		ParticleData & t = Particle[p];
		double & d = Link[n].d;					
		int & f = Link[n].f;  int r = Re[f];		
		int & i = Link[n].i,  & j = Link[n].j;

		int i0 = i +Ex[f], j0 = j +Ey[f];
		int i1 = i -Ex[f], j1 = j -Ey[f];
		int i2 = i1-Ex[f], j2 = j1-Ey[f];

		double Wx = t.Vx + t.Va*(t.Sy-j-Ey[f]*d);
		double Wy = t.Vy - t.Va*(t.Sx-i-Ex[f]*d);

		if( d < 0.5 )   //Quadratic interpolation boundary condition;
		{
			Dist[i][j][r]  = d*(1.0+d*2)*Dist[i0][j0][f] + (1.0-d*d*4)*Dist[i][j][f] - d*(1.0-d*2)*Dist[i1][j1][f];
			Dist[i][j][r] -= 6.0*Alpha[f]*(Ex[f]*Wx+Ey[f]*Wy); 
		}
		else
		{
			Dist[i][j][r]  = 1.0/d/(d*2+1)*Dist[i0][j0][f] + (d*2-1)/d*Dist[i1][j1][r] - (d*2-1)/(d*2+1)*Dist[i2][j2][r];
			Dist[i][j][r] -= 6.0/d/(d*2+1)*Alpha[f]*(Ex[f]*Wx+Ey[f]*Wy);
		}

		if( Dist[i][j][r] <= 0 )
		{
			cout<<"boundary condition error:  "<<NowStep<<"   "<<i<<"  "<<j<<"  "<<f<<"  :  "<<Dist[i][j][r]<<endl;
		}
	}
} 

//*************************************************************************************************
void ParticleForce()
{	
	for( int p=0; p<ParticleNum; ++p )
	{
		Particle[p].Fx = Particle[p].Fy = Particle[p].ParTor = 0;
	}

	for( int n=0; n<LinkNum; ++n )
	{			
		int & p = Link[n].p;
		ParticleData & t = Particle[p];
		double & d = Link[n].d;					
		int &  f = Link[n].f;  int r = Re[f];		
		int &  i = Link[n].i,  & j = Link[n].j;
		int i1 = i+Ex[f],  j1 = j+Ey[f];

		double Lx = D(i)+d*Ex[f]-t.Sx,  Ly = D(j)+d*Ey[f]-t.Sy;
		double fx, fy, Wx = t.Vx-t.Va*Ly, Wy = t.Vy+t.Va*Lx;

		//Galilean invariant momentum exchange method, [Wen et al., J. Comput. Phys., V266, 2014];
		//The hydrodyanamic force is evaluated after the streaming step;
		fx = Dist[i1][j1][f]*(D(Ex[f])-Wx) - Dist[i][j][r]*(D(Ex[r])-Wx);
		fy = Dist[i1][j1][f]*(D(Ey[f])-Wy) - Dist[i][j][r]*(D(Ey[r])-Wy);

		t.Fx += fx;   t.Fy += fy;   t.ParTor += Lx*fy-Ly*fx;
	}	

	for( int p=0; p<ParticleNum; ++p )
	{
		ParticleData & t = Particle[p];
		t.Ax = t.Fx / t.ParMass;
		t.Ay = t.Fy / t.ParMass;
		t.Aa = t.ParTor/t.ParMom;

		t.Vx += t.Ax;   t.Vy += t.Ay;   t.Va     += t.Aa;
		t.Sx += t.Vx;   t.Sy += t.Vy;   t.ParFai += t.Va;

		t.ParFai = (t.ParFai>=PI*2) ? t.ParFai-PI*2 : (t.ParFai<-PI*2) ? t.ParFai+PI*2 : t.ParFai;
	}
}

//*************************************************************************************************
void ParticleGeometry()
{
	int Hold, Born, Gone, Line;
	LinkNum = Hold = Born = Gone = Line = 0;

	//Mark the old particle nodes;
	FOR_iDX_jDY
	{
		if( Grid[i][j].Type == PARTICLE )   Grid[i][j].Type = MARKED;
	}

	//Check the nodes inside all particles;
	for( int p=0; p<ParticleNum; ++p )
	{
		ParticleData & q = Particle[p];
		double CosF = cos(q.ParFai), SinF = sin(q.ParFai);
		double HalfX = sqrt(Sq(q.ParA*CosF)+Sq(q.ParB*SinF));
		double HalfY = sqrt(Sq(q.ParA*SinF)+Sq(q.ParB*CosF));

		for( int i=int(q.Sx-HalfX-0.01); i<=int(q.Sx+HalfX+1.01); ++i )
		for( int j=int(q.Sy-HalfY-0.01); j<=int(q.Sy+HalfY+1.01); ++j )
		{
			GridData & t = Grid[i][j];
			if( InsideParticle(p,D(i),D(j)) )
			{
				++Hold;
				if( t.Type != PARTICLE )
				{
					++Gone;
					t.Type = PARTICLE;
					t.Den  = t.Vx = t.Vy = 0;
					for( int f=0; f<DQ; ++f )   Dist[i][j][f] = Temp[i][j][f] = 0;
				}
			}
		}
	}

	//Reconstruct the newborn fluid nodes;
	FOR_iDX_jDY
	{
		GridData & t = Grid[i][j];
		if( t.Type == MARKED )
		{
			++Born;							//newborn fluid node, outside the particles;
			int Num = 0;
			for( int f=0; f<DQ; ++f )   Dist[i][j][f] = 0;

			for( int f=1; f<DQ; ++f )		//find the fluid nodes from neighbors
			{
				int i1=i+Ex[f], j1=j+Ey[f], i2=i+Ex[f]*2, j2=j+Ey[f]*2, i3=i+Ex[f]*3, j3=j+Ey[f]*3;

				if( Grid[i1][j1].Type == FLUID )
				{
					++Num;
					if( Grid[i2][j2].Type == FLUID )
					{
						if( Grid[i3][j3].Type == FLUID )
							for( int h=0; h<DQ; ++h )   Dist[i][j][h] += Dist[i1][j1][h]*3 - Dist[i2][j2][h]*3 + Dist[i3][j3][h];
						else
							for( int h=0; h<DQ; ++h )   Dist[i][j][h] += Dist[i1][j1][h]*2 - Dist[i2][j2][h];
					}
					else
					{
						for( int h=0; h<DQ; ++h )   Dist[i][j][h] += Dist[i1][j1][h];
					}
				}
			}

			//calculate the macroscopic quantities;
			t.Den = t.Vx = t.Vy = 0;
			for( int f=0; f<DQ; ++f )
			{
				Dist[i][j][f] = Temp[i][j][f] = Dist[i][j][f] / Num;
				t.Den += Dist[i][j][f];
				t.Vx  += Dist[i][j][f] * Ex[f];
				t.Vy  += Dist[i][j][f] * Ey[f];
			}

			t.Vx /= t.Den;   t.Vy /= t.Den;
			if( fabs(t.Den-1) > 0.1 )   cout<<"new fluid node error:  "<<i<<"   "<<j<<"   "<<t.Den<<endl;
		}
	}

	FOR_iDX_jDY
	{
		if( Grid[i][j].Type == MARKED )   Grid[i][j].Type = FLUID;
	}

	//Compute the delta values of all particles;
	for( int p=0; p<ParticleNum; ++p )
	{
		ParticleData & q = Particle[p];
		double CosF = cos(q.ParFai), SinF = sin(q.ParFai);
		double HalfX = sqrt(Sq(q.ParA*CosF)+Sq(q.ParB*SinF));
		double HalfY = sqrt(Sq(q.ParA*SinF)+Sq(q.ParB*CosF));

		for( int i=int(q.Sx-HalfX-0.01); i<=int(q.Sx+HalfX+1.01); ++i )
		for( int j=int(q.Sy-HalfY-0.01); j<=int(q.Sy+HalfY+1.01); ++j )
		{
			if( Grid[i][j].Type == FLUID )
			{
				for( int f=1; f<DQ; ++f )
				{
					int i1=i+Ex[f], j1=j+Ey[f];
					if( Grid[i1][j1].Type == PARTICLE && InsideParticle(p,D(i1),D(j1)) )
					{
						if( LinkNum >= LinkMax )
						{
							cout<<"Link number is larger than LinkMax:  "<<LinkNum<<endl;
							return;
						}
						++Line;
						Link[LinkNum++] = LinkData(p,i,j,f,GetDelta(p,i,j,f));
					}
				}
			}
		}
	}
}

