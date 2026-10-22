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
extern LinkData (*Link);

extern int NowStep, LinkNum;
extern double ParMass, ParDen, ParMom, ParRad, ParTor, ParFai, ParMove, ParForce;		
extern double Sx, Sy, Fx, Fy, Vx, Vy, Va, Ax, Ay, Aa;

//*************************************************************************************************
double GetDelta( int i, int j, int f ) 
{
	double aa = Sq(Ex[f])+Sq(Ey[f]);						
	double bb = ( (D(i)-Sx)*Ex[f] + (D(j)-Sy)*Ey[f] ) * 2;	
	double cc = ( Sq(D(i)-Sx) + Sq(D(j)-Sy) ) - Sq(ParRad);			
	double dd = bb*bb-4*aa*cc;  	{ if( dd < 0 ) dd=0; } 	
	double Delta = (-bb-sqrt(dd))/aa/2; 

	if( Delta>1 && Delta<1+1E-10 )	   Delta = 1;
	if( Delta<0 || Delta>1 )
	{
		cout<<NowStep<<"  delta is error:  "<<setprecision(16)<<Delta<<"   "<<i<<"   "<<j<<"   "<<f<<endl;
	}

	return Delta;
}

//*************************************************************************************************
void ParticleCollide()
{
	#pragma omp parallel for
	for( int n=0; n<LinkNum; ++n )
	{
		double & d = Link[n].d;					
		int & f = Link[n].f;  int r = Re[f];		
		int & i = Link[n].i,  & j = Link[n].j;

		int i0 = i +Ex[f], j0 = j +Ey[f];
		int i1 = i -Ex[f], j1 = j -Ey[f];
		int i2 = i1-Ex[f], j2 = j1-Ey[f];

		double Wx = Vx + Va*(Sy-j-Ey[f]*d);	
		double Wy = Vy - Va*(Sx-i-Ex[f]*d);

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
	Fx = Fy = ParTor = 0;

	for( int n=0; n<LinkNum; ++n )
	{			
		double & d = Link[n].d;					
		int &  f = Link[n].f;  int r = Re[f];		
		int &  i = Link[n].i,  & j = Link[n].j;
		int i1 = i+Ex[f],  j1 = j+Ey[f];

		double Lx = D(i)+d*Ex[f]-Sx,  Ly = D(j)+d*Ey[f]-Sy;
		double fx, fy, Wx = Vx-Va*Ly, Wy = Vy+Va*Lx;  	

		//Galilean invariant momentum exchange method, [Wen et al., J. Comput. Phys., V266, 2014];
		//The hydrodyanamic force is evaluated after the streaming step;
		fx = Dist[i1][j1][f]*(D(Ex[f])-Wx) - Dist[i][j][r]*(D(Ex[r])-Wx);
		fy = Dist[i1][j1][f]*(D(Ey[f])-Wy) - Dist[i][j][r]*(D(Ey[r])-Wy);

		Fx += fx;   Fy += fy;   ParTor += Lx*fy-Ly*fx;
	}	

	Ax = Fx / ParMass;   
	Ay = Fy / ParMass;
	Aa = ParTor/ParMom;
	
	Vx += Ax;   Vy += Ay;   Va     += Aa;
	Sx += Vx;   Sy += Vy;   ParFai += Va;

	ParFai = (ParFai>=PI*2) ? ParFai-PI*2 : (ParFai<-PI*2) ? ParFai+PI*2 : ParFai;
}

//*************************************************************************************************
void ParticleGeometry()
{
	int Hold, Born, Gone, Line;
	LinkNum = Hold = Born = Gone = Line = 0;

	//Check the newborn fluids and the gone fluids
	for( int i=int(Sx-ParRad-0.01); i<=int(Sx+ParRad+1.01); ++i )
	for( int j=int(Sy-ParRad-0.01); j<=int(Sy+ParRad+1.01); ++j )
	{
		GridData & t = Grid[i][j];

		if( Sq(Sx-i)+Sq(Sy-j) <= Sq(ParRad) )   //the nodes inside the particle;
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
		else												
		{
			if( t.Type == PARTICLE )					
			{
				++Born;							//newborn fluid node, outside the particle;
				t.Type = MARKED;

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
	}

	//compute the delta values of the particle;
	for( int i=int(Sx-ParRad-0.01); i<=int(Sx+ParRad+1.01); ++i )
	for( int j=int(Sy-ParRad-0.01); j<=int(Sy+ParRad+1.01); ++j )
	{
		if( Grid[i][j].Type == MARKED )   Grid[i][j].Type = FLUID;	//new fluid node, parallel free;

		if( Grid[i][j].Type == FLUID )
		{
			for( int f=1; f<DQ; ++f )			
			{				
				if( Grid[i+Ex[f]][j+Ey[f]].Type == PARTICLE )		//find the fluid-solid links;
				{
					++Line;
					Link[LinkNum++] = LinkData(i,j,f,GetDelta(i,j,f));
				}				
			} 
		}
	}
}

