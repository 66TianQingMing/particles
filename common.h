/**********************************************************************
*                                                                     *
*           Lattice Boltzmann Method for Particle Migration           *
*                                                                     *
*           Binghai Wen  @ Guangxi Normal University, China           *
*                                                                     *
**********************************************************************/

#pragma once 
#pragma warning( disable : 4018 )  //signed & unsigned;
#pragma warning( disable : 4267 )  //signed & unsigned;
#pragma warning( disable : 4786 )  //for stl;
#pragma warning( disable : 4996 )  //for sprintf;

#define MIGRATION

//*************************************************************************************************
const int DX	= 1001; 
const int DY	= 101;  

//*************************************************************************************************
const char FLUID	= 1;		
const char BOUNDARY	= 2;		
const char PARTICLE	= 3;		
const char MARKED	= 4;		
const char UNKNOWN	= 97;	

//*************************************************************************************************
const int HBB		= 0;		
const int BFL		= 1;
const int LL2		= 2;

//*************************************************************************************************
const int DQ = 9;
const int Ex[DQ] = { 0,  1,  0, -1,  0,  1, -1, -1,  1 };
const int Ey[DQ] = { 0,  0,  1,  0, -1,  1,  1, -1, -1 };	
const int Re[DQ] = { 0,  3,  4,  1,  2,  7,  8,  5,  6 };	
const double Alpha[DQ]={ 4.0/9, 1.0/9, 1.0/9, 1.0/9, 1.0/9, 1.0/36, 1.0/36, 1.0/36, 1.0/36 };
const double PI = 3.14159265358979323846264338327950288;
const double S2 = 1.4142135623730950488016887242097;   //sqrt(2.0);

//*************************************************************************************************
#define D(x)  (double(x))
#define FOR_iDX_jDY	 for(int i=0; i<DX;++i) for(int j=0;j<DY;++j)
#define FOR_iDX_jDY_Fluid  for(int i=0; i<DX;++i) for(int j=0;j<DY;++j) if(Grid[i][j].Type == FLUID)

//Square of a number;
template<class T> inline T Sq( T x )
{
	return x*x;
}

//Determine if two doubles are equal;
inline bool Eq( double x, double y, double d=1E-12 )
{
	return ( x>y ? x-y : y-x ) < d;	
}

//Equilibrium distribution function;
inline double Feq( int f, double Den, double Vx, double Vy ) 
{
	double DotMet = Vx*Ex[f] + Vy*Ey[f];	
	return Den*Alpha[f]*(1.0+3.0*DotMet+4.5*DotMet*DotMet-1.5*(Vx*Vx+Vy*Vy));
}

struct GridData
{
	int Type;
	double Den, Vx, Vy;
};

struct LinkData
{
	double d;
	int i, j, f;
	LinkData() {  };
	LinkData( int i, int j, int f, double d )
	{
		this->i=i;  this->j=j;  this->f=f;  this->d=d;
	}
};