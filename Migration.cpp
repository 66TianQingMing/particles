#include <conio.h>
#include <cstdlib>
#include <direct.h>
#include <errno.h>
#include <fstream> 
#include <iomanip>
#include <iostream> 
#include <windows.h>
#include "common.h"
using namespace std;

//*************************************************************************************************
double	  (*Dist)[DY][DQ];
double	  (*Temp)[DY][DQ];
GridData  (*Grid)[DY];
ParticleData Particle[ParticleNum];
LinkData  (*Link);
const char * BcName[] = { "HBB", "BFL", "LL2" };

ofstream TrackFile;
double ExactUx[DY], ExactUxMod, NormalErr1, NormalErr2;
double Mass, Density, Tau, ReTau, Viscosity, Pressure, DenLeft, DenRight, Radius, Reynolds;		
int    TaskNum=1, No, NowStep, AllStep, ShowStep, SaveStep, BeginTime, StepTime, LastTime, LinkNum, BcType;
int    OutputType=TRACK_OUTPUT;

double ParMove;

//*************************************************************************************************
void Initialize();
void SetFlowField();
void ShowData();
void SaveOutput();
void SaveData();
void SaveFlowField();
void SaveTask();

void GlobalCollide();
void BoundaryCollide();
void InletAndOutlet();
void MacroCalculate();	
void ParticleCollide();
void ParticleForce();	
void ParticleGeometry();							
void LatticeRemesh();									
int  TimeInterval()  { StepTime=GetTickCount()-LastTime;  LastTime=GetTickCount();  return StepTime; }


//*************************************************************************************************
int main( int argc,char *argv[] )
{
	for( No=0; No<TaskNum; ++No )
	{
		Initialize();
		SetFlowField();
		ShowData( );
		SaveOutput( );

		for( NowStep=1; NowStep<=AllStep; ++NowStep )	
		{
			GlobalCollide();			
			BoundaryCollide();	
			InletAndOutlet();

#ifdef MIGRATION
			ParticleCollide();
			ParticleForce();	
			ParticleGeometry();							
			LatticeRemesh();									
#endif
			MacroCalculate();

			if( NowStep%SaveStep == 0 )  { SaveOutput(); }
			if( NowStep%ShowStep == 0 )  { TimeInterval();  ShowData(); }
		}

		SaveTask();
	}

	TrackFile.close();
	delete [] Grid;   delete [] Link;
	delete [] Dist;   delete [] Temp;

	cout<<endl<<" Press Enter key to quit ...... ";   cin.get();
	return 0;
}


//*************************************************************************************************
void Initialize()
{
	if( _mkdir("result")!=0 && errno!=EEXIST )
	{
		cout<<"Create result directory error ..."<<endl<<flush;
		cin.get();   exit(1);
	}

	Density	  = 1;
	Tau		  = 0.75;
	BcType    = LL2; 
	Pressure  = 8E-4;
		
	AllStep	  = 1000 * 10000;
	SaveStep  = 1000;
	ShowStep  = 1000;
	NowStep   = StepTime = 0;
	BeginTime = LastTime = GetTickCount();

#ifdef MIGRATION
	ParMove = 0;
	for( int p=0; p<ParticleNum; ++p )
	{
		ParticleData & t = Particle[p];
		t.ParA      = 12.5;
		t.ParB      = 6.25;
		t.ParDen    = Density;
		t.ParMass   = PI*t.ParA*t.ParB*t.ParDen;
		t.ParMom    = t.ParMass*(Sq(t.ParA)+Sq(t.ParB))/4;
		t.Sx        = DX/2+(D(p)-D(ParticleNum-1)/2)*t.ParA*4;
		t.Sy        = D(p+1)*D(DY-1)/D(ParticleNum+1);
		t.ParTor = t.ParFai = t.Fx = t.Fy = t.Vx = t.Vy = t.Va = t.Ax = t.Ay = t.Aa = 0;
	}
#endif

	if( No == 0 )
	{
		Dist = new double [DX][DY][DQ];
		Temp = new double [DX][DY][DQ]; 
		Grid = new GridData [DX][DY];
		Link = new LinkData [LinkMax];

		if( Dist==0 || Temp==0 || Grid==0 || Link==0 )
		{
			cout<<"Memory allocation is error ..."<<endl<<flush;
			cin.get();   exit(1);
		}
	}

	fill( (double*)Dist, (double*)Dist+DX*DY*DQ, 0 );		
	fill( (double*)Temp, (double*)Temp+DX*DY*DQ, 0 );		
	memset( (void*)Grid,  0, DX*DY*sizeof(GridData));
}


//*************************************************************************************************
void SetFlowField()		
{
	Radius    = D(DY-1)/2;
	ReTau     = double(1)/Tau;
	Viscosity = (Tau*2-1)/6;
	DenLeft   = Density + Pressure/2*3;
	DenRight  = Density - Pressure/2*3;
	Reynolds  = (Pressure/(DX-1)*Sq(Radius)/Viscosity/3) * (DY-1)/Viscosity;

	FOR_iDX_jDY
	{
		GridData & t = Grid[i][j];
		if( j==0 || j==DY-1 )
		{
			t.Type = BOUNDARY;		t.Den = 0;
		}
		else
		{
			t.Type = FLUID;			t.Den = Density;			
		}
		
		t.Vx  = t.Vy = 0;
		for( int f=0; f<DQ; ++f )	
		{
			Dist[i][j][f] = Temp[i][j][f] = (t.Type==BOUNDARY) ? UNKNOWN : Feq(f, t.Den, t.Vx, t.Vy);
		}
	} 

#ifdef MIGRATION
	ParticleGeometry();
#endif

	Mass = 0;
	FOR_iDX_jDY_Fluid   Mass += Grid[i][j].Den;

	//Compute the analytical resolution of Poiseuille flow;
	ExactUxMod = 0;
	ExactUx[0] = ExactUx[DY-1] = 0;
	for( int j=1; j<DY-1; ++j )
	{
		double r = D(j) - D(DY-1)/2;
		ExactUx[j] = Pressure/(DX-1)/Viscosity*(Sq(Radius)-Sq(r))/2;	//fluid mechanics,P297;
		ExactUxMod += Sq(ExactUx[j]);
	}
	ExactUxMod = sqrt(ExactUxMod);

	cout<<endl<<" Migration: DX="<<DX<<"   DY="<<DY<<"   Tau="<<Tau<<"   BC="<<BcName[BcType]<<"   Pressure="<<Pressure<<"   Reynolds="<<Reynolds<<endl;
	cout<<"**************************************************************************************"<<endl;
}

//*************************************************************************************************
void GlobalCollide()		
{
#pragma omp parallel for 
	FOR_iDX_jDY_Fluid
	{
		for( int f=0; f<DQ; ++f )
		{	
			int ii = i+Ex[f];	if(ii < 0) ii+=DX;	else if(ii >= DX) ii-=DX; 
			int jj = j+Ey[f];	if(jj < 0) jj+=DY;	else if(jj >= DY) jj-=DY; 
			Temp[ii][jj][f] = Dist[i][j][f] - ReTau*(Dist[i][j][f] - Feq(f, Grid[i][j].Den, Grid[i][j].Vx, Grid[i][j].Vy));
		}
	}

	double (*p)[DY][DQ];
	p = Dist;   Dist = Temp;   Temp = p;
}

//*************************************************************************************************
void BoundaryCollide()	
{
	const double d = 1.0;
	const int fb[6][2] = { {1,4}, {1,7}, {1,8}, {DY-2,2}, {DY-2,5}, {DY-2,6} };

	#pragma omp parallel for
	for( int i=0; i<DX; ++i )		    
	{
		for( int n=0; n<6; ++n )
		{
			int TheBcType = BcType;
			int j = fb[n][0], f = fb[n][1];   int r = Re[f];   

			int i0 = i +Ex[f], j0 = j +Ey[f];
			int i1 = i -Ex[f], j1 = j -Ey[f];
			int i2 = i1-Ex[f], j2 = j1-Ey[f];
			if(i0 < 0) i0+=DX;	else if(i0 >= DX) i0-=DX; 

			if( TheBcType==LL2 && (i2<0 || i2>=DX) )   TheBcType = BFL;	 
			if( TheBcType==BFL && (i1<0 || i1>=DX) )   TheBcType = HBB;	 

			switch( TheBcType )
			{
			case HBB:	//Halfway bounce back;
				Dist[i][j][r] = Dist[i0][j0][f]; 
				break;
				 
			case BFL:   //Bouzidi interpolation boundary condition;
				if( d < 0.5 )
					Dist[i][j][r]  = d*2*Dist[i0][j0][f] + (1.0-d*2)*Dist[i][j][f];
				else
					Dist[i][j][r]  = 1.0/(d*2)*Dist[i0][j0][f] + (2*d-1)/(d*2)*Dist[i1][j1][r];
				break;

			case LL2:	//Quadratic interpolation boundary condition;
				if( d < 0.5 )
					Dist[i][j][r]  = d*(1.0+d*2)*Dist[i0][j0][f] + (1.0-d*d*4)*Dist[i][j][f] - d*(1.0-d*2)*Dist[i1][j1][f];
				else
					Dist[i][j][r]  = 1.0/d/(d*2+1)*Dist[i0][j0][f] + (d*2-1)/d*Dist[i1][j1][r] - (d*2-1)/(d*2+1)*Dist[i2][j2][r];
				break;
			} 

			if( Dist[i][j][r]<=0 || Dist[i][j][r]>=1 )
			{
				cout<<"Boundary condition error:  ("<<i<<","<<j<<","<<f<<")  :  "<<Dist[i][j][r]<<endl;
			}			
		}
	} 
} 

//*************************************************************************************************
void MacroCalculate(void)
{
	double TmpMass = 0;

#pragma omp parallel for reduction(+:TmpMass) 
	FOR_iDX_jDY_Fluid
	{
		GridData & t = Grid[i][j];
		double *  Df = Dist[i][j];

		t.Den = Df[0]+Df[1]+Df[2]+Df[3]+Df[4]+Df[5]+Df[6]+Df[7]+Df[8];
		t.Vx  = (Df[1]+Df[5]+Df[8]-Df[3]-Df[6]-Df[7]) / t.Den;
		t.Vy  = (Df[2]+Df[5]+Df[6]-Df[4]-Df[7]-Df[8]) / t.Den;

		if( t.Den!=t.Den || fabs(t.Den-Density)>0.1 )
		{
			cout<<"Density: "<<NowStep<<",  ("<<i<<","<<j<<")   Den-1.0 = "<<(t.Den-1)<<endl;
		}

		TmpMass += t.Den;
	}

	Mass = TmpMass;
}

//*************************************************************************************************
void InletAndOutlet()
{
	//Pressure boundary condition;	
	const double W1 = 1.0/9;
	for( int j=1; j<DY-1; ++j )   
	{
		double * t = Dist[0][j];
		double UxLeft  = 1.0 -(t[0]+t[2]+t[4]+2.0*(t[3]+t[6]+t[7]))/DenLeft; 
		t[1] = t[3]+6*DenLeft*W1*UxLeft; 
		t[5] = t[7]-0.5*(t[2]-t[4])+(0.5-3*W1)*DenLeft*UxLeft;        
		t[8] = t[6]+0.5*(t[2]-t[4])+(0.5-3*W1)*DenLeft*UxLeft;        

		t = Dist[DX-1][j];
		double UxRight = -1.0+(t[0]+t[2]+t[4]+2.0*(t[1]+t[5]+t[8]))/DenRight;
		t[3] = t[1]-6*DenRight*W1*UxRight;  
		t[6] = t[8]-0.5*(t[2]-t[4])-(0.5-3*W1)*DenRight*UxRight;   
		t[7] = t[5]+0.5*(t[2]-t[4])-(0.5-3*W1)*DenRight*UxRight;  
	}
}

//*************************************************************************************************
void LatticeRemesh()
{
	const int Span = 2;
	double CenterX = 0;
	for( int p=0; p<ParticleNum; ++p )   CenterX += Particle[p].Sx;
	CenterX /= ParticleNum;

	if( CenterX-DX/2 > Span )
	{
		for( int p=0; p<ParticleNum; ++p )   Particle[p].Sx -= Span;
		ParMove += Span;
		for( int n=0; n<LinkNum; ++n )  Link[n].i -= Span;

		for( int i=5; i<DX-5-Span; ++i )
		{
			memcpy( (void*)&Grid[i], (const void*)&Grid[i+Span], DY*sizeof(GridData ) );
			memcpy( (void*)&Dist[i], (const void*)&Dist[i+Span], DY*DQ*sizeof(double) );
		}

		for( int i=DX-5-Span; i<DX-5-1; ++i )	//smooth the tail part;
		{
			memcpy( (void*)&Grid[i], (const void*)&Grid[DX-5-1], DY*sizeof(GridData) );
			memcpy( (void*)&Dist[i], (const void*)&Dist[DX-5-1], DY*DQ*sizeof(double) );
		}
	}
}

//*************************************************************************************************
void ShowData( )		
{
	//Normal error on the cross section;
	NormalErr1 = NormalErr2 = 0;
	for( int j=1; j<DY-1; ++j )
	{
		NormalErr1 += Sq(Grid[10   ][j].Vx-ExactUx[j]);
		NormalErr2 += Sq(Grid[DX-11][j].Vx-ExactUx[j]);
	}
	NormalErr1 = sqrt(NormalErr1)/ExactUxMod;
	NormalErr2 = sqrt(NormalErr2)/ExactUxMod;	

#ifdef MIGRATION
	cout<<"   "<<NowStep<<"   "<<setiosflags(ios::fixed)<<setprecision(10)<<Mass<<setprecision(3);
	for( int p=0; p<ParticleNum; ++p )
		cout<<"   P"<<p<<"("<<Particle[p].Sx+ParMove-DX/2<<","<<Particle[p].Sy<<","<<Particle[p].ParFai<<")";
	cout<<setprecision(2)<<"   "<<NormalErr1*100<<"%"<<"   "<<NormalErr2*100<<"%"<<"   "<<StepTime<<endl;
#else
	cout<<"   "<<NowStep<<"     "<<setiosflags(ios::fixed)<<setprecision(10)<<Mass<<setprecision(2)<<"     "<<NormalErr1*100<<"%"<<"     "<<NormalErr2*100<<"%"<<"     "<<StepTime<<endl;
#endif
}

//*************************************************************************************************
void SaveOutput()
{
	if( OutputType == FIELD_OUTPUT )	SaveFlowField();
	else							SaveData();
}

//*************************************************************************************************
void SaveData()			
{
#ifdef MIGRATION
	if( NowStep == 0 )
	{ 
		char FileName[128];
		sprintf( FileName, "result\\Migration_%d_%d_N%d_A%3.2f_B%3.2f.txt", DX-1, DY-1, ParticleNum, Particle[0].ParA, Particle[0].ParB );
		TrackFile.open( FileName );
		TrackFile<<"Step   Particle   Mass   Sx   Sy   Vx   Vy   Va   Fai   Fx   Fy   Tor"<<endl;
	}
	for( int p=0; p<ParticleNum; ++p )
	{
		ParticleData & t = Particle[p];
		TrackFile<<NowStep<<"  "<<p<<"  "<<Mass<<"  "<<(t.Sx+ParMove-DX/2)/100<<"  "<<t.Sy/100<<"  "<<t.Vx<<"  "<<t.Vy<<"  "<<t.Va<<"  "<<t.ParFai<<"  "<<t.Fx<<"  "<<t.Fy<<"  "<<t.ParTor<<endl;
	}
#else
	char FileName[128];
	ofstream File;
	sprintf( FileName, "result\\Velocity_%s_T%3.2f_R%3.2f.txt", BcName[BcType], Tau, Radius*2 );
	File.open( FileName );

	File<<"Scale   ExactUx   Ux   Ux1   Ux2   Density"<<endl;
	File<<"0   0   0   0   0   0"<<endl;
	for( int j=1; j<DY-1; ++j )
	{
		File<<setprecision(12)<<double(j)/Radius/2<<"   "<<ExactUx[j]<<"   "<<Grid[DX/2][j].Vx<<"   "<<Grid[0][j].Vx<<"   "<<Grid[DX-1][j].Vx<<"   "<<Grid[DX/2][j].Den<<endl;
	}
	File<<"1   0   0   0   0   0"<<endl<<endl;
	File.close();
#endif

}

//*************************************************************************************************
void SaveFlowField()
{
	char FileName[128];
	ofstream File;
	sprintf( FileName, "result\\FlowField_%08d.txt", NowStep );
	File.open( FileName );

	File<<setprecision(12);
	File<<"# Step "<<NowStep<<endl;
	File<<"# GridSize "<<DX<<" "<<DY<<endl;
	File<<"# ParticleNum "<<ParticleNum<<endl;
	for( int p=0; p<ParticleNum; ++p )
	{
		ParticleData & t = Particle[p];
		File<<"# Particle "<<p<<" Sx "<<t.Sx+ParMove<<" Sy "<<t.Sy
			<<" A "<<t.ParA<<" B "<<t.ParB<<" Fai "<<t.ParFai<<endl;
	}
	File<<"# X Y Type Density Vx Vy"<<endl;

	FOR_iDX_jDY
	{
		GridData & t = Grid[i][j];
		File<<D(i)+ParMove<<" "<<j<<" "<<t.Type<<" "<<t.Den<<" "<<t.Vx<<" "<<t.Vy<<endl;
	}
	File.close();
}

void SaveTask()
{
}
