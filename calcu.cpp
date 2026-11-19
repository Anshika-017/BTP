// this program is for estimation of optimal fitting parametres wallmu, tau and ss
// xp and yp value of tracer partivle is in r and theta form

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <fstream>
#include <iostream>
using namespace std;
#include <stdlib.h>
#include "reconF.h"

#define PI 3.14159265

//long double reconF(long double xp,long double yp,long double zp,long double ss,long double wallmu,long double tau);
long double *ranF(long double idum); //It generates random numbers and returns a pointer (memory address) to access them.


long double ir1[97]; //array storing 97 values 
int iff1; //it is integer iff1 = 0 → Random-number generator has not been initialized. 
// iff1 = 1 → Random-number generator has been initialized.

int main()
{
	ofstream fwrite1,fwrite2,fwrite3,fwrite4; //Used to write data into files.
	ifstream fread1,fread2; //Used to read data from files.

	long double wallmu,ss,tau,wallmumin,wallmumax,taumin,taumax,ssmax,ssmin,e=1.0,M[108],C,Ms,W,R,theta;
	long double r1,r2,r3,error,phi,idum=-1,xp,yp,zp,pos[1][3],countp,t,sum=0;
	long double *ran1,*ran2,*ran3;
	int i,Q=1,j;
	
	wallmumin=0.0986856;
	wallmumax=0.0986856;
	taumin=1.0E-06;
	taumax=6.0E-06;
	ssmin=0.00050;
	ssmax=0.00060;

	
	fwrite1.open("ss1_det1.txt");
	fwrite2.open("wallmu1_det1.txt");
	fwrite3.open("tau1_det1.txt");
	fwrite4.open("error1_det1.txt");
	fread2.open("det1.txt");			// value of counts recorded by respective detector 

	for(i=0;i<108;i++)					// i is number of calibration point 
	{

		// simply open the det1.txt 
		fread2>>M[i]; // put values in array M
		sum=sum+M[i]; //calculate sum of all counts
		
		//printf (" wallmu %f\n",sum);

	}
	

	Ms=sqrt(sum); //. It calculates a weight for each position
	

	fread2.close();



	long double best=1e30,bss=0,btau=0,bwall=0; int bq=0; int MAXQ=getenv("MAXQ")?atoi(getenv("MAXQ")):2000;
	while (e>0.006 && Q<=MAXQ)
	{
		ran1=ranF(idum); //generate a random function
		r1=ran1[0];  // extract the randomly generated number The random number is intended to lie between 0 and 1.
		r1=ran1[0];  
		idum=ran1[1]; //The second element contains the updated random generator state.Because the next time we call ranF(), we want it to continue generating random numbers from its updated state rather than restart from the beginning.
		
		//printf (" r1 %f\n",r1);

		wallmu=wallmumin+(r1*(wallmumax-wallmumin)); // This is a formula to generate a random parameter between its minimum and maximum limits.
		// wallmu is actually fixed in this program. It is not changing randomly because its minimum and maximum values are identical.
		
		
		
		
		//generate another number 
		ran2=ranF(idum);
		r2=ran2[0];
		idum=ran2[1];
		
		//printf (" wallmu %f\n",wallmu);

		ss=ssmin+(r2*(ssmax-ssmin));




		// repeate for another number

		ran3=ranF(idum);
		r3=ran3[0];
		idum=ran3[1];
		
		//printf (" ss %f\n",ss);

		tau=taumin+(r3*(taumax-taumin));
		//printf (" tau %f\n",tau);



		printf (" no of trail %d\n",Q);

		error=0;
	//	fread2.open("det2.txt");
		fread1.open("gridp.txt");
		
	    for(i=0;i<108;i++)
		{
			for (j=0;j<3;j++)
			{
				//This creates an array with 1 row and 3 columns.
				fread1>>pos[0][j];
			}

			R=pos[0][0];
			theta=pos[0][1];

			xp=R*cos(theta*PI/180);
			yp=R*sin(theta*PI/180);

			zp=pos[0][2];
		    //printf("XP %f\n",i);
			//printf("ZP %f\n",zp);
			countp=reconF(xp,yp,zp,ss,wallmu,tau);
			C=countp;
		
		//	fread2>>M;
		  //printf("MES %f\n",M);
		  //printf("CES %f\n",C);

			W=sqrt(M[i])/Ms;  // this gives us the weight of each calibration point  if measured count is more its weight is more 

			t=(C-M[i])/((C+M[i])/2); // This tells us the relative difference between the predicted and measured values.
			
			phi=W*(pow(t,2));  // phi is the weighted error contribution from the current calibration point.
			
			error=error+phi;
			//printf("%f\n",M[i]);
		}
		fwrite1 << ss << endl;
		fwrite2 << wallmu << endl;
		fwrite3 << tau << endl;
		fwrite4 << error <<endl;

		fread1.close();
		//fread2.close();


		//printf("%f\n",error);

		e=error;
		if(error<best){best=error;bss=ss;btau=tau;bwall=wallmu;bq=Q;printf("  new best error %.6Lf  (trial %d) ss=%.9Lg tau=%.4Lg\n",best,Q,ss,tau);fflush(stdout);}
		Q=Q+1;
	}
	fwrite1.close();
	fwrite2.close();
	fwrite3.close();
	fwrite4.close();
	printf("\nFINISHED after %d trials. last error=%.6Lf  threshold=0.006\n",Q-1,e);
	printf("BEST: error=%.6Lf at trial %d  ss=%.9Lg  tau=%.6Lg  wallmu=%.7Lg\n",best,bq,bss,btau,bwall);
	//fread1.close();
	//fread2.close();
}

long double *ranF(long double idum)
{
	long double m,ia,ic,rm,k,rtemp; static long double iy=0;
	static long double ran1[3];
	int j;
	m=714025;
    ia=1366;
    ic=150889;
    rm=1/m;
   // iff=0;

	if((idum < 0) || (iff1==0))
	//if (idum < 0)
	{
		iff1=1;
		idum=fabs(fmod((ic-idum),m));

		for (j=0;j<97;j++)
		{
			idum=fmod(ia*idum+ic,m);
			ir1[j]=idum;
		}

		idum=fmod(ia*idum+ic,m);

		iy=idum;
		//printf("idum %f\n",idum);

	}
	k=1+(96*iy)/m;
	j=floor(k);

	if((j > 96)||(j < 0))
	{
		printf ("I Am In % d\n",j);
	}

	iy=ir1[j];

	rtemp=iy*rm;
		 
	//printf("rtemp %f\n",rtemp);
	//printf("iy %f\n",iy);

	idum=fmod(ia*idum+ic,m);
	
	//printf("idum1 %f\n",idum);

	ran1[0]=rtemp;
	ran1[1]=idum;
	ran1[2]=iff;

	ir1[j]=idum;
	 //printf("rtemp1 %f\n",rtemp);
	 
	return (ran1);
	 

}