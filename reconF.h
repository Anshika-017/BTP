// this program is for one dimensional movement of tracer particle (Monte Carlo program)

#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
using namespace std;
#include <stdlib.h>
#include <cstdlib>
#include <time.h>


long double montecarlo(long double xp,long double yp,long double zp,long double xd,long double yd,long double zd,long double xdd,long double ydd,int k,int i,long double ss,long double tau,long double wallmu);
long double *hrho(long double xp,long double yp,long double zp,long double xd,long double yd,long double zd,long double xdd,long double ydd,int k);
long double depthtwo(long double tetcri,long double theta,long double alpha,long double rho,long double h,int ihist);
long double depthone(long double tetcri,long double theta,long double alpha,long double rho,long double h,long double oa,long double ob,int ihist);
long double *ran(long double idum);
long double *abscof(long double egamma);
long double photopeak(long double egamma);

//long double recon(int j,int i,long double xp,long double yp,long double zp)
long double ir[97]; // array
int iff; // It is used as a flag in the random-number generator to keep track of whether its initialization has already happened.


long double reconF(long double xp,long double yp,long double zp,long double ss,long double wallmu,long double tau)
{
	//ofstream fwrite;   //ofstream is used to write data into a file.
	ifstream fread1;   // for reading the data


//	long double xd,yd,xdd,ydd,zd,count,countp,xp,yp,zp;
	long double count,countp,xd,yd,xdd,ydd,zd;   // xdd, ydd	Additional detector geometry coordinates
	long double m[1][5]; //array of  row and  columns 
	int l, k=1;
	int i=0;


	fread1.open("crysd1.txt");
	//fread2.open("gridp1.txt");
	//fwrite.open("count1.txt");
	for(l=0;l<5;l++)
	{
		fread1>>m[0][l];
	}
	xd=m[i][0];
	xdd=m[i][1];
	yd=m[i][2];
	ydd=m[i][3];
	zd=m[i][4];
	//printf("detector co-ordinate is-%15.5f%10.1f%10.1f%10.1f%10.1f",xd,yd,zd,xdd,ydd);

	//for (ndim=0;ndim<73;ndim++)
	//{
	//	for (q=0;q<3;q++)
	//	{
	//		fread2>>pos[0][q];
	//	}
	//	xp=pos[0][0];
	//	yp=pos[0][1];
	//	zp=pos[0][2];
		
	//	printf(" doing for point....... %d\n",ndim);
		//printf("zp %f\n",zp);
		count=montecarlo(xp,yp,zp,xd,yd,zd,xdd,ydd,k,i,ss,tau,wallmu);
		//printf("%f",count);
		countp=count;
	//	fwrite<<countp<<endl;
	//}

	//fread1.close();

	return (countp);
	//return 0;

	//fread2.close();
	//fwrite.close();
	
}
// dd-> Distance between the two detector points ; hrho1[2]-> Array storing two calculated values ; 
// h-> Distance along the detector line (projection length) ; rho-> Perpendicular distance from the source point to the detector line
long double *hrho(long double xp,long double yp,long double zp,long double xd,long double yd,long double zd,long double xdd,long double ydd,int k)
{
	static long double hrho1[2]; long double dd,h,rho;
	dd=sqrt(pow((xd-xdd),2)+pow((yd-ydd),2)+pow((zd-zd),2));
	h=fabs((xp-xdd)*(xd-xdd)+(yp-ydd)*(yd-ydd)+(zp-zd)*(zd-zd))/dd;
	rho=sqrt(pow((xp-xdd),2)+pow((yp-ydd),2)+pow((zp-zd),2)-pow(h,2));

//	printf(" i m in hrho");
	hrho1[0]=h;
	hrho1[1]=rho;
	//printf("rho %f\n",rho);
//	printf("h %f\n",h);
//	printf("dd %f\n",dd);


	return hrho1;
}

long double montecarlo(long double xp,long double yp,long double zp,long double xd,long double yd,long double zd,long double xdd,long double ydd,int k,int i,long double ss,long double tau,long double wallmu)
{
	long double tt=5,ray=25.4,perech=0.02,pi=3.1415926,*hrho2,idum=-1;

	long double solang=0,inteff=0,solangsqr=0; // solang = solid angle ; inteff = Integrated efficiency ; solangsqr	Square of solid angle
	
	long double alfmax,walf,alpha,oa,tetmax,tetmin,ob,tetcri,length;
	// alfmax-> Maximum alpha angle ; walf	-? Weight or interval associated with alpha ; alpha	-> alpha angle ; oa	-> Geometry-related quantity ; tetmax-> Maximum theta angle ; tetmin-> Minimum theta angle ; ob->Geometry-related quantity ;tetcri->Critical theta angle ; length->A length or path distance
	
    long double theta,wthet,prob,effint,solang1,term,count,*abscof4;
	// theta-> Angular variable ; wthet-> Weight or interval associated with theta ; prob-> Probability ; effint-> Efficiency-related intermediate value ; solang1-> Intermediate solid-angle value ; term-> Intermediate calculation term ; count-> Calculated count ; abscof4-> Pointer to a long double variable

	long double rho,r,h,munai,d1,d2,ppic,munaip,munaic;
	// rho-> Geometrical distance/quantity ; r-> Distance or radius ; h-> Geometrical length or height ; munai-> Attenuation-related quantity ; d1,d2-> Intermediate distances or lengths ; ppic-> Intermediate calculated value ; munaip,munaic-> Additional attenuation-related quantities


	long double *ran3,*ran4,*ran5,*ran6;  //pointers to random numbers 

	int gperd=2,ihist=12000,p,iff=0; //ihist=12000 means the program is set up to perform 12,000 Monte Carlo histories/samples.
	
	//ofstream fwrite1,fwrite2,fwrite3,fwrite4,fwrite5,fwrite6,fwrite7;
	ifstream fread;

//	fread.open("setpara.txt");
//	for(n=0;n<3;n++)
//	{
//		fread>>para[0][n];
//	}
//	wallmu=para[i][0];
//	ss=para[i][1];
//	tau=para[i][2];
//	wallmu=0.0150255;
//	ss=0.000153703;
//	tau=1.98E-06;
	hrho2=hrho(xp,yp,zp,xd,yd,zd,xdd,ydd,k);
	h=hrho2[0];
	rho=hrho2[1];
	abscof4=abscof(1);

	// munai = munaip + munaic;  // Total = Photoelectric + Compton coefficients
	munai=abscof4[0];
	munaip=abscof4[1];
	munaic=abscof4[2];
	

	ppic=photopeak(1);
	

//	printf("%f\n",rho);
//	printf("rho %f\n",rho);
//	printf("h %f\n",h);

	length=2*ray;

	if (rho > ray) //This condition checks whether the source point lies outside the detector's radius in the perpendicular direction.
	{
		alfmax=asin(ray/rho);
		walf=alfmax/pi;  //This calculates the normalized angular weight associated with the allowed range of angles
		
		//monte carlo iterations beginss now

		for (p=0;p<ihist;p++)
		{
			ran3=ran(idum);
			r=ran3[0];
			idum=ran3[1];
			//printf("r%f\n",r);
			
			//srand ( time(NULL) );
			//r = rand() % 100000  ;
			//r=r/100000;



			//ob = distance from the source to the nearer intersection.
			// oa = distance from the source to the farther intersection.
			
			//theta tells us the inclination of the gamma-ray path relative to the detector's longitudinal axis.
			//alpha is the angle that determines the sideways direction of a ray from the source as it enters the circular detector cross-section.


			alpha=alfmax*(2*r-1);
			oa= rho*cos(alpha)+sqrt(pow(ray,2)-pow((rho*sin(alpha)),2));   // This calculates one of the two intersection distances between the selected ray and the circular detector boundary.

			tetmax=atan(oa/h); /// max theta angle 

			ob= rho*cos(alpha)-sqrt(pow(ray,2)-pow((rho*sin(alpha)),2));
			tetmin= atan(ob/(h+length));
			tetcri=atan(ob/h);

			
	    	ran4=ran(idum);
			r=ran4[0];
			idum=ran4[1];

			//printf("r2 %f\n",r);
			//printf("idum %f\n",idum);
	    	//srand ( time(NULL) );

			//r = rand() % 100000  ;
			//r=r/100000;

			theta= acos(cos(tetmin)-r*(cos(tetmin)-cos(tetmax)));
			wthet=0.5*(cos(tetmin)-cos(tetmax));

			//printf("h %f\n",h);

			d1=depthone(tetcri,theta,alpha,rho,h,oa,ob,ihist);  //Calculates the distance travelled by the ray inside the detector.

			prob=1-(exp(-munai*d1)); // This calculates the probability that the gamma ray undergoes at least one interaction inside the detector.
		
			prob=prob*exp(-wallmu*tt/cos(theta));  // This multiplies the interaction probability by the transmission probability through the wall.
			
			solang=solang+wthet*walf;
			effint=prob*wthet*walf; //weighted interaction probability for the current ray.
			inteff=inteff+effint; //accumulated weighted interaction probability over all simulated rays.
			solangsqr=solangsqr+pow((wthet*walf),2); 
		}
		//printf("d1 %f\n",d1);
		//printf("alpha %f\n",alpha);
		//printf("theta %f\n",theta);
		//printf("alfmax %f\n",alfmax);
		
		

		solang1=solang/ihist;
		inteff=inteff/ihist;

		//ss = source strength (activity), likely expressed in curies (Ci).
		//tau = tau is a parameter that controls the count-rate reduction due to dead-time effects in the detector or counting system.

		term=gperd*((3.7E10)*ss)*ppic*inteff;  //The result, term, represents the estimated count rate before the final correction.
		
		count=perech*term/(1+tau*term);  // perech acts as a scaling factor of 0.02 (2%) applied to the calculated count.

		//printf("count %f\n",count);
		//printf("inteff %f\n",inteff);
		
	}
	else  // rho < = ray
	{
		tetmax=atan((ray+rho)/h);
        tetcri=atan((ray-rho)/h);
        tetmin=0;

		wthet=0.5*(cos(tetmin)-cos(tetmax));
		for (p=0;p<ihist;p++)
		{
			ran5=ran(idum);
			r=ran5[0];
			idum=ran5[1];

			
	    	//srand ( time(NULL) );

			//r = rand() % 100000  ;
			//r=r/100000;

			theta=acos(cos(tetmin)-r*(cos(tetmin)-cos(tetmax)));
			if (theta >= tetcri)
			{
				alfmax=acos((pow(rho,2)+pow((h*tan(theta)),2)-pow(ray,2))/(2*h*rho*tan(theta)));
                walf=alfmax/pi; //walf represents the fraction of the full angular range that is allowed.
			}
			else
			{
				alfmax=pi;
                walf=1;
			}

			ran6=ran(idum);
			r=ran6[0];
			idum=ran6[1];

			
			//srand ( time(NULL) );

			//r = rand() % 100000  ;
			//r=r/100000;

			alpha=alfmax*(2*r-1);

			d2=depthtwo(tetcri,theta,alpha,rho,h,ihist);
			prob=1-exp(-munai*d2);
            prob=prob*exp((-wallmu*tt)/(cos(theta)));

			solang=solang+wthet*walf;
            effint=prob*wthet*walf;
            inteff=inteff+wthet*walf*prob;
            solangsqr=solangsqr+pow((wthet*walf),2);
		}

		solang1=solang/ihist;
		inteff=inteff/ihist;
		term=gperd*((3.7E10)*ss)*ppic*inteff;
		count=perech*term/(1+tau*term);
		//printf("d2 %f\n",d2);
		//printf("alpha %f\n",alpha);
		//printf("theta %f\n",theta);
		//printf("alfmax %f\n",alfmax);
		//printf("count %f\n",count);


	}
	

	
	//printf("%f",count);
	return (count);
}


long double depthtwo(long double tetcri,long double theta,long double alpha,long double rho,long double h,int ihist)
{
	long double ray=25.4,length,pi=3.1415926;
	long double lambda,depth2,d2;

	length=2*ray;

	//x1p=h*tan(theta)*sin(alpha);
	//y1p=h*tan(theta)*cos(alpha);
	//z1p=0;

	lambda=sqrt(pow(((h+length)*tan(theta)),2) + pow(rho,2) - 2*(h+length)*tan(theta)*rho*cos(alpha));

	if (lambda <= ray)
	{
		depth2=length/cos(theta);
	}
	else
	{
		depth2=(-h/cos(theta))+(rho*cos(alpha)+sqrt(pow(ray,2)-pow((rho*sin(alpha)),2)))/sin(theta);
	}

	d2=depth2;
//	printf("d2 %f",d2);

	return (d2);
}


long double depthone(long double tetcri,long double theta,long double alpha,long double rho,long double h,long double oa,long double ob,int ihist)
{
	long double ray=25.4,length,pi=3.1415926;
	long double x1p,y1p,z1p,depth1,d1;

	length=2*ray;
//	printf("h2 %f\n",h);

	if(theta < tetcri)  // if(theta < tetcri)
	// In this case, the ray enters the detector through its curved cylindrical surface.
	{
		
		x1p=ob*sin(alpha);
		y1p=ob*cos(alpha);
		z1p=h-ob/tan(theta);

		if (oa/tan(theta) >= (h+length))  //This checks whether the ray reaches the far end of the detector before exiting through the opposite curved surface.
		{
			depth1=-ob/sin(theta)+(h+length)/cos(theta);
		}
		else
		{
			depth1=(oa-ob)/sin(theta);
		}
	}
	else //In this case, the ray enters through the front circular face of the detector rather than through its curved surface.

	{
		x1p=h*tan(theta)*sin(alpha);
		y1p=h*tan(theta)*cos(alpha);
		z1p=0;

		if (oa/tan(theta) >= (h+length))
		{
			depth1=length/cos(theta);
		}
		else
		{
			depth1=oa/sin(theta)-h/cos(theta);
		}

	}

	d1=depth1;
//	printf("d1 %f",d1);

	return (d1);
}


long double *ran(long double idum)
{
	long double m,ia,ic,rm,k,rtemp; static long double iy=0; static int iff=0;
	static long double ran1[3];
	int j;
	m=714025;
    ia=1366;
    ic=150889;
    rm=1/m;
   // iff=0;

	if((idum < 0) || (iff==0))
	//if (idum < 0)
	{
		iff=1;
		idum=fabs(fmod((ic-idum),m));

		for (j=0;j<97;j++)
		{
			idum=fmod(ia*idum+ic,m);
			ir[j]=idum;
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

	iy=ir[j];

	rtemp=iy*rm;
		 
	//printf("rtemp %f\n",rtemp);
	//printf("iy %f\n",iy);

	idum=fmod(ia*idum+ic,m);
	
	//printf("idum1 %f\n",idum);

	ran1[0]=rtemp;
	ran1[1]=idum;
	ran1[2]=iff;

	ir[j]=idum;
	 //printf("rtemp1 %f\n",rtemp);
	 
	return (ran1);
	 

}


long double *abscof(long double egamma)
{

	// munaip: Photoelectric absorption coefficient — represents gamma-ray photons completely absorbed by atoms.
// munaic: Compton interaction coefficient — represents gamma-ray photons scattered by electrons after losing energy.
	static long double abscof1[3]; long double ray=25.4,munaip,munaic,munai;

	if(egamma >= 0.55)
	{
		munaip=0.016-0.0283*egamma+0.0138*(pow(egamma,2));
	}

	else
	{
		munaip=0.115-0.426*egamma+0.42*(pow(egamma,2));
	}

	if(egamma >= 0.55)
	{
		munaic=0.0355-0.0222*(egamma)+0.00772*(pow(egamma,2))-0.00102*(pow(egamma,3));
	}

	else
	{
		munaic=0.051-0.0731*(egamma)+0.0507*(pow(egamma,2));
	}

	munai=munaip+munaic;
	abscof1[0]=munai;
	abscof1[1]=munaip;
	abscof1[2]=munaic;

	return (abscof1);
}


long double photopeak(long double egamma)
{

	// It estimates the probability that a gamma ray deposits its full energy inside the detector crystal, contributing to the photopeak.

	long double ray=25.4,chord,ebar,pphot,pcomp,pphot1,pcomp1,ppic,munaip,munaic,munai,eint,*abscof2,*abscof3;


	// egamma-> Energy of the incoming gamma ray
	// ray-> Detector radius used in the model
	// chord-> Characteristic distance travelled inside the detector
	// ebar-> Estimated representative energy of the Compton-scattered photon
	// pphot,pcomp-> Relative photoelectric and Compton contributions at energy ebar
	// pphot1,pcomp1-> Corresponding contributions at the original energy
	// ppic-> Final calculated photopeak probability
	// munai,munaip,munaic-> Total, photoelectric and Compton attenuation coefficients
	// eint-> Probability of at least one interaction over the specified chord
	// abscof2,abscof3-> Pointers used to receive values from abscof()


	chord=2*ray/3;
	ebar=egamma*(1.0-0.985*exp(-.425/egamma));

	abscof2=abscof(ebar);
	munai=abscof2[0];
	munaip=abscof2[1];
	munaic=abscof2[2];


	pphot=munaip/munai; // contribution due to photoelectric 

	pcomp=munaic/munai; // contri of compton

	eint=1-exp(-chord*munai);


	// calculating the final photopeak
	abscof3=abscof(egamma);
	munai=abscof3[0];
	munaip=abscof3[1];
	munaic=abscof3[2];

	pphot1=munaip/munai;
	pcomp1=munaic/munai;
	ppic=pphot1+pcomp1*pphot*eint/(1-pcomp*eint);
	// Calculate photopeak probability:
// 1. pphot1: Direct photoelectric absorption of the incoming gamma ray.
// 2. The second term accounts for Compton scattering followed by photoelectric absorption.
//    Repeated Compton interactions are also considered using a geometric series:
//    A + A*x + A*x^2 + ... = A/(1-x)
//    where x = pcomp*eint (probability contribution of each additional Compton interaction).
//    This is an approximate model using representative scattered energy and fixed chord length.

ppic = pphot1 + pcomp1*pphot*eint/(1-pcomp*eint);

	return (ppic);
}
















			





 

	


