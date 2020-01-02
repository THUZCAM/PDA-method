//PDA method for smFRET
//by changing the gaussian distribution through (aveR,sigmaR), the distance distribution between smFRET donor and acceptor is converted into the smFRET efficiency and can be compared with experimental measurements.
 
#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#define NUM 125
#define TOTAL 40
double avr=0.3,avg=0.3,beta=1.17,alpha=0.11,R0=51;//emperical values for smFRET setup
double pi=3.14;

//calculate the value of F!
double fact(int F)
{
	double product=1;
	for(int i=1;i<=F;i++) product=product*i;
	return product;
}

//calculate a gaussian distribution
double Pepsilong(double epsilong, double aveR, double sigmaR)
{
	double tempP=R0/6/sqrt(2*pi)/sigmaR*pow(beta,1.0/6)/pow(1-epsilong,2)*pow(1.0/(1-epsilong)-(1+alpha),-7.0/6);
	tempP=tempP*exp(-pow(R0*pow(beta,1.0/6)*pow(1.0/(1-epsilong)-(1+alpha),-1.0/6)-aveR,2)/2/pow(sigmaR,2));
	return tempP;
}

//calculate a bionomial distribution
double Pcondition(double epsilong, int F, int FRT)
{
	if(FRT<0) return 0;
	double tempP=fact(F)/fact(FRT)/fact(F-FRT)*pow(epsilong,FRT)*pow(1-epsilong,F-FRT);
	return tempP;
}

//main code
void main()
{
	int i,j,k,l,N,num;
	FILE *fp1;
	double Pr,Pg,Pe,Psignal[TOTAL+2];
	int PN[NUM];
	double sum1,sum2,eff,epsilong;
	int total=TOTAL;
	double aveR=68,sigmaR=3;//the average and variance of the gaussian distribution assumed for distance between donor and acceptor

	i=1;
	fp1=fopen("intensity.txt","r");//input background photon distribution
	while (!feof(fp1))
	{
		fscanf(fp1,"%d",&num);
		PN[i++]=num;
	}
	fclose(fp1);

	sum2=0;epsilong=0.99;//the maximal FRET efficiency
	for(l=0;l<90;l++)//calculate apparent smFRET efficiency
	{
		Pe=Pepsilong(epsilong,aveR,sigmaR);
		epsilong=epsilong-0.01;eff=0;
		for(k=0;k<=total;k++)
		{
			sum1=0;
			for(N=30;N<=100;N++)//total number of photons
			{
				for(i=0;i<=5;i++)
				{
					Pr=pow(avr,i)*exp(-avr)/fact(i);//calculate a Poisson distriubtion
					for(j=0;j<=5;j++)
					{
						Pg=pow(avg,j)*exp(-avg)/fact(j);
						sum1=sum1+PN[N+1]*Pr*Pg*Pcondition(epsilong,N-i-j,int(eff*N)-i)*Pe;
					}
				}
			}
			eff=eff+1.0/total;
			if(l==0) Psignal[k+1]=sum1;
			else Psignal[k+1]=Psignal[k+1]+sum1;
			sum2=sum2+sum1;
		}
	}
 
	fp1=fopen("PDA.txt","w");//output the PDA results
	for(k=0;k<=total;k++)
	{
		fprintf(fp1,"%e %e\n", double(k)/total, Psignal[k+1]/sum2);
	}
	fclose(fp1);
}
