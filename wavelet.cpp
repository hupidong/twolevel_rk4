#include <math.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include "complex.h"
const double pi=acos(-1.0);
void main()
{
	using namespace std;
	int i=0,j=0,k=0;
	char s1[200]="E:\\My Dbank\\My Buffer Files\\C++\\TwoLevelHHG\\ODE45_DensityMatrix_Gauss_withoutchirp_xi=4\\";
	char s2[200]="0";
	char *q;
    q=s2;
	strcpy(q,s1);
	double omega_0, omega_1,rabbi_1,mu, T, tao, tstart, dt;
	long N;
	fstream para;
	strcat(q,"relative_parameters.dat");
	para.open(q,ios::in|ios::binary);
	if(!para)
	{
		cout<<"Can not open file "<<q<<endl;
		system("pause");
		exit(0);
	}
	para.read((char *)&omega_0,sizeof(double));
	para.read((char *)&omega_1,sizeof(double));
	para.read((char *)&rabbi_1,sizeof(double));
	para.read((char *)&mu,sizeof(double));
	para.read((char*)&N,sizeof(long));
	para.read((char*)&dt,sizeof(double));
	para.read((char*)&T,sizeof(double));
	para.read((char *)&tstart,sizeof(double));
	para.close();
	strcpy(q,s1);
	N -= 1;
	tao=6.0;
	double *t;
	t=new double[N];
	for(i=0;i<N;i++)
	{
		t[i]=tstart+i*dt;
	}

	double wti,wtf,wdt,wwi,wwf,wdw;//小波分析的起始时间，结束时间，步长，最小频率，最大频率，频率间隔
	wti=-3*T;
	wtf=-wti;
	wdt=0.02*T;
	wwi=0.0*omega_1;
	wwf=80*omega_1;
	wdw=0.2*omega_1;
	int mw,nt;
	nt=int((wtf-wti)/wdt)+1;
	mw=int((wwf-wwi)/wdw)+1;
    cout<<"nt="<<nt<<endl;cout<<"mw="<<mw<<endl;
    double *tcenter,*omg;
	tcenter=new double[nt];
	omg=new double[mw];
	fstream waveletplot;
	double unit=0;
	strcat(q,"wtplot.dat");

	waveletplot.open(q,ios::out|ios::binary);
	if(!waveletplot)
	{
		cout<<"Can not open file "<<q<<endl;
		system("pause");
		exit(0);
	}
	waveletplot.write((char *)&nt,sizeof(int));
	waveletplot.write((char *)&mw,sizeof(int));
	for(i=0;i<nt;i++)
	{
		tcenter[i]=wti+i*wdt;
		unit=tcenter[i]/T;
	    waveletplot.write((char *)&unit,sizeof(double));
	}
	for(i=0;i<mw;i++)
	{
		omg[i]=wwi+i*wdw;
		unit=omg[i]/omega_1;
	    waveletplot.write((char *)&unit,sizeof(double));
	}
	waveletplot.close();
	strcpy(q,s1);
	complex **wavelet;
    wavelet=(complex **)malloc(nt*sizeof(complex *));
	for(i=0;i<nt;i++)
	{
		wavelet[i]=(complex *)malloc(mw*sizeof(complex));
	}
	double *da,*datime;
	da=new double[N];
	datime=new double[N];

	fstream dipole;
	strcat(q,"dipole.dat");

	dipole.open(q,ios::in|ios::binary);
	if(!dipole)
	{
		cout<<"Can not open file "<<q<<endl;
		system("pause");
		exit(0);
	}
	for(i=0;i<N;i++)
	{
//		cout<<i<<endl;
		datime[i]=(tstart+i*dt)/T;
		dipole.read((char *)&da[i],sizeof(double));
	}
	dipole.close();
	strcpy(q,s1);
/*	Engine *ep;
	mxArray *tt=NULL,*result=NULL,*ts;
	if(!(ep=engOpen("\0")))
	{
		cout<<"Can not start matlab engine"<<endl;
	}
	tt=mxCreateDoubleMatrix(1,N,mxREAL);
	ts=mxCreateDoubleMatrix(1,N,mxREAL);
	memcpy((void *) mxGetPr(tt),(void*)datime,N*sizeof(double));
	engPutVariable(ep,"tt",tt);
	memcpy((void *) mxGetPr(ts),(void*)da,N*sizeof(double));
	engPutVariable(ep,"ts",ts);
	engEvalString(ep,"plot(tt,ts);");
*/
	ofstream out;
	cout<<"ok"<<endl;
	strcat(q,"\\dipoleawt.txt");
	out.open(q);
	if(!out)
	{
		cout<<"Can not open file "<<q<<endl;
		system("pause");
		exit(0);
	}
	double omgt,z;
    for(i=0;i<nt;i++)
	{
		for(j=0;j<mw;j++)
		{
			wavelet[i][j]=0.0;
			for(k=0;k<N;k++)
			{
				omgt=omg[j]*(t[k]-tcenter[i]);
				wavelet[i][j]=wavelet[i][j]+da[k]*sqrt(omg[j]/tao)*complex(cos(omgt),sin(omgt))*exp(-omgt*omgt/2.0/tao/tao)*dt;
			}
			z=mod(wavelet[i][j]);
			out<<z<<"  ";
		}
		out<<endl;
	}//*/
	dipole.close();////*
	out.close();
	free(wavelet);
	delete[] t,da,omg,tcenter;//*/
}