//
//  main.m
//  useful
//
//  Created by Phil Ahrenkiel on 3/7/20.
//  Copyright © 2020 Phil Ahrenkiel. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#include "arr.hpp"

darr1 getVect(double mX,double mY,uiarr2 &coefA);
darr1 getVect(double mX,double mY,uiarr2 &coefA)
{
	std::size_t nCoef=coefA.size(0);
	darr1 d(nCoef);
	for(std::size_t k=0;k<nCoef;k++)
	{
		uint nX=coefA(k,0),nY=coefA(k,1);
		double fX=1,fY=1;
		if(nX>0)fX=xpowy(mX,nX);
		if(nY>0)fY=xpowy(mY,nY);
		d(k)=fX*fY;
	}
	return d;
}


darr1 getVect(double mX,uiarr1 &coefA);
darr1 getVect(double mX,uiarr1 &coefA)
{
	std::size_t nCoef=coefA.size();
	darr1 d(nCoef);
	for(std::size_t k=0;k<nCoef;k++)
	{
		uint nX=coefA(k);
		double fX=1;
		if(nX>0)fX=xpowy(mX,nX);
		d(k)=fX;
	}
	return d;
}

void test1();
void test1()
{
	long nRadPoints=6,nAngPoints=4;
	long nPow=3;
	double r0=1,rRange=0.4;
	long nDataPoints=nRadPoints*nAngPoints;

	//double cPhi=cos(Pi*phi/180),sPhi=sin(Pi*phi/180);
	long i,j,k=0,l;
	darr2 phi(2,nDataPoints);
	for(i=0;i<nRadPoints;++i)
	{
		double mRad=r0*(1.+rRange*(i-(nRadPoints-1.)/2.));
		for(j=0;j<nAngPoints;++j)
		{
			double mAng=2*Pi*j/nAngPoints;
			double mX=mRad*cos(mAng);
			double mY=mRad*sin(mAng);
			phi(0,k)=mX;
			phi(1,k)=mY;
			k++;
		}
	}
	std::cout<<"phi:\n"<<phi<<"\n";
	long nCoef=sqr(nPow+1);
	uiarr2 coef(nCoef,2);
	k=0;
	unsigned int nX,nY;
	for(nY=0;nY<=nPow;nY++)
		for(nX=0;nX<=nPow;nX++)
		{
			coef(k,0)=nX;
			coef(k,1)=nY;
			k++;
		}
	std::cout<<"coef:\n"<<coef<<"\n";

	double m=0,m2=0;
	for(l=0;l<nDataPoints;l++)
	{
		m+=phi(0,l)+phi(1,l);
		m2+=sqr(phi(0,l))+sqr(phi(1,l));
	}
	m/=2*nDataPoints;m2/=2*nDataPoints;
	double xyRange=sqrt(m2-m*m);

	//
	darr2 d(nCoef,nDataPoints);
	for(l=0;l<nDataPoints;l++)
	{
		darr1 dVect=getVect(phi(0,l)/xyRange,phi(1,l)/xyRange,coef);
		for(k=0;k<nCoef;k++)
			d(k,l)=dVect(k);
	}
	std::cout<<"d:\n"<<d<<"\n";

	darr2 u,s,v;
	bool res=d.SVD(u,s,v);
	darr2 dp=u*s*v.transpose();
	std::cout<<"d':\n"<<dp<<"\n";


	darr2 invd=d.pseudoinverse(1e8);
	//std::cout<<"id:\n"<<invd<<"\n";
	darr2 p=d*invd;
	std::cout<<"unit:\n"<<p<<"\n";
	
	//find matrix
	double diffShift0=100,diffShiftRange=0.4,diffShiftAng=0;//Pi/2;
	darr2 diffShift(2,nDataPoints);
	k=0;
	for(i=0;i<nRadPoints;++i)
	{
		double mRad=diffShift0*(1.+diffShiftRange*(i-(nRadPoints-1.)/2.));
		for(j=0;j<nAngPoints;++j)
		{
			double mAng=2*Pi*j/nAngPoints+diffShiftAng;
			double mX=mRad*cos(mAng);
			double mY=mRad*sin(mAng);
			diffShift(0,k)=mX;
			diffShift(1,k)=mY;
			k++;
		}
	}
	darr2 mComp=diffShift*invd;
	std::cout<<"matrix:"<<mComp<<"\n";

	darr2 diffShiftP=mComp*d;
	std::cout<<"diffShift:"<<diffShift<<"\n";
	std::cout<<"diffShift':"<<diffShiftP<<"\n";
}

void test2();
void test2()
{
	std::size_t m=4;
	std::size_t n=32;
	
	darr2 d(m,n);
	std::size_t i,j;
	for(i=0;i<m;++i)
		for(j=0;j<n;++j)
			d(i,j)=2*rndom()-1;

	darr2 u,s,v;
	bool res=d.SVD(u,s,v);

	//darr2 di=d.pseudoinverse(1e8);
	//darr2 p=d*di;
	std::cout<<s<<"\n";
}


void test3();
void test3()
{
	long nPoints=9;
	long nPow=5;
	double r0=1,rRange=0.4;

	//double cPhi=cos(Pi*phi/180),sPhi=sin(Pi*phi/180);
	long i,j,k=0,l;
	darr2 phi(1,nPoints);
	for(i=0;i<nPoints;++i)
	{
		double mRad=r0*(1.+rRange*(i-(nPoints-1.)/2.));
		phi(0,k)=mRad;
		k++;
	}
	std::cout<<"phi:\n"<<phi<<"\n";
	long nCoef=nPow+1;
	uiarr1 coef(nCoef);
	k=0;
	unsigned int nX;
	for(nX=0;nX<=nPow;nX++)
	{
		coef(k)=nX;
		k++;
	}
	std::cout<<"coef:\n"<<coef<<"\n";

	double m=0,m2=0;
	double magPhi=0;
	for(l=0;l<nPoints;l++)
	{
		if(fabs(phi(0,l))>magPhi)magPhi=fabs(phi(0,l));
		m+=phi(0,l);
		m2+=sqr(phi(0,l));
	}
	m/=nPoints;m2/=nPoints;
	double xRange=magPhi;//sqrt(m2-m*m);

	//
	darr2 d(nCoef,nPoints);
	for(l=0;l<nPoints;l++)
	{
		darr1 dVect=getVect(phi(0,l)/xRange,coef);
		for(k=0;k<nCoef;k++)
			d(k,l)=dVect(k);
	}
	std::cout<<"d:\n"<<d<<"\n";

	darr2 u,s,v;
	//bool res=d.SVD(u,s,v);
	//std::cout<<"S:\n"<<s<<"\n";

	//darr2 dp=u*s*v.transpose();
	//std::cout<<"d':\n"<<dp<<"\n";


	darr2 invd=d.pseudoinverse(1e8);
	//std::cout<<"id:\n"<<invd<<"\n";
	darr2 p=d*invd;
	std::cout<<"unit:\n"<<p<<"\n";
	
	//find matrix
	double diffShift0=100,diffShiftRange=0.4,diffShiftAng=0;//Pi/2;
	darr2 diffShift(2,nPoints);
	k=0;
	for(i=0;i<nPoints;++i)
	{
		double mRad=diffShift0*(1.+diffShiftRange*(i-(nPoints-1.)/2.));
		diffShift(0,k)=mRad;
		diffShift(1,k)=-mRad;

		k++;
	}
	darr2 mComp=diffShift*invd;
	std::cout<<"matrix:"<<mComp<<"\n";

	darr2 diffShiftP=mComp*d;
	std::cout<<"diffShift:"<<diffShift<<"\n";
	std::cout<<"diffShift':"<<diffShiftP<<"\n";
}

void test4();
void test4()
{

	darr2 d(7,7);
	d=d.icol()+d.irow();
	//const darr2 A(3,3);
	//A=2;
	std::cout<<"d:\n"<<d<<"\n";
}

int main(int argc, const char * argv[]) {
	@autoreleasepool {
	    // Setup code that might create autoreleased objects goes here.
	}
	test3();
	return NSApplicationMain(argc, argv);
}
