#pragma once
#ifndef _STEM_H_
#define _STEM_H_
#include "organ.h"

class CStem: public COrgan
{
public:
	CStem();
	CStem(const TInitInfo ,int);
	~CStem();

	double get_length() {return length;}
	double get_diameter() {return diameter;}
	void set_length(double x) {length=x;}
	void set_diameter(double x) {diameter=x;}
	void set_deadmass(double x) {deadmass=x;}
	double get_deadmass(){return deadmass;}
	double get_deadCHOpool();
	double get_deadNmass(){return deadNmass;}
	double get_deadNpool(){double temp = deadNpool; deadNpool = 0; return temp;}

	void update(double area, int icur, CDevelopment*, const TInitInfo, const TWeather&, double Tlagleaf);
	void senescence(double area, CDevelopment* ,const TInitInfo info, double Tlagleaf);

private:
	int ID = 0;
	double length = 0.;
	double diameter = 0.;
	double deadmass = 0.; //total stem mass at time of senescence
	double deadCHOpool = 0.; // amount of stem mass to be put back into C_pool
	double deadNmass = 0.;// total N trapped in senesced stem and lost by plant
	double deadNpool = 0.;// temp N recovered from stem
	
};	
#endif