#include "StdAfx.h"
#include "thermaltime.h"
#include <math.h>

CThermalTime::CThermalTime(void)
{
	Tcur = 25.0; //set to air or soil temperature?
	Tbase = 0.0;
	Topt = 27.98; //as per potato leaf appearance rate
	Tmax = 40.46; // as per potato leaf appearance rate
	actualAge = sum = 0.0;
	dTmpr = 0.0;
	timeStep = 60.0;
}

double CThermalTime::Agefraction(const TInitInfo info, double age, double T, double duration, double Tlagleaf)
{
	/**********************/
	//returns fraction of organ(leaf) still capable of expanding
	// Based on modified Ng and Loomis leaf expansion model Fleisher and Timlin (2007)
	/*********************/
	double R = 0.,x = 0.;
	//if (age == 0) age = Page(info, T, Tlagleaf)*info.timeStep/(24*60);
	//if (13 < Tlagleaf < 17) age = 0.448;
	R = 0.5586*exp(-1.3563*age)+0.4431*exp(-0.4813*age);
	if (age >= duration) R = 0.;
	if (R > 1) R = 1.;
	if (R < 0) R = 0.;
	if (age <= 0.5) R = 1.;
	x=1.;
	return R;
	//if (age<7)
	//	dblAgefraction = 1;
	//else if (age<12)
	//	dblAgefraction = 1 - (0.9/5) * (age-7);
	//else
	//	dblAgefraction = 0.1 - (0.1/3) * (age-12);
	//if (dblAgefraction < 0) dblAgefraction = 0;
}


double CThermalTime::beta_fn(double t, double R_max, double t_m, double t_e)
{
//Generalized Temperature Response Model
	double deltaT1=0., deltaT2 = 0., alpha = 0.;
	const double t_b = 0.;
	deltaT1 = t_m - t_b;
	deltaT2 = t_e - t_m;
	alpha = (deltaT1/deltaT2);
	if (deltaT1 < 0 || deltaT2 < 0 ) return -1;
	else return R_max*pow((t-t_b)/deltaT1,alpha)*(t_e-t)/deltaT2;
}
double CThermalTime::Cage(const TInitInfo info)
{
	//Compute chronological aging
	double dt = info.timeStep/(24.*60.);
	return dt;
}
double CThermalTime::Page(const TInitInfo info, double T, double Tlagleaf) 
//Compute physiological age based on yesterday's air temperature and 
// divides into day fraction
{
	double R = 0.;
	//Tlagleaf = 15.3;
	//4 parameter log-normal curve:
	//R = (0.4287 + 0.5915*exp(-0.5*pow(log(T/23.97)/0.201,2)))*info.timeStep/(24*60);
	//R = (0.0292*Tlagleaf + 0.0309) * info.timeStep/(24*60);
	R = (0.0292*Tlagleaf+0.0309)*info.timeStep/(24*60); 
	//if ((Tlagleaf >= 14) && (Tlagleaf <= 16)) R = R - 0.034*info.timeStep/(24*60);//Temporary tweek - at 17/12 treatment, underpredicting expansion rate

	if (R < 0) R = 0.;
	if (R > 1) R = 1.;
	if (T<4) R = 0.; // assumed no growth at very low T
	//if (T==12.7)
	//	R = 0.4512;
	//else if (T==16)
	//	R = 0.4656;
	//else if (T==18.3)
	//	R = 0.7289;
	//else if (T==20.)
	//	R = 0.8256; //linearly interpolated
	//else if (T==21.3) 
	//	R = 0.8884;
	//else if (T==26.3)
	//	R = 0.9845;
	//else 
	//	R = 0.6192;
	return R;
}
void CThermalTime::initialize(double step)
{
	/*************/
	/* MAIZEIM Code
	/* Not sure method and/or variables are used in SPUDSIM
	/***********/
	Tcur = 0.0;
	timeStep = step;
	actualAge = sum = 0.0;
	add(Tcur);
}

void CThermalTime::add(double x)
{
	Tcur = x;
	double dD = timeStep/MINUTESPERDAY;
	if (Tcur <= Tbase||Tcur >= Tmax)
	{
		dTmpr = 0.0;
	}
	else
	{
		dTmpr = (Tcur-Tbase);
	}

	sum += (dTmpr*dD);
	actualAge += dD;
}

double CThermalTime::Tfactor(const TInitInfo info, double T, double Tlagleaf)
/*************************/
//estimates temperature influence on fraction of leaf capable of expanding
/* Based on modified Ng and Loomis, Fleisher and Timlin (2007)
/*************************/
{
	
	double R;
	//Tlagleaf = 15.3;
	//if (info.timeStep == 1440) //1 day increment use the following
	{	//4 parameter log-normal curve:
		R = (0.3414 + 0.5927*exp(-0.5*pow(log(Tlagleaf/27.65)/0.3075,2.)));
		/* use interpolated values for 2009 paper except 14/10, 23/18 - too suppressed - need way to make expansion less sensitive!*/
		if (Tlagleaf <= 12.9) 
			R = (Tlagleaf-4)*0.435 / (12.9-4);
		else if (Tlagleaf <=15.7)
			R = (Tlagleaf-12.9)*(0.53-0.435)/(15.7-12.9)+0.435;
			//R = 0.8;
		else if (Tlagleaf <=19.4)
			R = (Tlagleaf-15.7)*(0.64-0.53)/(19.4-15.7)+0.53;
			//R = 0.8;
		else if (Tlagleaf <= 22.1)
			R = (Tlagleaf-19.4)*(0.72-0.64)/(22.1-19.4)+0.64;
		else if (Tlagleaf <= 26.8) 
			R = (Tlagleaf-22.1)*(0.85-0.72)/(26.8-22.1)+0.72;
		else if (Tlagleaf <= 31.9) 
			R = (Tlagleaf-26.8)*(0.91-0.85)/(31.9-26.8)+0.85;
		else
			R = (Tlagleaf-31.9)*(-0.91)/(40-31.9)+0.91;
	}
	/*if (info.timeStep == 5) //5 min time increment//Apparent nonlinear T dependency on time-step...
	{
		if (T <= 12.9) R = (T-4)*0.435 / (12.9-4);
		if (12.9 < T <= 15.7) R = (T-12.9)*(0.53-0.435)/(15.7-12.9)+0.435;
		if (15.7 < T <= 19.4) R = (T-15.7)*(0.64-0.53)/(19.4-15.7)+0.53;
		if (19.4 < T <= 22.1) R = (T-19.4)*(0.72-0.64)/(22.1-19.4)+0.64;
		if (22.1 < T <= 26.8) R = (T-22.1)*(0.85-0.72)/(26.8-22.1)+0.72;
		if (26.8 < T <= 31.9) R = (T-26.8)*(0.91-0.85)/(31.9-26.8)+0.85;
		if (T > 31.9) R = (T-31.9)*(-0.91)/(40-31.9)+0.91;
	}
	*/
	//R = (0.3414 + 0.5927*exp(-0.5*pow(log(Tlagleaf/27.65)/0.3075,2)));

	if (R<0) R = 0.;
	if (R>1) R = 1.;
	return R;
}



void CThermalTime::update(double Tmpr, double step)
{
	timeStep = step;
	add(Tmpr);
}



CThermalTime::~CThermalTime(void)
{
}
