#include "stdafx.h"
#include "tubers.h"
#include <math.h>
//#include <stdlib.h>
using namespace std;

CTubers::CTubers(void)
{
	etgt = dblSolids = tubmax = 0.;
	initiated = false;
	maintCoeff = 0.0206; //gCH2O g-1DM day-1 at 25C for mature potato leaves, Ng and Loomis; will decline with age as Burton, 1964
	this->set_CH2O(0.);
	this->set_drymass_at_emergence(0.);
	this->set_currentNitrogenAmount(0.);
	this->set_drymass(0.);
}

CTubers::~CTubers(void)
{
}

double CTubers::get_feedback(double Tubfraction)
{
	/*********************/
	/* Modified from Ng and Loomis*/
	/* Determins modifier for photosynthetic rate based on tuber sinkstrength*/
	/*  -sinkstrength related to fraction of CHO sent to tubers at preivous timestep*/
	/*********************/
	double x = 0.;
	if (Tubfraction < 0.1) x = 0.;
	else if (Tubfraction <= 0.25) x = (0.4-0)/(.25-.1)*(Tubfraction - 0.1)+0;
	else if (Tubfraction <= 0.55) x = (0.65-0.4)/(0.55-0.25)*(Tubfraction - 0.25)+0.4;
	else if (Tubfraction <= 1) x = (0.85-0.65)/(1-0.55)*(Tubfraction - 0.55)+0.65;
	else if (Tubfraction <= 1.5) x = (1.0-0.85)/(1.5-1)*(Tubfraction - 1)+0.85;
	else x = 1;
	return x;
}

double CTubers::get_maintResp(double age)
{
	double dif;
	dif = age;
	if (dif < 5) maintCoeff = 0.0206;
	else if ((dif >= 5) && (dif < 35)) maintCoeff = 0.0206-(0.0206-0.009245)/(5-34)*(5-dif);
	else if ((dif >= 35) && (dif < 66)) //TOO MANY ELSE IF?
	{
		maintCoeff = 0.009245-(0.009245-0.004622)/(35.-65.)*(35.-dif);
	}
	else if ((dif >= 66) && (dif <= 85))
	{
		maintCoeff = 0.004622-(0.004622-0.002311)/(66.-85.)*(66.-dif);
	}
	else maintCoeff = 0.001;

	return maintCoeff*get_drymass(); // override Ng and Loomis for now
}


double CTubers::Growth(CDevelopment* dv, double Tlag, const TInitInfo info)
{
	/********************************/
	/* Estimates tuber potential growth
	/* From Legacy SIMPOTATO Code from GROSUB routine*/
	/* Since individual tubers are not simulated in either bigleaf version, SIMPOTATO legacy code used for both bigleafs
	/********************************/

	double Ptubgro = 0.;
	double dt = info.timeStep / (24.*60.);
	if ((dv->get_istage()>=2)&&(dv->get_istage()<=4))// variables post tuber initiation
	{
		// max tuber growth rate from SIMGUI
		if ((Tlag > 5) && (Tlag < 14))
			etgt = -0.555 + 0.11 * Tlag;
		else if ((Tlag >= 14) && (Tlag <= 18))
			etgt = 1.;
		else if ((Tlag > 18) && (Tlag < 35))
			etgt = 2.059 - 0.0589 * Tlag;
		else if ((Tlag <= 5) || (Tlag >= 35))
			etgt = 0.;
		// calculate % dry matter in growing tubers based on dtt17 only should be based on water status, nitrogen status, and soil temp (from SIMGUI)
		if ((dv->get_DTT17() < 0.5) && (Tlag > 20)){
			dblSolids = dblSolids - (0.5 - dv->get_DTT17())*0.002;
			if (dblSolids < 0.1) dblSolids = 0.1;
		}
		else if (dv->get_DTT17() > 0.5) {
			dblSolids = dblSolids + 2. * (dv->get_DTT17() - 0.5) * 0.002;
			if (dblSolids > 0.26) dblSolids = 0.26;
		}
		// set maximum amount for daily tuber growth
		tubmax = min(info.G3 * etgt, info.G3 * pow(dv->get_partub(),2.));
		if (tubmax <0) tubmax = 0.;
		//tubmax = max(info.G3 * etgt, info.G3 * pow(dv->get_partub(),2)); // not getting enough post TI early tuber demand, particularly for high determinate crops 
		//tubmax = info.G3 * pow(dv->get_partub(),2);
	}	
	else { // pre tuber-initiation
			etgt = 0.;
			dblSolids = 0.12;
			tubmax = 0.;
	}
	Ptubgro = info.G1 * tubmax *dt; //adjust for timestep
	//Ptubgro = info.G1 * info.G3*dt;
	return Ptubgro;
}
