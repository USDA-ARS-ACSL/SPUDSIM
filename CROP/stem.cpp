#include "stdafx.h"
#include "stem.h"

CStem::CStem()
:COrgan(), length(0.0), diameter(0.0) {}

CStem::CStem(const TInitInfo info, int n): COrgan()
{
	ID = n;
	diameter = 0.0;
	deadmass = 0.;
	deadCHOpool = 0.;
	//set_drymass(0.2/4*(info.timeStep/(24*60)));//from Spar04 measurements (4 nodes have stem mass of 0.2 g)
	//set_CH2O(0.2/4*(info.timeStep/(24*60)));
	if (info.bigleaf == 0)
	{
		if (n <= 3) //initial stemmass higher at emergence than later in season
		{
			set_drymass(0.2/4.);
			set_drymass_at_emergence(0.2/4.);
			set_CH2O(0.2/4.);
		}
		else
		{	set_drymass(0.0025);
			set_CH2O(0.0025);
		}
		length = 0.1; // made up!?
	}
	else{
		set_drymass(0.2);
		set_drymass_at_emergence(0.2);
		set_CH2O(0.2);
		length = 0.4;
	}
	set_optimumNitrogenRatio(0.05);
	set_minimumNitrogenRatio(0.04);
	set_actualNitrogenRatio(0.04);
	set_optimumNitrogenRatio(0.045);//from SUBSTOR
	set_actualNitrogenRatio(0.045);//from SUBSTOR
	set_minimumNitrogenRatio(0.025);//from SUBSTOR
	set_currentNitrogenAmount(this->get_drymass()*this->get_actualNitrogenRatio()); //SIMPOTATO
}
CStem::~CStem() {}

double CStem::get_deadCHOpool()
{
	double x = 0.;
	x = this->deadCHOpool;
	this->deadCHOpool = 0.;
	return x; 
}

void CStem::update(double area, int icur, CDevelopment * dv, const TInitInfo info, const TWeather &weather, double Tlagleaf)
{
	/*************************/
	/* Method only called when simulating BigLeaf (SIMPOTATO) Approach
	/* Function is to simulate stem senescence only */
	/**************************/
	senescence(area, dv,info, Tlagleaf);
}

void CStem::senescence(double area, CDevelopment * dv, const TInitInfo info, double T)
{
	/*******************************/
	/* Method only called for bigleaf model only
	/* Consult SIMPOTATO legacy code for specific information*/
	/*******************************/
	double slan = 0.;
	double stloss = 0.;
	if (dv->get_istage() < 3) stloss = 0.;
	else if (dv->get_istage() == 3)
	{
		/*
		slan = area/(500*3/dv->get_xstage())/(24*60/info.timeStep);
		if (dv->get_xstage()>4) slan = area/(500/dv->get_xstage())/(24*60/info.timeStep);
		stloss = slan*2*info.G4;
		if (stloss > this->get_drymass()) stloss = this->get_drymass();
		*/ //doesn't work well - too harsh on senescene, need to reevaluate!
	}
	deadmass += stloss;
	deadCHOpool = stloss;
	deadNmass = 0.25*stloss*this->get_actualNitrogenRatio();
	deadNpool = 0.75*stloss*this->get_actualNitrogenRatio();
	this->set_drymass(this->get_drymass()-deadCHOpool);
	this->set_CH2O(this->get_CH2O()-deadCHOpool);
	this->set_currentNitrogenAmount(this->get_currentNitrogenAmount()-deadNmass-deadNpool); //pool will be added back in during plant.cpp update-node method

}