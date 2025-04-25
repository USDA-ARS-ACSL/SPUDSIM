#include "stdafx.h"
#include "math.h"
#include "organ.h"
#define CO2_MW 44.0098
#define C_MW 12.011
#define CH2O_MW 30.03


COrgan::COrgan()
{
	Temperature=25.0; //should be set to equal soil temperature
	CH2O=Drymass=drymass_at_emergence=0.;
	dblCage = dblPage = 0.;
	dblAgefraction = dblTfactor = 0.;
	Potgrowth = Potdrymass = 0.;
	GrowthDuration=15.; 
	Longevity=50.;
	GDD = NULL;
	GDD = new CThermalTime();
	maintCoeff = 0.0206; //gCH2O g-1DM day-1 at 25C as in Ng and Loomis
	initialized = false;
}
COrgan::COrgan(const TInitInfo& info)
{
	initInfo = info;
	Temperature=25.0; //should be set to equal soil temperature
	CH2O=Drymass=drymass_at_emergence=0.;
	dblCage = dblPage = 0.;
	dblAgefraction = dblTfactor = 0.;
	Potgrowth = Potdrymass = 0.;
	GrowthDuration=15.; 
	Longevity=50.;
	GDD = NULL;
	GDD = new CThermalTime();
	maintCoeff = 0.0206; //gCH2O g-1DM day-1 at 25C as in Ng and Loomis
	initialized = false;
}

COrgan::~COrgan()
{
	if (GDD != NULL) delete GDD;
}

void COrgan::initialize()
{
	if (GDD == NULL) GDD = new CThermalTime();
	//GDD->initialize(initInfo.timeStep);
	GDD->initialize(60.);
	dblCage = dblPage = 0.;

}
//void COrgan::update() //is this routine used?
//{
//	GDD->add(temperature);
//	age = GDD->get_actualAge();
//	physAge=GDD->get_sum();
//	mass = CH2O*(C_MW/CH2O_MW)/0.4; // C content = 40%, hard coded for now
//}

double COrgan::get_maintResp()
{
	return maintCoeff*(Drymass-drymass_at_emergence);
}


void COrgan::update(int icur, const TInitInfo info, const TWeather& weather, double Tlagleaf) //update physiological age of nodal unit
{
	/***********************************/
	/* Method called during leaf.update*/
	/* Calls up member functions from Thermal time
	/* to estimate chronological and physiological aging
	/*  and related parameters that influence leaf expansion rate*/
	/***********************************/
	double x;
	dblCage += GDD->Cage(info);
	x=this->dblCage;
	//if (x==0) dblPage += GDD->Page(info, weather.airT, Tlagleaf); //only update Page at each 24 h of leaf life
	//if ((x - floor(x)) <= info.timeStep/(24*60)) dblPage += GDD->Page(info, weather.airT, Tlagleaf);
	dblPage +=GDD->Page(info, weather.airT, Tlagleaf);
	dblAgefraction = GDD->Agefraction(info, dblPage,weather.airT, GrowthDuration, Tlagleaf);
	dblTfactor = GDD->Tfactor(info, weather.airT, Tlagleaf);
}

void COrgan::import_CH2O(double dCH2O)
{
	CH2O += dCH2O;
	Drymass = CH2O*(C_MW/CH2O_MW)/0.4; // C content = 40%, hard coded for now ?
}

void COrgan::import_N(double dN)
{
//	N += dN;
}

void COrgan::respire()
// this needs to be worked on
// currently not used at all
{
	double Rm = 0.02; //maintenance respiration
	double Ka = 0.1; //growth respiration
	CH2O -= Ka*CH2O + Rm*CH2O;
}




