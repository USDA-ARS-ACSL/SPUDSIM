#pragma once
#ifndef _ORGAN_H_
#define _ORGAN_H_
#include "weather.h"
#include "thermaltime.h" 
#include "initInfo.h"
#include "nitrogenstatus.h"
#include "development.h"

struct TElement
{
	TElement() {CHO = 0; water = 0; nitrogen=0;}
	double CHO;
	double water;
	double nitrogen;
};

class COrgan
{
public:
	COrgan();
	COrgan(const TInitInfo&);
	virtual ~COrgan();
	virtual void import_CH2O(double);// import CHO from the reserve, virtual common metabolic reserve 
	virtual void import_N(double); // import N from the reserve
//	virtual double export_CHO(); //export CHO to the reserve, virtual common metabolic reserve 
//	virtual double export_N();
//	virtual void grow(double);
	virtual void respire();

	double get_Cage() {return dblCage;}
	double get_PAge() {return dblPage;}
	double get_drymass() {return Drymass;}
	double get_drymass_at_emergence() {return drymass_at_emergence;}
	double get_CH2O() {return CH2O;}
	bool get_Initialized() {return initialized;}
	double get_currentNitrogenAmount() {return currentNitrogenAmount;} // g N organ-1
	double get_actualNitrogenRatio() {return actualNitrogenRatio;} /*returns g g-1 organ-1*/
	double get_optimumNitrogenRatio() {return optimumNitrogenRatio;}
	double get_minimumNitrogenRatio() {return minimumNitrogenRatio;}
//	TElement * get_element() {return element;}
	double get_temperature() {return Temperature;}
	double get_longevity() {return Longevity;}
	double get_maintResp();
	double get_growthDuration() {return GrowthDuration;}
	double get_agefraction() {return dblAgefraction;}
	double get_Tfactor() {return dblTfactor;}
	TInitInfo get_initInfo() {return initInfo;}
	//CThermalTime * get_GDD() {return GDD;}
	//CThermalTime * get;

	virtual void set_age(double x) {dblCage = x;}
	virtual void set_physAge(double x) {dblPage=x;}
	virtual void set_drymass(double x) {Drymass=x;}
	virtual void set_drymass_at_emergence(double x) {drymass_at_emergence = x;}
	virtual void set_CH2O(double x) {CH2O=x;}
	virtual void set_maintCoeff(double x){maintCoeff = x;}
	virtual void set_currentNitrogenAmount(double x) {currentNitrogenAmount = x;}
	virtual void set_actualNitrogenRatio(double x) {actualNitrogenRatio = x;}
	virtual void set_minimumNitrogenRatio(double x) {minimumNitrogenRatio = x;}
	virtual void set_optimumNitrogenRatio(double x) {optimumNitrogenRatio = x;}
	virtual void set_initialized() {initialized = true;}
//	virtual void set_element(TElement * x) {element=x;}
	virtual void set_temperature(double x) {Temperature=x;}
	void set_longevity(double x) {Longevity=x;}
	virtual void set_growthDuration(double x) {GrowthDuration=x;}
	virtual void initialize();
	virtual void update(int icur, const TInitInfo info, const TWeather &weather, double Tdaylag);

private:
	COrgan(const COrgan&);
	TInitInfo initInfo;
	TNitrogenStatus NitrogenStatus;
	TElement * element;
	CThermalTime * GDD;
	bool initialized = 0; //set to true when emerged by 2DSOIL to allow pre-emergence simulation of soil transformations, etc

	double dblCage = 0; // chronological age of an organ, days
	double dblPage = 0; // physiological age accouting for temperature effect (in reference to endGrowth and lifeSpan, days)
	double dblAgefraction = 0; //fraction of leaf capable of growth
	double dblTfactor = 0; //T influence
	double maintCoeff = 0; //maintenance respiration coefficient for stem, root, and tuber only g CH2O g-1dm d-1
	double maintResp = 0; //mainteance respiration, g CH2o d-1
	
	double Drymass = 0; // biomass, g (dry weight)
	double drymass_at_emergence = 0; // initial biomass at emergence (g CH2O)
	double CH2O = 0; //glucose, MW = 180.18 / 6 = 30.03 g
	double Temperature = 0; // organ temperature, C
	double Longevity = 0; // life expectancy of an organ in days at optimal temperature (fastest growing temp), days
	double GrowthDuration = 0; // physiological days to reach the end of growth (both cell division and expansion) at optimal temperature, days
	double Potgrowth = 0; // Potential incremental growth computed for current time increment (g dry weight per time)
	double Potdrymass = 0; // Potential total dry mass at current time increment allowing for potential increase (g dry weight)
	
	//Nitrogen information
	double currentNitrogenAmount = 0; //total amount of N in g plant-1 for the organ
	double actualNitrogenRatio = 0; // actual current N amount, g g-1
	double minimumNitrogenRatio = 0; // minimum N requried to support new growth, g N g CHO-1
	double optimumNitrogenRatio = 0; // critical N required to support new growth at non-stressed level, g N g CHO-1
	// Add stress factors here

};
#endif
