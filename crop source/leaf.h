#pragma once
#ifndef _LEAF_H_
#define _LEAF_H_
#include "organ.h"
#include "weather.h"
#include "development.h"

class CLeaf: public COrgan
//leaf blade, call it "leaf" for convenience
{
public:
	CLeaf();
	CLeaf(const TInitInfo, int, int); // take current leaf rank and total leaf number to calculate potentialArea
	~CLeaf();
	bool isYoung() {return young;}
	bool isOld() {return old;}
	bool isInitiated() {return initiated;}
	bool isGrowing() {return growing;}
	bool isProlific() {return prolific;}
	bool isAging() {return aging;}
	bool isTerminated() {return terminated;}
	bool isDropped() {return dropped;}
	double get_area(){return area;}
	double get_deadmass(){return deadmass;}
	double get_deadNmass(){return deadNmass;}
	double get_deadCHOpool();
	double get_deadNpool(){double temp = deadNpool; deadNpool = 0; return temp;}
	double get_greenArea() {return greenArea;}
//	double get_senescentArea() {return senescentArea;}
	double get_psi_leafexpansion_stress() { return psi_leafexpansion_stress; }
	double get_potentialArea(){return potentialArea;}
	double get_potentialDWinc(){return potentialDWinc;}
	double get_nostresspotentialDWinc() {return nostresspotentialDWinc;}
	double get_potGrowth24(){return potGrowth24;}
	double get_area0(){return area0;}
	double get_maintResp(const TInitInfo info);
	double get_ageEffect(){return ageEffect;}
	double get_SLAleaf(){return SLAleaf;}
	double get_SLAmin() {return SLAmin;}
	int get_ID(){return ID;}
	//void set_environment(const TWeather & weather, double);//set pressure, PAR, CO2, T at leaf surface
	double get_length() {return length;}
	double get_SLA() {return SLA;}
	void set_deadmass(double x) {deadmass=x;}
	void set_area(double x) {area=x;}
	void set_area0(double x) {area0=x;}
	void set_length(double x) {length=x;}
	void set_SLA(const TInitInfo, double, double);
	void set_SLAleaf(double x) {SLAleaf=x;}

	double GTI(double);
	
	//Growth methods
	void expand(int icur, const TInitInfo, int potential, const TNitrogenStatus, double CH2O, CDevelopment*, double pdlwp, const TWeather&);
	void set_potentialArea(const TInitInfo);
	double get_LWPeffect(double predawn_psil, const TWeather&, const TInitInfo info); //from MAIZSIM that reduces potential leaf expansion based on pre-dawn water potential

	//Development methods
	void update(int icur, CDevelopment*, const TInitInfo, const TWeather&, const TNitrogenStatus, double Tlagleaf, int potential, double CH2O, double pdlwp, double swdf1, double C_pool_room);
	void senescence(int icur, CDevelopment* ,const TInitInfo info, const TWeather&, const TNitrogenStatus, double Tlagleaf, double swdf1, double C_pool_room);
	void senescence2(const TInitInfo info); //routine accounts for termination due to canopy overshading
	void calcLongevity();
	void set_ageEffect(const TInitInfo info); //set age effect on pgross parameters

private:
	TInitInfo initInfo;
	bool young, old, initiated, growing, prolific, aging, terminated, dropped;
	//int rank; //is this needed since nodal unit already has this value?
	int ID = 0; // used ID instead to keep consistent with nodal unit terminology
	double potGrowth24 = 0; //potential leaf area growth cfor current 24hours using yesterday's average Tair
	double potentialArea = 0; // potential leaf area with incremental increase
	double potentialAreainc = 0; // potential increase in leaf area
	double potentialDWinc = 0; //potential incremental increase in dry weight
	double potentialDW = 0; //potential leaf dry weight with incremental increase
	
	double nostresspotentialArea = 0; //as above, except these are potential increases with zero water stress
	double nostresspotentialAreainc = 0;
	double nostresspotentialDWinc = 0;
	double nostresspotentialDW = 0;


	double psi_leafexpansion_stress = 0; //effect of leaf water potential (0 to 1) on leaf expansion rate
	double senageArea = 0; // potential senescence due to leaf aging, cm2 (24 h total)
	double senstressArea = 0; //potential senescence due to leaf stress, cm2 (24 h total)
	double area = 0; // actual leaf area (even if leaf is dead, includes total leaf area at current time)
	double area0 = 0; //area at beginning of current day - this is used to estimate new leaf area production at each timestep
	double DWinc = 0; //actual incremental increase in dry weight
	double deadmass = 0; //total leaf mass at time of senescence
	double deadNmass = 0; //total leaf N 'trapped' in senesced leaf not recoverable
	double deadCHOpool = 0; // amount of leaf mass to be put back into C_pool
	double deadNpool = 0; // amount of senseced leaf N to be put back into N_available
	double greenArea = 0;  // 'green' leaf area
	double length = 0;
	double width = 0;
	double SLA = 0;
	double GDD = 0;
	double SLAnom = 0, SLAmin = 0, SLAmax = 0, SLAleaf = 0; //define acceptable variations in mass to weight ratios
	double Rmax_leaf = 0; //nonlinear beta function
	double T_base = 0, T_opt = 0, T_ceil = 0; //nonlinear beta function
	double dblMaxgrowth = 0; //parameter - max daily growth rate, cm2 cm2 leaf d-1
	double PAR = 0; //amount of par at leaf surface
	double leafT = 0; //leaf temperature
	double maintCoeff = 0; // maintenance respiration coefficient, g ch2o g-1 dm-1 d-1
	double maintResp = 0; //g CH2O d-1 lost to maintenance respiration
	double ageEffect = 0; //unitless (0 to 1) effect of leaf age on Pg parameters

	//Photosynthetic Parameters, values, coefficients, etc.
	double R_abs = 0, AirT = 0, CO2 = 0, RH = 0, Press = 0;

	//Water stress parameters from MAIZSIM

   
	double iter1 = 0, iter2 = 0,LeafHeight = 0;
};
#endif