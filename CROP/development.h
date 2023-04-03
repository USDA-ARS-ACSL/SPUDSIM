#pragma once
#ifndef _DEVELOPMENT_H_
#define _DEVELOPMENT_H_
#include "weather.h"
#include "initinfo.h"
#include "nitrogenstatus.h"

enum EPhase
  {
      //Seed, Juvenile, Inductive, preSilking, Reproductive, Maturity
	  Vegetative, Tuberlinear, Tuberdominant, Maturity, Preplanting, Germination, Emergence
  };

struct TEvent
{
public:
	TEvent() {day = 0;  done=false;}
	double day;
	bool done;
};

class CDevelopment
{
public:
    CDevelopment();
	CDevelopment(const TInitInfo&);
	~CDevelopment();
	
	// general functions //	
	void setParms();
	void update(int iCur, const TWeather&, const TInitInfo, const TNitrogenStatus, double area, double greenLAI, double maxLAI, double dayinc, double Tdaymax, double Tdaymin, double Tdayave, double Tdaylag, double Sradave, double Photoave, double Creserve, double Cpool);  // updates whole plant

	// member functions for thermal responses//
	double beta_fn(double t, double R_max, double t_m, double t_e);
	//double calcGTI(double, bool);  //not used?
	//double calcGDD(double); // not used?
	void response(int, const TInitInfo, const TWeather&, double, double, double); //whole plant thermal responses
	void tubinit(int, const TInitInfo, const TWeather&, const TNitrogenStatus, double, double, double, double, double, double); // compute plant tuber initiation response
	double fdtt17(double, double, double, double); //thermal time base 17 (SIMGUI)
	double fdtt22(double, double, double, double); //thermal time base 22 (SIMGUI)
	double fstt(double, double, double, double); // thermal time for soil (SIMGUI)
	
	// functions to get private data//
	bool Germinated() {return germination.done;}
	bool Emerged() {return emergence.done;}
	double dEmerged() {return emergence.day;}
	bool Vegetative() {return vegetative.done;}
	double dVegetative() {return vegetative.day;}
	bool Flowered() {return flowered.done;}
	double dFlowered() {return flowered.day;}
	virtual void set_Flowered(double x) {flowered.done = true;flowered.day = x;}

	bool Tuberlinear() {return tuberlinear.done;}
	double dTuberlinear() {return tuberlinear.day;}
	bool Tuberdominant() {return tuberdominant.done;}
	bool Matured() {return maturity.done;}
	TInitInfo get_initInfo() {return initInfo;}
	int get_finalLeafNo() {return finalLeafNo;}
	int get_istage() {return iStage;}
	double get_xstage() {return dblXstage;}
	double get_DTT17() {return dblDtt17;}
	double get_CDD17() {return dblCdtt17;}
	double get_DTT22() {return dblDtt22;}
	double get_CDTT22() {return dblCdtt22;}
	double get_tind() {return dblTind;}
	double get_partub() {return dblPartub;}
	double get_phyllochronsFromTI() {return phyllochronsFromTI;}
	double get_LvsAtTI() {return LvsAtTI;}
	double get_LvsInitiated(){return LvsInitiated;}
	double get_LvsAppeared(){return LvsAppeared;} // tip appreance rate is most conservative across cultivars and temperature regimes, persoanl comm with Dr. Tollenaar
	double get_LAR() {return dblLAR;}
	double get_BAR() {return dblBAR;}
	double get_Tamp_factor() {return ampfac;}
	double get_Area_factor() {return areafac;}
	double get_Daylength_factor() {return dlfac;}
	double get_Nitrogen_factor() {return nitfac;}
	double get_Radiation_factor() {return solfac;}
	double get_Tave_factor() {return tempfac;}
	double get_TuberIndex() {return tday;}


	virtual void set_LvsAppeared(int x){LvsAppeared += x;}

private:
	CDevelopment(const CDevelopment&); // supressing shallow copy constructor
	
	// thermal time and temperature response information
	double Rmax_leaf = 0; //nonlinear beta function for leaf (nodal unit) appearance
	double Rmax_stem = 0; //nonlinear beta function for basal branch appearance
	double T_base = 0, T_opt = 0, T_ceil = 0; //nonlinear beta function
	double dblDtt17 = 0, dblDtt22 = 0, dblStt = 0, dblCstt = 0, dblCdtt17 = 0, dblCdtt22 = 0; // T response variables from SIMGUI
	double dblLAR = 0; //leaf appearance rate per time increment
	double dblBAR = 0; //basal appearance rate per time increment

	// developmental stage tracking variables
	int iDAS = 0; //das after sowing
	double dblDAS = 0;
	int iStage = 0; //tracks istage as in SIMGUI
	int iPhase = 0; // keeps track of phenology during transition stages, equals 1 if there is a phase change
	int iEmerge = 0; // 0 or 1 to indicate emergence from soil
	int iJinit = 0; // julian day at which tuber initiation begins
	int iDflag = 0; // 0 or 1 to indicate appearance of apical stems
	int iMatjd = 0; // julian day of maturity
	int iRet = 0; // equals 1 when plant is at maturity
	double dblXstage = 0; //transition within given stage
	double dblFlowering = 0; // date (dae) of flower appearance
	TInitInfo initInfo;
	//TNitrogenStatus NStatus;
	TEvent preplanting;
	TEvent germination;
	TEvent emergence;
	TEvent vegetative;
	TEvent flowered;
	TEvent tuberlinear;
	TEvent tuberdominant;	
	TEvent maturity;
	
	// pre-emergence
	double dblSprlth = 0; //sprout length prior to emergence from soil, [cm]
	double dblSprlthp = 0; //sprout length at planting, [cm]
	double dblSdepth = 0; //current depth in soil of sprout?, [cm]
	double dblTsprwt = 0; //total sprout weight, [g plant]

	// Tuber variables
	double dblTind = 0, dblPartub = 0; // tuber induction factor, tuber partitioning factor

	double totLeafNo = 0, juvLeafNo = 0, LvsAtTI = 0, phyllochronsFromTI = 0; //number of total, juvenile (genetic coeff) lvs, and lvs appeared at tassel initiation
	double P2 = 0; //photoperiod sensitivity as used in CERES-Maize
	double GerminationRate = 0, EmergenceRate = 0, LvsInitiated = 0, LvsAppeared = 0, LvsExpanded = 0, Anthesis = 0;
	int initLeafNo = 0;//initial number of leaves at emergence
	int finalLeafNo = 0, curLeafNo = 0;  //do these need to be used?
	
	double ampfac = 0, areafac = 0, dlfac = 0, nitfac = 0, solfac = 0, tday = 0, tmn = 0, tempfac = 0;

};
#endif