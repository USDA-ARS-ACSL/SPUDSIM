#pragma once
#include "initInfo.h"
double QuadSolnUpper (double a, double b, double c );
double QuadSolnLower (double a, double b, double c );

class CGas_exchange
{
public:
	CGas_exchange(void);
	~CGas_exchange(void);

private:
	TInitInfo initInfo;
	double PFD, R_abs, Tair, CO2, RH, wind, age, SLA, width, Press, N; 
	double psileaf; //leaf water potential, MPa
	double psileaf_stress; //0 to 1 factor for stomatal closure effect of leaf water potential
	double leaf_age;
	void GasEx();
	void GasEx_psil(double psileaf, double etsupply, const TInitInfo info);
	void Photosynthesis(double pressure, const TInitInfo);
	void EnergyBalance(double pressure);
	void getParms(double Tubmod, double Feedback, double LAIclass);
	double gsw(double pressure, const TInitInfo info);
	double gbw();
	double Es(double T);
	double Slope(double T);
	double Rd();
	double set_PSIleafeffect(double pressure, const TInitInfo info);

public:
  void SetVal(double PFD, const TInitInfo info, double Tair, double CO2, double RH, double wind, double Press, double width, double Tubmod, double Feedback, double leafage, double psil, double etsupply, double laiclass);
  double get_VPD() {return VPD;}
  double get_gs() {return gs;}
  double get_ci() {return Ci;}
  struct tparms
  {
	  double Eav, TPU25, Eap, Eaj, Jm25, Vcm25, Rd25, Sj, Hj, Sv, Hv, g0, g1, g2, d1, d2, beta_ABA, delta, a_ABA, lambda_r, lambda_l, K_max;
	 // double Vpm25, Vcm25, Jm25, Rd25, EaVp, EaVc, Eaj, Sj, Hj, Ear, g0, g1, g2;
  } Parameters;
	enum gsModel {BBW, L, H};
	enum CalMethod {Stepwise, Simultaneous};
  double A_gross, A_net, ET, Tleaf, Ci, gs, gb, Rdc, iter1, iter2, VPD;
  
};
