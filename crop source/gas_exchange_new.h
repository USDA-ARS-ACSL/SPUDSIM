#pragma once
#include "initInfo.h"
double QuadSolnUpper (double a, double b, double c );
double QuadSolnLower (double a, double b, double c );

/*! \class CGasExchange
* \brief Class for gas exchange calculations\n
* This class simulates gas exchange in plant leaves for C3 and C4 plants. \n
* \par Usage
	- Use <b>SetParams</b> to initialize the GasExchange object with parameters for a specific variety of plant
	- Use <b>SetVal</b> to pass environmental variables for a simulation and return a structure with output. \n
* See \ref Interface for details on input and output


	*/

class CGas_exchange_new
{
public:
	CGas_exchange_new(void);
	~CGas_exchange_new(void);

private:
	TInitInfo initInfo;
	double PFD = 0, R_abs = 0, Tair = 0, CO2 = 0, RH = 0, wind = 0, age = 0, SLA = 0, width = 0, Press = 0, N = 0;
	double psileaf = 0; //leaf water potential, MPa
	double psileaf_stress = 0; //0 to 1 factor for stomatal closure effect of leaf water potential
	double leaf_age = 0;
	void GasEx_psil(double psileaf, double etsupply, const TInitInfo info);
	void Photosynthesis(double Ci, double psileaf, const TInitInfo);
	void EnergyBalance(double pressure);
	void getParms(double Tubmod, double Feedback, double Nitrogendeficiencyone);
	double SearchCi(double CO2i, double psileaf, const TInitInfo info);
	double EvalCi(double Ci, double psileaf, const TInitInfo info);

	double gsw(double pressure, const TInitInfo info);
	double gbw();
	double Es(double T);
	double Slope(double T);
	double Rd();
	double set_PSIleafeffect(double pressure, const TInitInfo info);
	double Ci_Ca = 0;  //!< ratio of internal to external CO2, unitless
	double errTolerance = 0; /*!< error tolerance for iterations */
	double eqlTolerance = 0; /*!< equality tolerance */
	int iter_total = 0;      //!< holds total number of iterations */
	int iter1 = 0, iter2 = 0;         //!< holds iteration counters
	int  iter_Ci = 0;   /*!< iteration value for Ci */
	bool isCiConverged = false; /*!< true if Ci iterations have converged */

public:
	void SetVal(double PFD, const TInitInfo info, double Nitrogendeficiencyone, double Tair, double CO2, double RH, double wind, double Press, double width, 
		double Tubmod, double Feedback, double leafage, double psil, double etsupply);
	double get_VPD() {return VPD;}
	double get_gs() {return gs;}
	double get_ci() {return Ci;}
	double get_psileaf_Pn_stress() {return psileaf_stress;}

	struct tparms
	{
		  double Eav, TPU25, Eap, Eaj, Jm25, Vcm25, Rd25, Sj, Hj, Sv, Hv, g0, g1, g2, d1, d2, beta_ABA, delta, a_ABA, lambda_r, lambda_l, K_max;
	} Parameters;
	enum gsModel {BBW, L, H};
	enum CalMethod {Stepwise, Simultaneous};
	double A_gross = 0, A_net = 0, ET = 0, Tleaf = 0, Ci = 0, gs = 0, gb = 0, Rdc = 0, VPD = 0;
  
};


	
	

 