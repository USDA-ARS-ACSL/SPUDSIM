// Updated Gas Exchange Class - Soo Kim, Updated Dennis Timlin

/*!  Coupled model of photosynthesis-stomatal conductance-energy balance for a maize leaf this unit simulates Maize leaf gas-exchange characteristics
  including photosynthesis, \n traspiration, boundary and stomatal conductances,
  and leaf temperature based \n on von Caemmerer (2000) C4 model, BWB stomatal
conductance (1987) and \n Energy balance model as described in Campbell and Norman (1998) 

Photosynthetic parameters were calibrated with PI3733 from
SPAR experiments at Beltsville, MD in 2002.

For potato, parameters come from Potato data in indoor chambers at Beltsville, MD in 2003
Stomatal conductance parameters were not calibrated

Nitrogenstress reduces Jmax and Vcmax in same proportion and is based on leaf N deficit

@authors Soo-Hyung Kim, Univ of Washington \n Dennis Timlin, USDA-ARS \n  David Fleisher, USDA-ARS \n
@version 1.0
@date August 2013

@note <b>-Bibliography </b>\n
Kim, S.-H., and J.H. Lieth. 2003. A coupled model of photosynthesis, stomatal conductance and transpiration for a rose leaf (Rosa hybrida L.). Ann. Bot. 91:771–781. \n
Kim, S.-H., D.C. Gitz, R.C. Sicher, J.T. Baker, D.J. Timlin, and V.R. Reddy. 2007. Temperature dependence of growth, development, and photosynthesis in maize under elevated CO2. Env. Exp. Bot. 61:224-236. \n
Kim, S.-H., R.C. Sicher, H. Bae, D.C. Gitz, J.T. Baker, D.J. Timlin, and V.R. Reddy. 2006. Canopy photosynthesis, evapotranspiration, leaf nitrogen, and transcription  \n
*/
// ws

#include "stdafx.h"
#include "gas_exchange_new.h"
#include <cmath>
using namespace std;

#define R 8.314 //ideal gas constant
#define maxiter 200 //maximum number of iterations
#define epsilon 0.97 //emissivity (Cambpell and Norman, 1998, pg 163
#define sbc 5.6697e-8 //stefan-boltzmann constant W m-2 k-4 - actually varies somewhat with temperature
#define scatt 0.15 //leaf reflectance + transmittance
#define f 0.15 //spectral correction
#define O 205.0 //gas units are mbar
#define Q10 2.0 //Q10 factor

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

inline double Square(double a) {return a*a;}
inline double Min(double a, double b, double c) {return (__min(__min(a,b),c));}

CGas_exchange_new::CGas_exchange_new()
{
	isCiConverged = false;
	errTolerance = 0.001;
	eqlTolerance = 1.0e-6;
}

CGas_exchange_new::~CGas_exchange_new()
{
}

void CGas_exchange_new::getParms(double Tubmod, double Feedback, double Nitrogendeficiencyone)
{
	/*********************/
	/* Assigns parameters to the 3 components of Farquar C3 model
	/* Includes T-dependencies
	/*********************/
	//Feedback = 1;
	double x = 0.;
	// These are C3 parameters from default potato information
	//  They include stepwise calculation method and BBW model for Gs
	if (Tubmod > 1) 
	{
		x = 1.3; //assume Pmax is doubled after TI
	}
	else x = 1.0;
	//x = 1.3;  //override increase in Pmax following TI
	x = 1.0;
	x = Tubmod;
	/**the following were derived from modified SPUDSIM RUE-T response functions, not published, DHF */
	/* THese are too low....
	Parameters.Vcm25 = 67; //*the following were derived from modified T response functions, not published, DHF 
	Parameters.Eav = 31100; //the Vcm25 was obtained from assuming Vcm20 = 72 as optimal
	Parameters.Hv = 149500;
	Parameters.Sv = 498;

	Parameters.Jm25 = 133; //the following were dervied from modified T response functions, not pulbished
	Parameters.Eaj = 31100; //The Jm of 133 was obtained from assuming Jm20 = 142 as optimal
	Parameters.Hj = 149500;
	Parameters.Sj = 498;
	*/
	//Previous Parameters Used in Older Edition - the T responses here were not correct for potato
	Parameters.Eav = 64800.;
	Parameters.Eaj = 37000.;
	Parameters.Hj = 220000.;
	Parameters.Sj = 710.;
	//Parameters.Vcm25 = 72; // default value, rubsico capacity for co2 at 25C, umol m-2 s-1
	//Parameters.Vcm25 = 86.5; //from Lawson et al, 2001 data
	//Parameters.Vcm25 = x*102.4; //from Soo's Rose paper
	Parameters.Vcm25 = Feedback*Nitrogendeficiencyone*x*108.4; //from Dennis data, 2002: Use as of 5/ 2008
	//Parameters.Jm25 = 142; //default value, pot rate of e transport, umol electrons m-2 s-1
	//Parameters.Jm25 = 188.8; //reduce rate of e transport in potato - need data or reference?
	//Parameters.Jm25 = x*162; //from Soo's rose paper
	Parameters.Jm25 = Feedback * Nitrogendeficiencyone *x*185.6; //from Dennis data, 2002: Use as of 5/ 2008
	//Parameters.Jm25 = x*142; //Matches Ng and Loomis and was default...
	//Parameters.TPU25 = x*13.5; //from Lawson et al.,. 2001
	Parameters.TPU25 = Feedback*x*11.55;//from soo's rose paper, rate of TP use umol m-2 s-1
	Parameters.Rd25 = 1.;// mitochondrial resp in the light at 25C, umol m-2 s-1
	Parameters.Eap = 47100.;
	//Parameters.g0 = 0.07; //0.096 from Kim et al. was original value; 0.07 from Liu et al., 2009 - measured data, mol m-2 s-1
	Parameters.g0 = 0.01; //Maizsim value and Vos and Oyarzum (reduces to 0!)
	//Parameters.g0 = 0.096;
	//Parameters.g1 = 16.57; //10.055 from Kim et al., was original value; 16.57 from Liu et al., 2009 - " with potato, unitless
	//Parameters.g1 = 10.055;
	Parameters.g1 = 12.0 ; //from Vos and Oyarzun, 1987 data

	Parameters.beta_ABA = 1.48e2; //Tardieu-Davies beta, Dewar (2002)
	Parameters.delta = -1.0;
	Parameters.a_ABA = 1.0e-4;
	Parameters.lambda_r = 4.0e-12; //Dewar
	Parameters.lambda_l = 1.0e-12;
	Parameters.K_max = 6.67e-3; // max. xylem conductance (mol m-2 s-1 MPa-1) from root to leaf, Dewar (2002)
	Parameters.d1 = 0.1468; Parameters.d2= 0.0103; // leaf age dependence parameters
}

void CGas_exchange_new::SetVal(double PFD, const TInitInfo info, double Nitrogendeficiencyone, double Tair, double CO2, double RH, double wind,  double Press, double width, double Tubmod, double Feedback, double leafage, double psil, double etsupply)
{
	/***********************************/
	/* Sets initial parameters and runs routines
	/* for the coupled biochemical C3 gas exchange, stomatal conductance, energy balance routine
	/* orignally coded by Soo Kim
	/***********************************/
	//const double	scatt = 0.15;
	this->PFD = PFD;
    double PAR = (PFD/4.55); //w par m-2
    double NIR = PAR; // If total solar radiation unavailable, assume NIR the same energy as PAR waveband
    this->R_abs = (1-scatt)*PAR + 0.15*NIR + 2*(epsilon*sbc*pow(Tair+273.,4.)); // times 2 for projected area basis
	// shortwave radiation (PAR (=0.85) + NIR (=0.15) solar radiation absorptivity of leaves: =~ 0.5
    this->CO2 = CO2;
    this->RH = __min(100.0, __max(RH, 25.0))/100;
    this->Tair = Tair;
    this->width = width;
	this->wind = wind;
    this->Press = Press;
	this->psileaf = psil;
	this->leaf_age = leafage;
    getParms(Tubmod, Feedback, Nitrogendeficiencyone);
	GasEx_psil(psileaf, etsupply, info);
}

void CGas_exchange_new::GasEx_psil(double psileaf, double etsupply, const TInitInfo info)
{
/**********************/
	/* Main looping routine for coupled models
	/*	Incorporate water stress effect
	/***********************/
	double Tleaf_old =0.;
	int iter = 1;
	iter_total = 0;
	Tleaf = Tair; 
	Tleaf_old = 0.;
	Ci = 0.7*CO2;
	gb = gbw(); //bounday layer conductance
	gs = gsw(psileaf, info); //stomatal conductance
	if (Tair < 1) //DHF - Assume enzymatic activity negligable at this point, so no photosynthesis.  Need literature?  At what point does plant die?
		// Loomis and Conner indicate C3 leaves can go to 0C without permanent damage, not sure what temporary affect is on gas exchange rates is.
	{
		A_net = -Rd();
		A_gross = __max(Rd(),0);
		gs = gsw(psileaf,info);
		EnergyBalance(etsupply); //assuming stomatal conductance, ET unaffected by very low temps
		return;
	}
	while ((fabs(Tleaf_old-Tleaf)>0.01) && (iter < maxiter))
	{
		Tleaf_old = Tleaf;
		Ci = SearchCi(Ci,psileaf,info);
		gs = gsw(psileaf,info);
		EnergyBalance(etsupply);
		iter2 =++ iter;	
	}

}

double CGas_exchange_new::gbw(void) //boundary layer conductance to vapor
{
	const double stomaRatio = 0.5; // for potato.
	double ratio;
	double d;
	ratio = Square(stomaRatio+1)/(Square(stomaRatio)+1);
    d = width/100*0.72; // characteristic dimension of a leaf, leaf width is converted from cm to m
  //  return 1.42; // total BLC (both sides) for LI6400 leaf chamber
    return (1.4*0.147*sqrt(__max(0.1,wind)/d))*ratio; 
	// multiply by 1.4 for outdoor condition, Campbell and Norman (1998), p109
	// multiply by ratio to get the effective blc (per projected area basis), licor 6400 manual p 1-9
}

double CGas_exchange_new::gsw(double pressure, const TInitInfo info)  // stomatal conductance for water vapor in mol m-2 s-1 
{
	double  Pn, aa, bb, cc, Ha, Hs, Cs, Ca, gg, gamma, Ds;
	double temp = set_PSIleafeffect(pressure, info);
    gsModel    myModel=BBW;     // need to put this in main program 
    Ca = CO2;
    Ha = RH;
    gamma = 36.9 + 1.88*(Tleaf-25)+0.036*Square(Tleaf-25);  // CO2 compensation point
	//gamma = 10;
	double P= Press /100;
	Cs = (Ca - (1.37*A_net/gb))*P; // surface CO2 in mole fraction
	if (Cs==gamma) Cs = gamma + 1;
	if (Cs<=gamma) Cs = gamma + 1;
	
	//if (f_age == 0) gg = 0.001; 
    //  else gg = Parameters.g0*f_age; // prevent division by zero
	gg = Parameters.g0;
	if (A_net <= 0) Pn = 0.00001;
       else Pn = A_net; // solving quadratic formula for surface humidity(hs)
    aa = temp*Parameters.g1*A_net/Cs;
    bb = gg+gb-(temp*Parameters.g1*Pn/Cs);
    cc = (-Ha*gb)-gg;
	Hs = QuadSolnUpper(aa,bb,cc);       // RH at leaf surface
	if (Hs > 1) Hs = 1.;
	if (Hs < 0) Hs = 0.;
    Ds = (1-Hs)*Es(Tleaf); // VPD at leaf surface
	if (A_net < 0) {
		return gg;
	}
	else{
		switch (myModel)
		{
		    case BBW: 
				return (gg+temp*Parameters.g1*(A_net*Hs/Cs));
				break;
			case L: 
				return gg + Parameters.g1*A_net/((Cs-gamma)*(1+Ds/Parameters.g2));
				break;
			case H:
				return (gg+Parameters.g1*(A_net*Ha/Ca));
				break;
			default:
				return (gg+temp*Parameters.g1*(A_net*Hs/Cs));
		}
	}
}

double CGas_exchange_new::set_PSIleafeffect(double pressure, const TInitInfo info)
{
	//Reduction in stomatal conductance using hourly bulk leaf water potential in MPa
	double sf, phyf;
	//sf = 2.3; phyf = -1.2;
	//sf = 6.1; phyf = -0.67; //sensitivity parameter Tuzet et al., (2003), tuned for potato using Vos and Oyarzun (1987) by DHFreference potential Tuzet et al., (2003), " by DHF
	sf = 4.0; phyf = -0.4; //as above, but using the Liu et al., 2009 parameters for potato for g0 and g1
	//if (pressure < -0.05) 
	//if (pressure < phyf)
	if (pressure < -0.05)
		this->psileaf_stress = __max((1+exp(sf*phyf))/(1+exp(sf*(phyf-pressure))),0.);
	else psileaf_stress = 1;
	//psileaf_stress = 1;
	double temp = this->psileaf_stress;
	if (info.Water_stress_off == 1) temp = 1; //remove waterstress if stress option is turned off
	if (info.Water_stress_simulation_type == 1) temp = 1;//don't use if SIMPOTATO approach employed

	return temp;
}

double CGas_exchange_new::SearchCi(double CO2i, double psileaf, const TInitInfo info)
{
	/*secant search to find optimal internal CO2 concentration*/
	int iter = 0;
	double fprime = 0., Ci1 = 0., Ci2 = 0., Ci_low = 0., Ci_hi = 0., Ci_m = 0.;
	double temp = 0.;
	Ci1 = CO2i;
	Ci2 = CO2i + 1.0;
	Ci_m = (Ci1 + Ci2) / 2.0;
	iter_Ci = 0;
	iter = 0;
	isCiConverged = true;
	do
	{
		iter++;
		//secant search method
		if (abs(Ci1-Ci2) <= errTolerance) {break;}
		if (iter >= maxiter)
		{
			isCiConverged = false;
			break;
		}
		fprime = (EvalCi(Ci2,psileaf,info) - EvalCi(Ci1,psileaf,info))/(Ci2-Ci1); //f'(Ci)
		if (fprime != 0.0)
		{
			Ci_m = max(errTolerance, Ci1-EvalCi(Ci1,psileaf,info)/fprime);
		}
		else
			Ci_m = Ci1;
		Ci1 = Ci2;
		Ci2 = Ci_m;
		temp = EvalCi(Ci_m,psileaf,info);
		double temp2 = maxiter;
	} while ((abs(EvalCi(Ci_m,psileaf,info)) >= errTolerance) || (iter < maxiter));

	//use bisectional search if above doesn't converge
	if (iter > maxiter)
	{
		Ci_low = 0.0;
		Ci_hi = 2.0*CO2;
		isCiConverged = false;
		while (abs(Ci_hi - Ci_low) <= errTolerance || iter > (maxiter*2.))
		{
			Ci_m = (Ci_low + Ci_hi)/2.;
			if (abs(EvalCi(Ci_low,psileaf,info)*EvalCi(Ci_m,psileaf,info)) <= eqlTolerance) break;
			else if (EvalCi(Ci_low,psileaf,info)*EvalCi(Ci_m,psileaf,info) < 0.0) Ci_hi = max(Ci_m,errTolerance);
			else if (EvalCi(Ci_m,psileaf,info)*EvalCi(Ci_hi,psileaf,info) < 0.0) Ci_low = max(Ci_m,errTolerance);
			else 
			{
				isCiConverged = false; break;
			}
		}
	}
	CO2i = Ci_m;
	Ci_Ca = CO2i / CO2;
	iter_Ci = iter_Ci + iter;
	iter_total = iter_total + iter;
	return CO2i;
}


double CGas_exchange_new::EvalCi(double Ci, double psileaf, const TInitInfo info)
{
	//Calculates a new value for Ci for the current values of photosynthesis and stomatal conductance
	//Determined from parameters from prior stem where energy balance was solved
	double newCi;
	Photosynthesis(Ci,psileaf,info);
	if(abs(gs) > eqlTolerance)
	{
		newCi = max(1.0,CO2-A_net * (1.6/gs + 1.37/gb)*(Press/100.0));
	}
	else
		newCi = max(1.0,CO2 - A_net * (1.6/eqlTolerance + 1.37/gb)*(Press/100.0));
	return (newCi-Ci);
}

void CGas_exchange_new::Photosynthesis(double Ci, double psileaf, const TInitInfo info)
{
	//C3 photosynthesis
	const double curvature = 0.999; // curvature factor of Av and Aj collimitation
	const double    theta = 0.7;		//for potato
	const int Kc25 = 404;
	const int Ko25 = 248;
	const double Eac = 59400.;
	const double Eao = 36000.;
	double alpha, Kc, Ko, gamma, Ia, Jmax, Vcmax, TPU, J, Av, Aj, Ap, Ac, Km, Ca, Cc, P,Tk;
	double f_age;
	f_age = this->leaf_age;
	f_age = 1;
	Tk = Tleaf + 273.0;
	gamma = 36.9 + 1.88*(Tleaf-25)+0.036*Square(Tleaf-25); //CO2 comp point in absense of mito respiration, in ubar
	Ia = PFD*(1-scatt);//absorbed irradiance
	alpha = (1-f)/2.; //apparent quantum efficiency, params adjusted to get value 0.3 for average C3 leaf
	A_net = 0;
	P  = Press/100.;
	Ca = CO2*P; //* conversion to partial pressure */ 
	Kc = Kc25*exp(Eac*(Tleaf-25.)/(298.*R*(Tleaf+273.)));
	Ko = Ko25*exp(Eao*(Tleaf-25.)/(298.*R*(Tleaf+273.)));
	Km = Kc*(1+O/Ko); //* effective M-M constant for Kc in the presence of O2 */
    Jmax = f_age*Parameters.Jm25*exp(((Tk-298.)*Parameters.Eaj)/(R*Tk*298.))*
		     (1+exp((Parameters.Sj*298.-Parameters.Hj)/(R*298.)))/
			 (1+exp((Parameters.Sj*Tk-Parameters.Hj)/(R*Tk))); // de Pury 1997
    //Vcmax = f_age*Parameters.Vcm25*exp(Parameters.Eav*(Tleaf-25)/(298*R*(Tleaf+273))); //old response from corn model
    Vcmax = f_age*Parameters.Vcm25*exp(((Tk-298.)*Parameters.Eav)/(R*Tk*298.))*
		     (1+exp((Parameters.Sv*298.-Parameters.Hv)/(R*298.)))/
			 (1+exp((Parameters.Sv*Tk-Parameters.Hv)/(R*Tk))); // Used peaked response, DHF
    TPU = f_age*Parameters.TPU25*exp(Parameters.Eap*(Tleaf-25.)/(298.*R*(Tleaf+273.)));
	Cc = Ci; //assume infinite gi
	gs = gsw(psileaf, info);
	gb = gbw();
	Av = (Vcmax*(Cc-gamma))/(Cc+Km);
    J =  (((alpha*Ia + Jmax) - sqrt(Square(alpha*Ia+Jmax) - 4.*alpha*Ia*(Jmax)*theta)) / (2.*theta)) ;
    Aj = J*(Cc-gamma)/(4.*(Cc+2.*gamma));
    Ap = 3.*TPU;
	Ac = ((Av+Aj)-sqrt(Square(Av+Aj)-4*curvature*Av*Aj))/(2.*curvature); //curvature account for collimitation between Av and Aj
	if (Cc > gamma)
		A_net = min(Ac,Ap)-Rd();
	else
	{
		A_net = Av -Rd();
	}
	A_gross = max(A_net+Rd(),0.0);
	gs = gsw(psileaf, info);
}

void CGas_exchange_new::EnergyBalance(double Jw)
// see Campbell and Norman (1998) pp 224-225
// because Stefan-Boltzman constant is for unit surface area by denifition,
// all terms including sbc are multilplied by 2 (i.e., gr, thermal radiation)
{
    const double lambda = 44000.; //J mol-1 at 25C
    const double psc = 6.66e-4;
    const double Cp = 29.3; // thermodynamic psychrometer constant and specific hear of air
	double gha = 0., gv = 0., gr = 0., ghr = 0., psc1 = 0., Ea = 0., thermal_air = 0., Ti = 0., Ta = 0.;
	double lastTi = 0., newTi = 0.;
	int iter = 0;

    Ta = Tair;
    Ti = Tleaf;
    //gha = gb*(0.135/0.147);  // heat conductance, gha = 1.4*.135*sqrt(u/d), u is the wind speed in m/s} Note: this was only true if stomatal ratio = 1
	gha = 1.4 * 0.135 * sqrt(__max(0.1,wind)/(width/100*0.72));
    gv = gs*gb/(gs+gb); 
    gr = (4.*epsilon*sbc*pow(273.+Ta,3.)/Cp)*2.; // radiative conductance, 2 account for both sides
    ghr = gha + gr;
    thermal_air = epsilon*sbc*pow(Ta+273.,4.)*2.; // emitted thermal radiation
    psc1 = psc*ghr/gv; // apparent psychrometer constant
    VPD = Es(Ta)*(1.-RH); // vapor pressure deficit
    Ea = Es(Ta)*RH; // ambient vapor pressure

	//iterative version
	newTi =-10;
	iter = 0;
	lastTi = Tleaf;
	double Res = 0., dRes = 0.;
	double thermal_leaf = 0.;
	while ((abs(lastTi-newTi)>0.001) &&(iter < maxiter))
	{
		lastTi = newTi;
		Tleaf = Ta + (R_abs - thermal_air - lambda*gv*VPD/Press)/(Cp*ghr + lambda * Slope(Tair)*gv);
		thermal_leaf = epsilon*sbc*pow(Tleaf+273.,4.)*2.;
		Res = R_abs - thermal_leaf - Cp*gha*(Tleaf-Ta)-lambda * gv*0.5*(Es(Tleaf)-Ea)/Press;
		dRes = -4.*epsilon*sbc*pow(273.+Tleaf,3.)*2. - Cp*gha*Tleaf - lambda*gv*Slope(Tleaf);
		newTi = Tleaf + Res/dRes;
		iter++;
	}
	Tleaf = newTi;
    ET = __max(0,1000.*gv*((Es(Tleaf)-Ea)/Press)/(1-(Es(Tleaf)+Ea)/(Press)));//1000 is to go from mol to mmol
    // accounting for additional transp. because of mass flow, see von Caemmerer and Farquhar (1981)
}



double CGas_exchange_new::Es (double T) //Campbell and Norman (1998), p 41 Saturation vapor pressure in kPa
{
    return (0.611*exp(17.502*T/(240.97+T)));
}

double CGas_exchange_new::Rd()   //Should be expanded to include other env. and physiological factors
{
	const double Ear =  66400.; //exponential rate of arrhenious function for mitochondrial respiration (J mol)
	return (Parameters.Rd25*exp(Ear*(Tleaf-25.)/(298.*R*(Tleaf+273.)))); 
	//return 0;
}

double CGas_exchange_new::Slope(double T) // slope of the sat vapor pressure curve: first order derivative of Es with respect to T
{
   const double b= 17.502; const double c= 240.97;
   return (Es(T)*(b*c)/Square(c+T)/Press); 
}


double QuadSolnUpper (double a, double b, double c )
{
   if (a==0) return 0;
	else if ((b*b - 4.*a*c) < 0) return -b/a;   //imaginary roots
    else  return (-b+sqrt(b*b-4.*a*c))/(2.*a);
}

double QuadSolnLower (double a, double b, double c )
{
    if (a==0) return 0;
	else if ((b*b - 4*a*c) < 0) return -b/a;   //imaginary roots
    else  return (-b-sqrt(b*b-4.*a*c))/(2.*a);
}
