//Coupled model of photosynthesis-stomatal conductance-energy balance for a maize leaf
// this unit simulates Maize leaf gas-exchange characteristics
// including photosynthesis, traspiration, boundary and stomatal conductances,
// and leaf temperature based on von Caemmerer (2000) C4 model, BWB stomatal
// conductance (1987) and Energy balance model as described in Campbell and
// Norman (1998) photosynthetic parameters were calibrated with PI3733 from
// SPAR experiments at Beltsville, MD in 2002 Stomatal conductance parameters
// were not calibrated
// Ver 1.0, S.Kim, 11/2002, Originally written in Pascal
// Translated into C++ by S.Kim June 11, 2003 following D.Timlin's translation of C3 model
// IMPORTANT: This model must not be released until validated and published


//DHF: This is now only the C3 routine mentioned above (i.e., C4 components deleted)

#include "stdafx.h"
#include "gas exchange.h"
#include  "math.h"
#include <stdlib.h>

#define R 8.314  // idealgasconstant
#define maxiter 100
#define epsilon 0.97
#define sbc 5.6697e-8
#ifdef _DEBUG
#define new DEBUG_NEW
#endif
inline double Square(double a) { return a * a; }
inline double Min(double a, double b, double c) {return (__min(__min(a,b),c));}

CGas_exchange::CGas_exchange()
{ 
}

CGas_exchange::~CGas_exchange()
{
}

void CGas_exchange::getParms(double Tubmod, double Feedback, double LAIclass)
{
	/*********************/
	/* Assigns parameters to the 3 components of Farquar C3 model
	/* Includes T-dependencies
	/*********************/
	//Feedback = 1;
	double x = 0;
	// These are C3 parameters from default potato information
	//  They include stepwise calculation method and BBW model for Gs
	if (Tubmod > 1) 
	{
		x = 1.3; //assume Pmax is doubled after TI
	}
	else x = 1.0;
	//x = 1.3;  //override increase in Pmax following TI
	x = 1.;
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
	Parameters.Eav = 64800;
	Parameters.Eaj = 37000;
	Parameters.Hj = 220000;
	Parameters.Sj = 710;
	//Parameters.Vcm25 = 72; // default value, rubsico capacity for co2 at 25C, umol m-2 s-1
	//Parameters.Vcm25 = 86.5; //from Lawson et al, 2001 data
	//Parameters.Vcm25 = x*102.4; //from Soo's Rose paper
	Parameters.Vcm25 = Feedback*x*108.4; //from Dennis data, 2002: Use as of 5/ 2008
	//Parameters.Jm25 = 142; //default value, pot rate of e transport, umol electrons m-2 s-1
	//Parameters.Jm25 = 188.8; //reduce rate of e transport in potato - need data or reference?
	//Parameters.Jm25 = x*162; //from Soo's rose paper
	Parameters.Jm25 = Feedback*x*185.6; //from Dennis data, 2002: Use as of 5/ 2008
	//Parameters.Jm25 = x*142; //Matches Ng and Loomis and was default...
	//Parameters.TPU25 = x*13.5; //from Lawson et al.,. 2001
	Parameters.TPU25 = Feedback*x*11.55;//from soo's rose paper, rate of TP use umol m-2 s-1
	Parameters.Rd25 = 1;// mitochondrial resp in the light at 25C, umol m-2 s-1
	Parameters.Eap = 47100;
	Parameters.g0 = 0.07; //0.096 from Kim et al. was original value; 0.07 from Liu et al., 2009 - measured data, mol m-2 s-1
	//Parameters.g0 = __min(0.07,0.07/LAIclass);
	//Parameters.g0 = __min(0.07/abs((LAIclass-2)),0.07);
	//if (LAIclass > 3) Parameters.g0 = 0.07 - 0.07*0.1*(LAIclass-3);
	
	Parameters.g0 = 0.01; //Maizsim value
	Parameters.g1 = 16.57; //10.055 from Kim et al., was original value; 16.57 from Liu et al., 2009 - " with potato, unitless
	Parameters.beta_ABA = 1.48e2; //Tardieu-Davies beta, Dewar (2002)
	Parameters.delta = -1.0;
	Parameters.a_ABA = 1.0e-4;
	Parameters.lambda_r = 4.0e-12; //Dewar
	Parameters.lambda_l = 1.0e-12;
	Parameters.K_max = 6.67e-3; // max. xylem conductance (mol m-2 s-1 MPa-1) from root to leaf, Dewar (2002)


	Parameters.d1 = 0.1468; Parameters.d2= 0.0103; // leaf age dependence parameters
}

void CGas_exchange::GasEx_psil(double psileaf, double etsupply, const TInitInfo info)
{
	/**********************/
	/* Main looping routine for coupled models
	/*	Incorporate water stress effect
	/***********************/
	double Tleaf_old;
	int iter = 1;
	Tleaf = Tair; 
	Tleaf_old = 0;
	Ci = 0.7*CO2;
	gb = gbw();
	gs = gsw(psileaf, info); //gs = gsw()
	A_net = (CO2-Ci)/((1.57/gs+1.37/gb)*Press/100);
	if (Tair < 4) //DHF - Assume enzymatic activity negligable at this point, so no photosynthesis.  Need literature?  At what point does plant die?
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
		Photosynthesis(psileaf,info);
		EnergyBalance(etsupply);
		iter2 =++ iter;	
	}
} 

void CGas_exchange::Photosynthesis(double pressure, const TInitInfo info)    //Incident PFD, Air temp in C, CO2 in ppm, RH in percent
{
    const double    f = 0.2;             //spectral correction
    const double    O = 205;             // gas units are mbar
    const double    theta = 0.7;		//for potato
    const double    scatt = 0.2;        //leaf reflectance + transmittance
	const int		Kc25 = 404; //*MM constant of rubisco for CO2 of C3 plants (de Pury and Farquar, 1997)
	const int		Ko25 = 248; //*MM constant of rubiscuo for O2 from above reference
    const long      Eac = 59400;  const long       Eao = 36000;   // activation energy values
	//double alpha, Kc, Ko, gamma, Ia, Tk, Jmax, Vcmax, TPU, J, Av, Aj, Ap, Assim, An, Km, Ca, Ci, da, P;
	double alpha, Kc, Ko, gamma, Ia, Tk, Jmax, Vcmax, TPU, J, Av, Aj, Ap, Assim, An, Km, Ca, da, P;
	double nextCi,lastCi, newci, dnewci, c_age, RMx;
	double f_age;
	int iter;
    An = 0;
    //c_age = (log(Parameters.d2)-log(Parameters.d1))/(Parameters.d2-Parameters.d1);     // age at the critical point}
    //RMx = 1/((1-exp(-Parameters.d1*(c_age)))*(exp(-Parameters.d2*c_age)));  // 1/maxima = setting the Y range of the function between 0 and 1}
    //f_age = RMx*(1-exp(-Parameters.d1*c_age))*(exp(-c_age*Parameters.d2)); // function for age effects modified from the original} - DHF - used chronological leaf age?
	f_age = this->leaf_age; //weighted leaf age effect on Pgross for entire canopy
	f_age = 1;
   	gamma = 36.9 + 1.88*(Tleaf-25)+0.036*Square(Tleaf-25);  // CO2 compensation point in the absence of mitochondirial respiration, in ubar}

//* Light response function parameters */
	Ia = PFD*(1-scatt);    //* absorbed irradiance */
	alpha = (1-f)/2; // apparent quantum efficiency, params adjusted to get value 0.3 for average C3 leaf

//* other input parameters and constants */
	Tk = Tleaf + 273.0;
	P  = Press/100;
	Ca = CO2*P; //* conversion to partial pressure */ 
	Kc = Kc25*exp(Eac*(Tleaf-25)/(298*R*(Tleaf+273)));
	Ko = Ko25*exp(Eao*(Tleaf-25)/(298*R*(Tleaf+273)));
	Km = Kc*(1+O/Ko); //* effective M-M constant for Kc in the presence of O2 */
    Jmax = f_age*Parameters.Jm25*exp(((Tk-298)*Parameters.Eaj)/(R*Tk*298))*
		     (1+exp((Parameters.Sj*298-Parameters.Hj)/(R*298)))/
			 (1+exp((Parameters.Sj*Tk-Parameters.Hj)/(R*Tk))); // de Pury 1997
    //Vcmax = f_age*Parameters.Vcm25*exp(Parameters.Eav*(Tleaf-25)/(298*R*(Tleaf+273))); //old response from corn model
    Vcmax = f_age*Parameters.Vcm25*exp(((Tk-298)*Parameters.Eav)/(R*Tk*298))*
		     (1+exp((Parameters.Sv*298-Parameters.Hv)/(R*298)))/
			 (1+exp((Parameters.Sv*Tk-Parameters.Hv)/(R*Tk))); // Used peaked response, DHF
    TPU = f_age*Parameters.TPU25*exp(Parameters.Eap*(Tleaf-25)/(298*R*(Tleaf+273)));
    Ci = Ca;
	nextCi = Ci;
	iter = 1;
	lastCi = 0;
	while ((fabs(Ci-lastCi) > 0.01) && (iter < maxiter))
	{
		/// iteration process to couple photosynthesis model with stomatal model : newton-raphson}
		/// let f(ci) = newci-ci where newci=g(ci), hence df(Ci)/dci=dg(Ci)/dci-1. Solve dg(Ci) for each region}
		lastCi = Ci;
		gs = gsw(pressure, info);
		Ci = __min(__max(0.0,Ca - A_net*(1.6/gs+1.37/gb)*P),2*Ca);//energy balance
        Av = (Vcmax*(Ci-gamma))/(Ci+Km);
        J =  (((alpha*Ia + Jmax) - sqrt(Square(alpha*Ia+Jmax) - 4*alpha*Ia*(Jmax)*theta)) / (2*theta)) ;
        Aj = J*(Ci-gamma)/(4*(Ci+2*gamma));
        Ap = 3*TPU;
       // if (((Square(alpha*Ia+Jmax) - 4*alpha*Ia*Jmax*theta) > 0) && (Ci > gamma))
		//{
           An = Min(Av,Aj,Ap)-Rd();
           // da = Vcmax*(Km+gamma)/((Ci+Km)**2) + 3*J*gamma/(4*(Ci+2*gamma)**2);
		  // if (Min(Av,Aj,Ap) == Ap) da = 0;
           //else 
			 //  if (Min(Av,Aj,Ap)== Av) da = Vcmax*(Km+gamma)/(Square(Ci+Km));
			  // else da = 3*J*gamma/(4*Square(Ci+2*gamma));
		//}
        A_net = An;
	    iter++;
	}
    iter1 = iter;
    //A_net = __max(An,0);
    A_gross = __max(A_net+Rd(),0);
    gs = gsw(pressure, info);
    //Rd = Rd();
}


void CGas_exchange::EnergyBalance(double Jw)
// see Campbell and Norman (1998) pp 224-225
// because Stefan-Boltzman constant is for unit surface area by denifition,
// all terms including sbc are multilplied by 2 (i.e., gr, thermal radiation)
{
    const long lambda = 44000; //J mol-1 at 25C
    const double psc = 6.66e-4;
    const double Cp = 29.3; // thermodynamic psychrometer constant and specific hear of air
	double gha, gv, gr, ghr, psc1, Ea, thermal_air, Ti, Ta;
    Ta = Tair;
    Ti = Tleaf;
    //gha = gb*(0.135/0.147);  // heat conductance, gha = 1.4*.135*sqrt(u/d), u is the wind speed in m/s} Note: this was only true if stomatal ratio = 1
	gha = 1.4 * 0.135 * sqrt(__max(0.1,wind)/(width/100*0.72));
    gv = gs*gb/(gs+gb); 
    gr = (4*epsilon*sbc*pow(273+Ta,3)/Cp)*2; // radiative conductance, 2 account for both sides
    ghr = gha + gr;
    thermal_air = epsilon*sbc*pow(Ta+273,4)*2; // emitted thermal radiation
    psc1 = psc*ghr/gv; // apparent psychrometer constant
    VPD = Es(Ta)*(1-RH); // vapor pressure deficit
    Ea = Es(Ta)*RH; // ambient vapor pressure
//	if (Jw == 0) //i.e., transpiration is zero
	{
        Tleaf = Ta + (psc1/(Slope(Ta) + psc1))*((R_abs-thermal_air)/(ghr*Cp)-VPD/(psc1*Press)); //eqn 14.6b linearized form using first order approximation of Taylor series
	}
//	else
	{
//		Tleaf = Ta + (R_abs - thermal_air - lambda*Jw)/(Cp*ghr); //use direct simulated value of latent heat loss (i.e., Jw is average water loss from the leaf, mol m-2 leaf s-1
	}
    ET = __max(0,1000*gv*((Es(Tleaf)-Ea)/Press)/(1-(Es(Tleaf)+Ea)/(Press)));//1000 is to go from mol to mmol
    // accounting for additional transp. because of mass flow, see von Caemmerer and Farquhar (1981)
}
 


double CGas_exchange::gsw(double pressure, const TInitInfo info)  // stomatal conductance for water vapor in mol m-2 s-1 
{
	double  Pn, aa, bb, cc, Ha, Hs, Cs, Ca, gg, gamma, Ds;
	double temp = set_PSIleafeffect(pressure, info);
    gsModel    myModel=BBW;     // need to put this in main program 
    Ca = CO2;
    Ha = RH;
    Cs = Ca - (1.37*A_net/gb); // surface CO2 in mole fraction
    gamma = 36.9 + 1.88*(Tleaf-25)+0.036*Square(Tleaf-25);  // CO2 compensation point
	if (Cs==gamma) Cs = gamma + 1;
	//if (f_age == 0) gg = 0.001; 
    //  else gg = Parameters.g0*f_age; // prevent division by zero
	gg = Parameters.g0;
	if (A_net <= 0) Pn = 0.00001;
       else Pn = A_net; // solving quadratic formula for surface humidity(hs)
    aa = temp*Parameters.g1*A_net/Cs;
    bb = gg+gb-(Parameters.g1*Pn/Cs);
    cc = (-Ha*gb)-gg;
    Hs = (-bb+sqrt(bb*bb-4*aa*cc))/(2*aa);       // RH at leaf surface
	if (Hs > 1) Hs = 1;
    Ds = (1-Hs)*Es(Tleaf); // VPD at leaf surface - 
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


double CGas_exchange::gbw(void) //boundary layer conductance to vapor
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

double CGas_exchange::Es (double T) //Campbell and Norman (1998), p 41 Saturation vapor pressure in kPa
{
    return (0.611*exp(17.502*T/(240.97+T)));
}

double CGas_exchange::Rd()   //Should be expanded to include other env. and physiological factors
{
	const long Ear =  66400; //exponential rate of arrhenious function for mitochondrial respiration (J mol)
	return (Parameters.Rd25*exp(Ear*(Tleaf-25)/(298*R*(Tleaf+273)))); 
	//return 0;
}

double CGas_exchange::set_PSIleafeffect(double pressure, const TInitInfo info)
{
	//Reduction in stomatal conductance using hourly bulk leaf water potential in MPa
	double sf, phyf;
	sf = 2.3; phyf = -1.2;
	//sf = 6.1; phyf = -0.67; //sensitivity parameter Tuzet et al., (2003), tuned for potato using Vos and Oyarzun (1987) by DHFreference potential Tuzet et al., (2003), " by DHF
	//sf = 4.0; phyf = -0.4; //as above, but using the Liu et al., 2009 parameters for potato for g0 and g1
	if (pressure < -0.05) 
		this->psileaf_stress = __max((1+exp(sf*phyf))/(1+exp(sf*(phyf-pressure))),0);
	else psileaf_stress = 1;
	//psileaf_stress = 1;
	double temp = this->psileaf_stress;
	if (info.Water_stress_off == 1) temp = 1; //remove waterstress if stress option is turned off
	if (info.Water_stress_simulation_type == 1) temp = 1;//don't use if SIMPOTATO approach employed
	return temp;
}


double CGas_exchange::Slope(double T) // slope of the sat vapor pressure curve: first order derivative of Es with respect to T
{
   const double b= 17.502; const double c= 240.97;
   return (Es(T)*(b*c)/Square(c+T)/Press); 
}

void CGas_exchange::SetVal(double PFD, const TInitInfo info, double Tair, double CO2, double RH, double wind,  double Press, double width, double Tubmod, double Feedback, double leafage, double psil, double etsupply, double LAIclass)
{
	/***********************************/
	/* Sets initial parameters and runs routines
	/* for the coupled biochemical C3 gas exchange, stomatal conductance, energy balance routine
	/* orignally coded by Soo Kim
	/***********************************/

	const  double scatt = 0.15;
	this->PFD = PFD;
    double PAR = (PFD/4.55); //w par m-2
    double NIR = PAR; // If total solar radiation unavailable, assume NIR the same energy as PAR waveband
    this->R_abs = (1-scatt)*PAR + 0.15*NIR + 2*(epsilon*sbc*pow(Tair+273,4)); // times 2 for projected area basis
	// shortwave radiation (PAR (=0.85) + NIR (=0.15) solar radiation absorptivity of leaves: =~ 0.5
    this->CO2 = CO2;
	//this->CO2 = 370;
    this->RH = __min(100.0, __max(RH, 25.0))/100;
    this->Tair = Tair;
    this->width = width;
	this->wind = wind;
    this->Press = Press;
	this->psileaf = psil;
	this->leaf_age = leafage;
    getParms(Tubmod, Feedback, LAIclass);
    //GasEx();
	GasEx_psil(psileaf, etsupply, info);
}


double QuadSolnUpper (double a, double b, double c )
{
    if (a==0) return 0;
	else if ((b*b - 4*a*c) < 0) return -b/a;   //imaginary roots
    else  return (-b+sqrt(b*b-4*a*c))/(2*a);
}

double QuadSolnLower (double a, double b, double c )
{
    if (a==0) return 0;
	else if ((b*b - 4*a*c) < 0) return -b/a;   //imaginary roots
    else  return (-b-sqrt(b*b-4*a*c))/(2*a);
}
