#include "stdafx.h"
#include "leaf.h"
#include "weather.h"
#include <math.h>

#define R 8.314  // idealgasconstant
#define maxiter 100 

inline double Square(double a) { return a * a; }
inline double Min4(double a, double b, double c, double d) {return __min(a,__min(b,__min(c,d)));}


CLeaf::CLeaf()
: COrgan()
{
}

CLeaf::CLeaf(const TInitInfo info, int n, int totalNo): COrgan()
{
	/************************/
	/* Adds to CORGAN initiation
	/* Initializes values of leaf
	/****************************/
	double dt = info.timeStep/(24*60); //converting minute to day decimal, 1= a day
	length=width= SLA = potentialArea = GDD = deadmass = 0.0;
	if (info.bigleaf == 0){
		area = 0.05; //from spar2004 measurements (leaf appearance measured as 0.05cm2 leaf at apearrance)
		if (n <= 3) //leaves at emergence more mass than leaves later in season
		{
			set_drymass(0.0025);
			set_drymass_at_emergence(0.0025);
			set_CH2O(0.0025);
		}
		else //new leaves after emergence have initial SLA of 100 cm2 g-1
		{	set_drymass(0.0005);
			set_CH2O(0.0005);
		}
	}
	else{
		area = 0.2;
		set_drymass(0.01); //from spar2004 measurements for 4 nodes
		set_drymass_at_emergence(0.01);//"
		set_CH2O(0.01);
	}
	greenArea = area;
	area0 = area;
	deadCHOpool = 0.;
	deadNpool = 0.;
	senstressArea=0.;
	senageArea=0.;
	psi_leafexpansion_stress = 1.;
	length = 1.; //initial leaf appeared defined by 1cm x 0.5 cm length by width
	width = 0.05;
	SLAleaf = SLAnom = 1.0/info.G4;//was 300, now based on G4 input parameter
	ID = n;
	//dblSLAmin = 123; // from daylit data at high T - may need to refine this value with more data
	//dblSLAmin = 120; // What are the minimum and maximum SLA values for individual leaves?
	//dblSLAmax = 300; // ?
	SLAmin = 175.; //was 175
	//SLAmin = 125;
	SLAmax = 225.;
	//dblMaxgrowth = 10;
	dblMaxgrowth = info.G2;
	if (dblMaxgrowth > 50) dblMaxgrowth = 10.; /*in case user forgets to adjust G2 parameter for single leaf instead of big leaf*/

	T_base = 0.0;  T_opt = 27.98; T_ceil = 40.46; Rmax_leaf = 0.979*dt;  //from Fleisher et al (2005) SPAR data
	old = initiated = prolific = aging = terminated = false;
	young = growing = true;
	maintCoeff = 0.0408; //gCH2O g-1DM day-1 at 25C for mature potato leaves, Ng and Loomis
	//set_optimumNitrogenRatio(0.08);//g N required per g new leaf CHO to get maximum growth (non-stressed), from SIMPOTATO 
	//set_actualNitrogenRatio(0.07); //
	//set_minimumNitrogenRatio(0.07); // from SIMPOTATO
	set_optimumNitrogenRatio(0.045);//from SUBSTOR
	set_actualNitrogenRatio(0.045);//from SUBSTOR
	set_minimumNitrogenRatio(0.025);//from SUBSTOR
	set_currentNitrogenAmount(this->get_actualNitrogenRatio() * this->get_drymass());
}


void CLeaf::expand(int icur, const TInitInfo info, int potential, const TNitrogenStatus nitrogen, double CH2O, CDevelopment * dv, double predawn_LWP, const TWeather &weather)
/***********************************/
/* Bigleaf =  0: 
/*	Uses leaf expansion routine (Fleisher and Timlin, 2007) - modified Ng and Loomis
/*	Estimates potential area expansion and C demand for current leaf each time-increment
/*	potential leaf area growth when potential = 0, actual when potential = 1 -- allows reuse of same routine
/*	For potential growth
/*   -Get 24 hour potential expansion based on prior 24 h temperature
/*	 -Convert to 1 hour time-step and estimate C demand using SL A
/*	For actual growth
/*   -Given C amount, convert to area using SLA
/*	Effect of leaf water potential on expansion from MAIZSIM is also included here
/*	Need to add N stress effect as well?  Does N direclty limit leaf expansion though or it as a result of reduced photosynthesis?
/* Bigleaf = 1:Note that if big leaf model is used:
	-potential growth for entire canopy is stored in nodalunit(0) and is strictly based on thermal time responses, and genetic parameter for canopy expansion
	-similar idea as above wrt to converting 24 h to 1 h growth basis
	-N stress factors in here in terms of reducing potential expansion - uses 24h basis as in SIMPOTATO
	-NEED to implement water stress effect
/***********************************/

{
	double dblGrowth, x, fage, fT,id;
	double dblnostressGrowth;
	//Original routine was based on 24 hour previous Temperature and leaf area
	//Therefore, we compute new 24 hour growth when each leaf grows for 24 hour on its own chronological time schedule
	//Also means must compute dblGrowth (below) only once at begining of each day//
	if (info.bigleaf==0)
	{
		if (potential == 0)
		{
			x = this->get_Cage();
			fage = get_agefraction();
			fT = get_Tfactor();
			id = this->get_ID();
			if (x == 0)  //this condition will never be met? leaf will have aged at least 1 time-step prior to update
			{
				set_area0(area);
				potGrowth24 = area*dblMaxgrowth*get_agefraction()*get_Tfactor() / (24.*60./info.timeStep);
				
			}
			if ((((x - floor(x)) <= (1.25*info.timeStep/(24.*60.)))&&((x-floor(x)) > 0.0001))) 
			//if ((x-floor(x))<=info.timeStep/(24*60))
			{
				set_area0(area);
				//potGrowth24 = area0*dblMaxgrowth*get_agefraction()*get_Tfactor()*get_LWPeffect(predawn_LWP) / (24*60/info.timeStep);
				potGrowth24 = area0*dblMaxgrowth*get_agefraction()*get_Tfactor() / (24.*60./info.timeStep); //psileaf effect to occur on hourly basis
			}
			if (potGrowth24 < 0) potGrowth24 = area0*dblMaxgrowth*get_agefraction()*get_Tfactor()/ (24.*60./info.timeStep);
			dblGrowth = potGrowth24 * get_LWPeffect(predawn_LWP, weather, info) ; //get growth per time increment and adjust for current water status
			if (dblGrowth < 0) dblGrowth = 0.;
			potentialArea = dblGrowth+area;
			potentialAreainc = dblGrowth;
			potentialDW = get_drymass() + potentialAreainc / SLAleaf;
			potentialDWinc = potentialAreainc / SLAleaf;
		}
		if (potential > 0) //compute actual growth: right now, assumes fixed SLA
		{
				DWinc = CH2O;
				//DWinc = potentialDWinc; //only use this line if non C-stressed
				area += DWinc * SLAleaf;
				greenArea = area;
				length = pow((area / 0.872 * 3.),0.5); //from 2005 leaf expansion paper - assume L:W of 3:1
				width = length / 3.; 
		}
	//			dblPgrowth = dblGrowth*dblArea; // potential increase in area, delta cm2
	//			dblSLA = dblPgrowth / CHO;  // resulting SLA for that increment based on CHO supply for leaf
				//if (dblSLAmax >= dblSLA >= dblSLAmin){ //SLA for leaf increment is above or equal to minimum value, so don't decrease expansion
	//			if ((dblSLAmax >= dblSLA) && (dblSLA >= dblSLAmin)){		
	//				dblPweight = dblPgrowth / dblSLA; // potential increase in weight, delta g CHO
	//				dblArea = dblArea*dblGrowth + dblArea; // total current leaf area after increase, cm2
	//				dblWeight = dblWeight + dblPweight; // total current leaf weight after increase, g
	//			}
				//if (CHO >= dblGrowth*dblArea/dblSLAnom) // is there enough CHO to fully increase leaf area by computed increment?
				//{
				//	dblPgrowth = dblGrowth*dblArea; // potential increase in area, delta cm2
				//	dblPweight = dblPgrowth / dblSLAnom; // potential increase in weight, delta g CHO
				//	dblArea=dblArea*dblGrowth + dblArea; // total current leaf area after increase, cm2
				//	dblWeight = dblArea / dblSLAnom; // total current leaf weight after increase, g
				//}
	//			else // if not, reduce current leaf growth to the max SLA based on CHO value
	//			{
	//				dblSLA = dblSLAmax;
	//				dblGrowth = CHO*dblSLA/dblArea;  // max possible growth is limited to available CHO, 
	//				dblPgrowth = dblGrowth*dblArea;
	//				dblPweight = dblPgrowth / dblSLA;
	//				dblArea=dblArea*dblGrowth + dblArea;
	//				dblWeight = dblWeight + dblPweight;
	//			}
	//		}
	}
	if (info.bigleaf == 1)
	{
		if ((icur == 1) || (fmod(icur,((24.0*60.0)/info.timeStep))==0))//update potential leaf expansion at 24h basis at beginning of each day
		{
			if ((dv->get_CDTT22()<=25)||(dv->get_istage()==1 && this->get_greenArea()*info.plantDensity/10000.<0.02)){
				//potGrowth24=info.G2/25*dv->get_CDTT22()*dv->get_DTT22() / (24*60/info.timeStep);//DHF - may be too restrictive on early leaf expansion?
				//potGrowth24=exp(0.5*dv->get_DTT22())*this->get_greenArea()-this->get_greenArea() / (24*60/info.timeStep); //SUBSTOR Approach
				potGrowth24=info.G2*dv->get_DTT22()/(24.0*60.0/info.timeStep); //remove early growth restriction as early lag phase in leaf area growth doesn't appear realistic - at least for chamber data
			}
			else{
				potGrowth24=info.G2*dv->get_DTT22()/(24.0*60.0/info.timeStep);
			}
			// To Do - add status to reduce potential growth - '1' below should be water stress factor
			if (info.Nitrogen_stress_off==0) //Right now, only influence on leaf expansion on hourly basis
			{
				//potGrowth24*=Min(1, weather.swdf2, nitrogen.Nitrogendeficiencytwo);
				//potGrowth24*=Min(1, get_LWPeffect(predawn_LWP, weather, info), nitrogen.Nitrogendeficiencytwo);
			}
			else
			{
				//potGrowth24*=Min(1,1,weather.swdf2);
				//potGrowth24*=Min(1,1,get_LWPeffect(predawn_LWP,weather,info));
			}

			if (potGrowth24 < 0) potGrowth24 = 0.;
		}
		if (potGrowth24 < 0) potGrowth24 = 0.;
		if (info.Water_stress_simulation_type==0) dblGrowth = potGrowth24*get_LWPeffect(predawn_LWP, weather, info) ; //get growth per time increment and adjust for current water status; //get overall canopy area growth per time increment
		else
			dblGrowth = potGrowth24*weather.swdf2;
		
		//dblGrowth = potGrowth24;	// currently too restrictive on growth - re-allocation of CHO away from shoot to root accomplishes this reduction indirectly?
									// AD 2012-07-12 dissabled it because too many leaves grow under water stress in SB2005 experiment
		dblnostressGrowth = potGrowth24;  //theoretical growth rate without water stress
		if (potential == 0) // potential growth
		{
			potentialArea = dblGrowth + greenArea;
			potentialAreainc = dblGrowth;
			potentialDW = get_drymass() + potentialAreainc / SLAleaf;
			potentialDWinc = potentialAreainc / SLAleaf; 
			nostresspotentialArea = dblnostressGrowth + greenArea;
			nostresspotentialAreainc = dblnostressGrowth;
			nostresspotentialDW = get_drymass() + nostresspotentialAreainc / SLAleaf;
			nostresspotentialDWinc = nostresspotentialAreainc / SLAleaf;
		}
		if (potential > 0) //compute actual growth: right now, assumes fixed SLA
		{
				DWinc = CH2O;
				area += DWinc * SLAleaf;
				greenArea += DWinc * SLAleaf;
				length = pow((area / 0.872 * 3.),0.5); //from 2005 leaf expansion paper - assume L:W of 3:1
				width = length / 3.; 
		}

	}
	
}


double CLeaf::get_maintResp(const TInitInfo info)
{
	double temp = 0.;
	if (info.bigleaf == 0)
	{
		double age = 20./15.*get_PAge(); //20/15 is scaling factor for Leaf physiological age between Ng and Loomis scale and scale used in SPUDSIM (15 days to Gdur here, 20 days in Ng and Loomis)
		// Routine comes from leaf age effect table from Ng and Loomis (1984) for Pg and Rm//
		if(age == 0) temp = 0.417;
		if(0 < age  && age <=10)
		{
			temp = age*(1-0.417)/10. + 0.417;
		}
		else if (10 < age && age <= 20)
		{
			temp = 1.;
		}
		else if (age > 20)
		{
			temp = -1.*(age-20.)/(83.-20.)+1.;
		}
		return temp*maintCoeff*(get_drymass()-get_drymass_at_emergence());
		//return maintCoeff*get_drymass(); // override Ng and Loomis for now
	}
	else
	{
		if (get_drymass() > 0) temp = maintCoeff*(get_drymass()-get_drymass_at_emergence());
		else
			temp = 0.;
		return temp;
	}
}
 

double CLeaf::get_LWPeffect(double predawn_psil, const TWeather &weather, const TInitInfo info) //From MAIZSIM 12-17-2010
//create a function which simulates the reducing in leaf expansion rate
//when predawn leaf water potential decreases. Parameterization of rf_psil
//and rf_sensitivity are done with the data from Boyer (1970) and Tanguilig et al (1987) YY
// predawn leaf potential needs to be in Bars, not MPa, so note conversion factor in equation DHF
// Assume same relative reduction in leaf expansion weather small or big-leaf

//Modified for SPUDSIM using Gandar and Tanner 1976 data in MPa based on current hourly psil and not predawn
{
	//double rf_psil=-1.87; //corn
	//double rf_sensitivity=1.92; //corn
	double rf_psil = -0.31; // potato threshold at which leaf expansion starts to decline
	double rf_sensitivity = 11.6; // based on best model fit
	
	//double rf_psil=-0.4;
	//double rf_sensitivity = 4.;
	
	
	double effect;
	//effect=(1+exp(rf_psil*rf_sensitivity))/(1+exp(rf_sensitivity*(rf_psil-predawn_psil*10)));
	//if (weather.psil_ < -0.05) {
	//if (weather.psil_ < rf_psil) {
	//	effect=max((1+exp(rf_psil*rf_sensitivity))/(1+exp(rf_sensitivity*(rf_psil-weather.psil_))),0);
	//}

	//if (predawn_psil/10 < rf_psil) {
	if (predawn_psil/10. < -0.05) {
		effect=__max((1+exp(rf_psil*rf_sensitivity))/(1+exp(rf_sensitivity*(rf_psil-predawn_psil/10.))),0.);
		psi_leafexpansion_stress = __max((1 + exp(rf_psil * rf_sensitivity)) / (1 + exp(rf_sensitivity * (rf_psil - predawn_psil / 10.))), 0.);
	}
	else {
		effect = 1.;
		psi_leafexpansion_stress = effect;
	}
	if (info.Water_stress_off == 1) {
		effect = 1.;
		psi_leafexpansion_stress = effect;
	}
	return effect;
}

void CLeaf::set_SLA(const TInitInfo info, double Sradave, double LAI)
{
	/*******************************************/
	/* Determines nominal specific leaf area
	/* Big leaf 1 model is big canopy simulation from SIMPOTATO and assumes nominal SLA is genetic
	/* Big leaf 0 uses a modified Ng and Loomis approach for individual leaf*/
	/*		-includes physiological aging and canopy shading effect of new increment of leaf area
	/* This probably should be modified with N? since potato tends to keep Pmax with smaller leaves in response to N stress
	/*	-need to determine if SLA is actually influenced by N stress, however...
	/********************************************/


	if (info.bigleaf == 0)
	{
		double effShade = 0, effAge = 0;
		double age = 0;
		double ratio = 0;
		age = this->get_PAge()*20./15.; //20/15 adjusts for aging scale diffs between SPUDSIM and Ng and Loomis
		//age = this->get_Cage();
		ratio = (Sradave * 1000000.)/41870. / LAI * 0.8; //convert MJ m-2 d-1 to LY m-2 d-1
		//Age affect as in Ng and Loomis
		if (age <= 8) effAge = 1.;
		else if (age <= 15) effAge = (0.8-1.)/(15.-8.)*(age-8)+1.;
		else if (age <= 20) effAge = 0.8;
		else if (age > 20) effAge = (-0.8)/(200.-20.)*(age-20.)+0.8;
		//Shading effect as in Ng and Loomis
		if (ratio <= 125) effShade = (0.933-1.)/125.*(ratio)+1;
		else if (ratio <= 325) effShade = (0.68-0.933)/(325-125)*(ratio-125)+0.933;
		else if (ratio <= 500) effShade = (0.653-0.68)/(500-325)*(ratio-325)+0.68;
		else if (ratio <= 10000000) effShade = (0.64-0.653)/(10000000.-500.)*(ratio-500.)+0.653;
		else effShade = 0.64;
	//	effShade = 1;
		set_SLAleaf(SLAnom*effShade*effAge);
		//set_SLAleaf(125);
		set_SLAleaf(SLAnom*effAge);
		//if (this->get_PAge() < 5) set_SLAleaf(300);
		//set_SLAleaf(SLAnom*effAge);
	}
	if (info.bigleaf == 1) this->set_SLAleaf(1./info.G4);
	//set_SLAleaf(SLAnom);//comment out this line
}

void CLeaf::set_potentialArea(const TInitInfo info)//* may want to use this routine to adjust G2 with leaf position
{
	//	initiated = true;
}

void CLeaf::update(int icur, CDevelopment * dv, const TInitInfo info, const TWeather &weather, const TNitrogenStatus nitrogen, double Tlagleaf, int potential, double CH2O, double predawnLWP, double swdf1, double C_pool_room)
{
	/*******************************/
	/* Update leaf aging, potential expansion, senescence, age effect on photosynthesis
	/*******************************/

	COrgan::update(icur, info, weather, Tlagleaf); //get aging and thermal-related growth responses
	if (get_PAge() <= 0.3*get_growthDuration()) young = true; //assume leaves 30% of final age are still dependent on import from other C sources
	else
		young = false;
	if (get_PAge() > 0.7*get_growthDuration()) old = true;
	else
		old = false;
	if (get_PAge() >= get_growthDuration()) growing = false;
	expand(icur, info, potential, nitrogen, CH2O, dv, predawnLWP, weather);
	senescence(icur, dv,info, weather, nitrogen, Tlagleaf, swdf1, C_pool_room);
	set_ageEffect(info);
}

void CLeaf::set_ageEffect(const TInitInfo info)
{
	// Following routine is the physiological age effect of each leaf on Pgross
	// Values come from Ng and Loomis (1984)
	// Curretnly modified response so that leaves < 15 Pdays still have max Pgross parameters
	// If big leaf model is used, there is no simulation of ageeffect
	if (info.bigleaf == 0)
	{
		double x,y;
		/*
		y = 1;
		x = get_PAge()*20/15;//20/15 adjusts for difference in aging scale between Spudsim and Ng and Loomis
		//if (x <= 0) y = 0.417;
		//else if (x <= 10) y = (1-0.417)/(10)*(x-0)+0.417;
		if (x <= 20) y = 1; //my modification of Ng and Loomis (2 lines above)
		//else if (x <= 20) y = 1;
		else if (x > 20) y = (0-1)/(83-20)*(x-20)+1;
		if (y < 0) y = 0;
		ageEffect = y;
		*/
		//ageEffect = 0.2242+(0.9951-0.2242)/(1+pow(this->get_PAge()/7.5939,5.7234));// 4 parameter logistic curve developed for potato, 2002 data
		x = get_PAge();
		if (x < 7.5) y = 0.0945*x+0.3152;
		else if (x <= 15) y = 1.;
		else if (x <= 35) y = -0.0433*x+1.665;
		else y = 0.15;
		ageEffect = y;
		//ageEffect = 1;//right now, remove ageEffect - this should be included in leaf gas exchange routine
	}
	else
	{
		ageEffect = 1.;
	}
}

double CLeaf::get_deadCHOpool()
{
	double x;
	x = this->deadCHOpool;
	this->deadCHOpool = 0.;
	return x;
}

void CLeaf::senescence2(const TInitInfo info)
{
	/*************************/
	/* Forced senescence of current leaf due to excessive canopy shading*/
	/* Called up during plant->update_nodes class*/
	/* For big leaf = 0, estimate of senescence done in plant.cpp, update_nodes method
	/* For big leaf = 1, esitmate is done in this routine */
	/*************************/
	double shadeloss = 0;

	if (info.bigleaf == 0)
	{
		this->terminated = true;
		deadmass = 0.5*get_drymass();
		deadCHOpool = 0.5*get_drymass();
		deadNpool = get_drymass() * (get_actualNitrogenRatio() - get_minimumNitrogenRatio());
		deadNmass = get_drymass() * get_minimumNitrogenRatio();
		this->set_drymass(0.);
		this->set_CH2O(0.);
		this->set_area(0.);
		this->set_currentNitrogenAmount(0.);
		greenArea = 0.;
	}
	else //already accounted for in senescence routine for bigleaf = 1!
	{
		/*
		shadeloss = this->get_greenArea() * (1-(1-0.008*((this->get_greenArea()*info.plantDensity/10000)-4)))/(24*60/info.timeStep);
		if (shadeloss > 0) {
			deadmass += 0.5*shadeloss*info.G4;
			deadCHOpool += 0.5*shadeloss*info.G4;// programming note: needs to be += here since first senescence routine is called earlier in time-step
			//deadNpool = deadCHOpool * (get_actualNitrogenRatio() - get_minimumNitrogenRatio());
			deadNpool = deadCHOpool * get_actualNitrogenRatio() * 0.75;
			//deadNmass = deadCHOpool * get_minimumNitrogenRatio();// programming note: as above
			deadNmass = deadCHOpool * get_actualNitrogenRatio() * 0.25;
			greenArea -= shadeloss;
			this->set_drymass(this->get_drymass()-shadeloss*info.G4);
			this->set_CH2O(this->get_CH2O()-shadeloss*info.G4);
			this->set_currentNitrogenAmount(this->get_currentNitrogenAmount()-deadNmass);
		}
		*/
	}
}

void CLeaf::senescence(int icur, CDevelopment * dv, const TInitInfo info, const TWeather& weather, const TNitrogenStatus nitrogen, double Tlagleaf, double swdf1, double C_pool_room)
{
	/*******************************/
	/* Simulate individual leaf senescence (big leaf = 0) or whole canopy senesence (bigleaf = 1)
	/* If big leaf = 0 (individual organ model:):
		/* If leaf reaches 35 days of physiological age it is senesced*/
		/* Remove 50% of leaf CHO and move to deadCHOpool
		/* Remove minimum N content required for structure and return rest of N to deadNpool
		/* Set greenarea and dry mass and N content of leaf to 0 since there is no longer green leaf material
	/* If big leaf model is used:
		-Senescence is based on N and water stress, plant age, current greenleaf area, and thermal responses as in SIMGUI
		-Derived from portion of GROSUB routine - comments in that section were from SIMPOTATO
		-Note, SIMGUI assumes 25% N of the 50% part of dead leaf lost to senescence is lost too
	/******************************/

	double slan;//thermal time based senescence from SIMGUI (cm2 leaf area)
	double plas;//total senescence from SIMGUI for time increment (cm2 leaf area)
	double slfw, slfn, slfc, slft; //water, n, c, and T stress effects on senescence from SIMGUI (0 to 1 unitless)
	double Clost = 0.5; //how much of leaf mass is lost when senesced

	if (info.bigleaf == 0)
	{
		if (get_PAge() >= 35 && !this->isTerminated())
		{
			this->terminated = true;  //assumed 35 physiological days of aging for leaf senescence
			deadmass = 0.5*get_drymass();
			deadCHOpool = 0.5*get_drymass();
			deadNpool = get_drymass()*(get_actualNitrogenRatio() - get_minimumNitrogenRatio());
			deadNmass = get_drymass()*get_minimumNitrogenRatio();
			this->set_drymass(0.);
			this->set_CH2O(0.);
			this->set_area(0.);
			this->set_currentNitrogenAmount(0.);
			greenArea = 0.;
		}
	}
	else //big leaf canopy - values originally calculated based on 24 hour increment.
	{
		// from SIMGUI, values originally were based on 24 h time increment
		// to adjust here, use greenleaf area at time of next day, then divide by 24//
		slan = 0.; plas = 0.; slfw = 1.; slfn = 1.; slfc = 1.; slft = 1.;
		if (fmod(icur,((24.0*60.0)/info.timeStep))==0)
		{
			//natural aging (thermal time) stress as in SIMGUI
			
			if (dv->get_istage() < 3)
				senageArea = dv->get_CDTT22()*this->get_greenArea()/10000.; // SIMGUI
			else if (dv->get_istage() == 3)
			{
				senageArea = this->get_greenArea()/(500.0*3.0/dv->get_xstage());
				if (dv->get_xstage() > 4)
				{
					senageArea = this->get_greenArea() / (500. / dv->get_xstage());/*original*/
					senageArea = this->get_greenArea() / ((500. - (60. * (dv->get_xstage() - 4.))) / dv->get_xstage());
					/*idea here is to scale increasingly with age over original routine, wasn't able to get rapid enough senesnce*/
				}
				if (dv->get_xstage()>10)
				{ 
					senageArea = this->get_greenArea()/(150./dv->get_xstage()); /*DHF added to hasten senescence towards end of season*/
				}
			}		
			//natural aging stress as in SUBSTOR and modified by Fleisher et al.,2003, exponent value was -1.6, but too rapid decline in senescence rate
			/*
			if (dv->get_istage() < 2) 
				senageArea = dv->get_CDTT22()*this->get_greenArea()/10000;
			else if (dv->get_istage() < 4)
			{
				//int temp = dv->get_xstage();
				senageArea =  dv->get_CDTT22()*this->get_greenArea()/10000*exp(-2+0.8*dv->get_xstage())*pow(info.G1,0.5);	

			}
			*/
		// calculate water, nitrogen, carbohydrate and temperature stress effects on leaf senescence
		
			slfn = 0.95 + 0.05*nitrogen.Nitrogendeficiencytwo;
			slfw = 0.95 + 0.05*swdf1;
			
			//during stage 3, allow water and N stress to senesce leaf area more rapidly
			if (dv->get_istage() > 2) slfn = 0.9 + 0.1*nitrogen.Nitrogendeficiencytwo;
			if (dv->get_istage() > 2) slfw = 0.9 + 0.1*swdf1;

			//slfn = 1;
			if (this->get_greenArea()*info.plantDensity/10000. > 4.) slfc = 1.-0.008*(this->get_greenArea()*info.plantDensity/10000.-4.);
			//total leaf kill at 0C is too drastic.  This will cause gradual leaf kill down to -7C.  Based on conversation with Gene Eastin.
			//however, John Hess metnions a year when light freeze in August stopped tuber growth without killing leaves	
			if (Tlagleaf <= 6)
			{
				//slft = 1 - (6 - Tlagleaf)/6; //from SUBSTOR, except Lagleaf used instead of Tmin
				if (Tlagleaf <= 0) 
				{
					slft = 1. - 0.02*pow(Tlagleaf,2.);
				}
			}
//			slfw = 1; //test
			senstressArea = (this->get_greenArea()*(1.-Min4(slfw,slfc,slft,slfn)));
			if (senstressArea > senageArea) senageArea = senstressArea;
			
		}
		slan = senageArea/(24.0*60.0/info.timeStep);
		

		//determine if there is room available for C to be mobilized from dropped leaves
		double tempC = Clost*slan*info.G4;
		if (tempC >  (0.9*C_pool_room)) {
			Clost = 1.; //that is, if there's no room left, no CHO is mobilized
		}
		deadmass += Clost*slan*info.G4; //Clost% whole leaf is dropped from plant
		deadCHOpool = slan*info.G4*(1.-Clost); //but (100-cLost)% is mobilized back into C_deadpool for immediate use - this gets accounted for later
		//deadNpool = deadCHOpool*(get_actualNitrogenRatio()-get_minimumNitrogenRatio());
		deadNpool = deadCHOpool*get_actualNitrogenRatio()*0.75;//gets returned to plant
		//deadNmass = deadCHOpool * get_minimumNitrogenRatio();
		deadNmass = deadCHOpool*get_actualNitrogenRatio()*0.25;//lost to ground litter
		greenArea -=(slan);
		this->set_drymass(this->get_drymass()-slan*info.G4);
		this->set_CH2O(this->get_CH2O()-slan*info.G4);
		//this->set_drymass(this->get_drymass()-deadCHOpool*2);
		//this->set_CH2O(this->get_CH2O()-deadCHOpool*2);

		this->set_currentNitrogenAmount(this->get_currentNitrogenAmount()-deadNmass - deadNpool); //deadNpool will be added back into plant in update_nodal plant.cpp routine
	}

}



void CLeaf::calcLongevity()
// see Lizaso et al. (2003)
{
	const double L0 = 150.0;
	const double Lx = 850.0;
//	double LN_l = 3.59 + 0.498*totalLeaves; // nodal position of most longevous leaf
//	double Wi = (1.0/3.0)*totalLeaves; // Width function of the bell shape
//	double LL =L0 + Lx * exp(-pow(rank-LN_l,2)/(2*pow(Wi,2)));
//	set_longevity(LL);
}

double CLeaf::GTI(double T_avg)
// general thermal index
// improved calculation of GDD, Steward et al. (1998)
{
	double b1 = 0.0432;
	double T_opt = 32.2;
	return b1*T_avg*T_avg*(1-0.6667*T_avg/T_opt);
}
CLeaf::~CLeaf() {}
