//
#include "stdafx.h"
#include "plant.h"
#include "gas_exchange_new.h"
#include "radtrans.h"
#include "timer.h"
#include "math.h"
#include <vector>
#include <iostream>
#include <sstream>
#define PRIMORDIA 4   //Leaf numbers visible (appeared vs emerged) at emergence (from msmts 2003-05 BARC, MD
#define CO2_MW 44.0098
#define C_MW 12.011
#define CH2O_MW 30.03

using namespace std;

CPlant::CPlant()
{
}

CPlant::CPlant(const TInitInfo info )
{
	nodalUnit = NULL;
	tubers = NULL;
	roots = NULL;
	develop = NULL;
	initInfo = info;
	C_seed = initInfo.Seedreserve;
	C_reserve = initInfo.Seedreserve;
	totalMass = C_reserve; // weight of tuber seed piece as input from initialization file
	EffectiveCanopySLA = 0.;
	C_content = 0.4; // 40% C
	assimilate_old = 0.;
	C_pool = 0.; 
	C_pool_root = 0.;
	C_deadpool = 0.;
	C_deadpool_season = 0.;
	C_pool_used = C_seed_used = 0.;
	C_pool_room = 0.;
	C_demand = C_supply = 0.0;
	C_maxPoolSize = C_reserve;
	PhotosyntheticFeedback = 1.;
	instantResp = 0.;
	maintRespiration = 0.0;
	growthRespiration = 0.0;
	Total_soil_water_uptake = 0.0;
	Hourly_soil_water_uptake = 0.0;
	
	N_seed = initInfo.Seedreserve*0.05;//assumes 5% N content in seed tuber (DW), with g N plant-1
	//this->set_AvailableNitrogen(N_seed / 2 );//assume 1/2 of N_seed is immediately available for growth, in g N plant-1
	N_reserve = N_seed;// N content in seedpiece which may be available for growth, g N plant-1
	//N_reserve = 10000;
	N_seedused = 0.;
	TotalNitrogen = totalLeafNitrogen = totalStemNitrogen = totalRootNitrogen = totalTuberNitrogen = totalDeadNitrogen = 0.;

	dblheatstressVegetative = dblheatstressReproductive = 1.;

	maxLAI = 0.;
	reference = 0.;
	Rgleaf = 0.6;//from ng and Loomis (1984), g CH2O g-1 DM
	Rgorgan =0.7;
	sowingDay = 1.0;
	age = 0.0;
	Bnewnode = false;

	roots = new CRoots(initInfo);
	tubers = new CTubers();
	develop = new CDevelopment(initInfo);
	inodeNumber = igreennodeNumber = 0;
	finalNodeNumber = 0;
	//nodalUnit = new CNodalUnit[100]; // create enough leaf nodes for now, to be replaced by dynamic collection
	//for (int i=0; i <= PRIMORDIA; i++)
	//{
	//	nodalUnit[i].initialize(i,0,0,i,develop);
    //   	nodeNumber = i;
	//}
	//finalNodeNumber = info.genericLeafNo;
	photosynthesis_net =photosynthesis_gross = transpiration = assimilate =  0.0;
	leafArea =greenLeafArea = deadLeafArea = potentialLeafArea = 0.0;
	stemMass = leafMass = tuberMass = rootMass = deadMass = 0.0;
	basalstemNo = apicalstemNo = mainleafNo = basalleafNo = apicalleafNo=0;
	rootPart = shootPart = rootPart_old = shootPart_old = tubPart = tubPart_old = 0.0;

	temperature = 20.; //no meaning for now
	Tcanopyave = 0.;
	counter = 0;
	Tubfraction = Tubmod = 0.;
	leafageEffect = 1.;
	VPD = conductance =  0.; 
	carbonRatio = 1.; //at beginning, assume that carbon supply matches with carbon demand, YY
	swdf1 = 1.;
}

CPlant::~CPlant() 
{
	if (nodalUnit != NULL) delete [] nodalUnit;
	if (roots != NULL) delete roots;
	if (develop != NULL) delete develop;
	if (tubers != NULL) delete tubers;
}

void CPlant::calcGasExchange(const TWeather & weather, const TInitInfo info)
{
	/****************************************/
	//All COMMENTED OUT STATEMENTS ARE FROM SOO"S ORIGINAL CORN MODEL WITH SUNLIT AND SHADED LEAF CLASSES
	/* Simulates leaf level gas exchange for sunlit and shaded leaf areas
	/* Radtrans routine is part of lightenv.dll that estimates sunlit and shaded lai and ppf
	/*   using Campbell and Norman approaches
	/* This information is then used with the Gas_Exchange class to estimate Pg, Pn, Tr and leaf surface for the shaded
	/*  and sunlit leaves
	/* This gets modified by plant tuber sink strength effect (as in Ng and Loomis), age factors, etc.
	/*   Rather than simulating an increase in RUE after tuber initation as in SIMPOTATO, I adopted the Ng and Loomis approach 
	/*   to adjust Pg based on tuber sink strength.
	/* Then Pnet is added to Assimilate pool once adjusted to per unit ground area per hour basis */
	/* NEEDS:
	/*	-adjustment for N stress
	/*  -still have issues matching real-world responses.  Should I just use the leaf model to replace big-leaf area in
	      SSIMPOTATO, but keep all other routines?
	/*******************************************/
	const double tau = 0.50; // atmospheric transmittance, to be implemented as a variable
	const double LAF = 2.2; //1.7 to 2.47 for potato as in Campbell and Norman (1998) - DHF

	const double leafwidth = 5.; //estimate of projected width of leaf to wind direction, value is converted to m in gas exchange method, so use cm width here. 
	const double atmPressure= 101.3; //kPa, to be predicted using altitude
	//double activeLeafRatio = greenLeafArea/leafArea;
	double LAI = get_greenLAI();
	//double activeLeafRatio = 1;
	//double LAI = leafArea*initInfo.plantDensity / (100.0*100.0);
	double poolmod = 1.; //value to reduce photosynthetic parameters based on if demand < photosynthesisx0.8
	int redo = 1; //flag to tell model to loop again in gas exchange if assimilate is much greater than demand
	CGas_exchange_new * sunlit = new CGas_exchange_new();
	CGas_exchange_new * shaded = new CGas_exchange_new();

	CSolar* sun = new CSolar();
	CRadTrans * light = new CRadTrans();
	Timer timer;
	int mm, dd, yy;
	timer.caldat(weather.jday, mm, dd, yy);
	int jday = timer.julday(1, 1, yy);
	/*note that solrad should be PAR + NIR in mw m-2*/
	//radTrans2(weather.jday, weather.time, initInfo.latitude, initInfo.longitude, weather.solRad, weather.PFD, LAI, LAF);
	sun -> SetVal(weather.jday - jday + 1, weather.time, initInfo.latitude, initInfo.longitude, initInfo.altitude, weather.solRad);
	light -> SetVal(*sun, LAI, LAF);
	sunlit_PFD = light ->Qsl();
	shaded_PFD = light->Qsh();
	sunlit_LAI = light->LAIsl();
	shaded_LAI = light->LAIsh();

	do {
		sunlit->SetVal(sunlit_PFD, info, NitrogenStatus.Nitrogendeficiencyone, weather.airT, weather.CO2, weather.RH,weather.wind, atmPressure, leafwidth, poolmod, PhotosyntheticFeedback, leafageEffect, weather.psil_, weather.ET_supply);
		shaded->SetVal(shaded_PFD, info, NitrogenStatus.Nitrogendeficiencyone, weather.airT, weather.CO2, weather.RH,weather.wind, atmPressure, leafwidth, poolmod, PhotosyntheticFeedback, leafageEffect, weather.psil_, weather.ET_supply);

		if (develop->get_istage()!=1) //post TI
		{
			if (initInfo.bigleaf == 0)
			{
				CGas_exchange_new* sunlit2 = new CGas_exchange_new();
				CGas_exchange_new * shaded2 = new CGas_exchange_new();
				sunlit2->SetVal(sunlit_PFD, info, NitrogenStatus.Nitrogendeficiencyone, weather.airT, weather.CO2, weather.RH,weather.wind, atmPressure, leafwidth, 2., PhotosyntheticFeedback, leafageEffect, weather.psil_, weather.ET_supply);
				shaded2->SetVal(shaded_PFD, info, NitrogenStatus.Nitrogendeficiencyone, weather.airT, weather.CO2, weather.RH,weather.wind, atmPressure, leafwidth, 2., PhotosyntheticFeedback, leafageEffect, weather.psil_, weather.ET_supply);
				photosynthesis_gross = (sunlit->A_gross*sunlit_LAI + shaded->A_gross*shaded_LAI)*(1-Tubmod)+Tubmod*(sunlit2->A_gross*sunlit_LAI + shaded2->A_gross*shaded_LAI);
				photosynthesis_gross = (sunlit->A_gross*sunlit_LAI + shaded->A_gross*shaded_LAI);
				if ((sunlit_LAI+shaded_LAI)>0) reference = (sunlit->A_gross*sunlit_LAI*3./(sunlit_LAI+shaded_LAI) + shaded->A_gross*shaded_LAI*3./(sunlit_LAI+shaded_LAI))*(1-Tubmod)+Tubmod*(sunlit2->A_gross*sunlit_LAI*3./(sunlit_LAI+shaded_LAI) + shaded2->A_gross*shaded_LAI*3./(sunlit_LAI+shaded_LAI));
					else reference = 1.;
				reference = 1.;
				photosynthesis_net = (sunlit->A_net*sunlit_LAI + shaded->A_net*shaded_LAI)*(1-Tubmod)+Tubmod*(sunlit2->A_net*sunlit_LAI + shaded2->A_net*shaded_LAI);//
				photosynthesis_net = (sunlit->A_net*sunlit_LAI + shaded->A_net*shaded_LAI);
				transpiration = (sunlit->ET*sunlit_LAI + shaded->ET*shaded_LAI)*(1-Tubmod) + Tubmod*(sunlit2->ET*sunlit_LAI + shaded2->ET*shaded_LAI);//mmol h2o m-2 ground s-1
				//transpiration = sunlit->ET; //mmol h20 m-2 leaf s-1
				//temperature = (sunlit->Tleaf*sunlit_LAI + shaded->Tleaf*shaded_LAI)/LAI*(1-Tubmod)+Tubmod*(sunlit->Tleaf*sunlit_LAI + shaded->Tleaf*shaded_LAI)/LAI;
				temperature = (sunlit->Tleaf*sunlit_LAI + shaded->Tleaf*shaded_LAI)/LAI;
				set_sunPFD(sunlit_PFD);
				set_shadePFD(shaded_PFD);
				set_PFD(weather.PFD);
				set_sunPg(sunlit->A_gross*(1-Tubmod) + Tubmod*(sunlit2->A_gross));
				set_shadePg(shaded->A_gross*(1-Tubmod) + Tubmod*(shaded2->A_gross));
				set_VPD(sunlit->get_VPD());
			}
			else
			{
				photosynthesis_gross = (sunlit->A_gross*sunlit_LAI+shaded->A_gross*shaded_LAI);
				photosynthesis_net = (sunlit->A_net*sunlit_LAI+shaded->A_net*shaded_LAI);
				if (info.Water_stress_simulation_type==1) {
					photosynthesis_gross *= weather.swdf1;
					photosynthesis_net *= weather.swdf2;
				}
				transpiration = (sunlit->ET*sunlit_LAI+shaded->ET*shaded_LAI);
				temperature = (sunlit->Tleaf*sunlit_LAI+shaded->Tleaf*shaded_LAI)/LAI;
				set_sunPFD(sunlit_PFD);
				set_shadePFD(shaded_PFD);
				set_PFD(weather.PFD);
				set_sunPg(sunlit->A_gross);
				set_shadePg(shaded->A_gross);
			}
		}
		else //prior to TI, same routine for single and bigleaf canopy
		{
			photosynthesis_gross = (sunlit->A_gross*sunlit_LAI + shaded->A_gross*shaded_LAI);//umol co2 m-2 ground s-1
			if ((sunlit_LAI+shaded_LAI)>0) reference = (sunlit->A_gross*sunlit_LAI*3./(sunlit_LAI+shaded_LAI) + shaded->A_gross*shaded_LAI*3./(sunlit_LAI+shaded_LAI));
			else reference = 1.;
			photosynthesis_net = (sunlit->A_net*sunlit_LAI + shaded->A_net*shaded_LAI);//
		
			if (info.Water_stress_simulation_type == 1) {
				photosynthesis_gross *= weather.swdf1;
				photosynthesis_net *= weather.swdf2;
			}
			transpiration = (sunlit->ET*sunlit_LAI + shaded->ET*shaded_LAI);//mmol h2o m-2 ground s-1
			if (LAI != 0) 
			{
				temperature = (sunlit->Tleaf*sunlit_LAI + shaded->Tleaf*shaded_LAI)/LAI;
			}
			else 
			{
				temperature = sunlit->Tleaf;
			}
			set_sunPg(sunlit->A_gross);
			set_shadePg(shaded->A_gross);
		}
		// adjust for N deficiency based on SPUDSIM - now takes place in Gasexchange class reducing Jmax
		//photosynthesis_gross *= NitrogenStatus.Nitrogendeficiencyone;
 		//photosynthesis_net *= NitrogenStatus.Nitrogendeficiencyone;
		//transpiration *= NitrogenStatus.Nitrogendeficiencyone;
		if (LAI != 0)
		{
			this->conductance = __max(0.,((sunlit->get_gs() * sunlit_LAI + shaded->get_gs() * shaded_LAI)/LAI));//average stomatal conductance
			//this->conductance = (sunlit->get_gs() * sunlitLAI() + shaded->get_gs() * shadedLAI())/LAI);//average stomatal conductance
			this->internal_CO2 = __max(0.,((sunlit->get_ci() * sunlit_LAI + shaded->get_ci()*shaded_LAI)/LAI)); //average internal CO2 concentrtion
		}
		else this->conductance = 0.;
		assimilate = (photosynthesis_gross*CO2_MW/1.0e6)*(60.0*initInfo.timeStep)/initInfo.plantDensity;//g CO2 plant-1 inc-1

		/* Limit photosynthetic rate to twice growh demand from last time-step
		/* Otherwise, we get issue with too much C_pool build up */
		if ((shootPart_old + rootPart_old + tubPart_old)>0){
			if(assimilate*1*CH2O_MW/CO2_MW > 2 * (shootPart_old + rootPart_old + tubPart_old)/Rgorgan) {
 				redo = 1;
				//assimilate = 1.2 * (shootPart_old + rootPart_old + tubPart_old)/Rgorgan/(1*CH2O_MW/CO2_MW);
				poolmod *= 0.98;
				if (poolmod < 0.35) {
					redo = 0;
				}
			}
			else redo = 0;
		} 
		else redo = 0;
	} while (redo == 1);
	assimilate = (photosynthesis_gross*CO2_MW/1.0e6)*(60.0*initInfo.timeStep)/initInfo.plantDensity;//g CO2 plant-1 inc-1
	assimilate_old = assimilate;
	reference = (reference*CO2_MW/1.0e6)*(60.0*initInfo.timeStep)/initInfo.plantDensity;//g CO2 plant-1 inc-1
	if (reference <= 0) reference = 0.01;
	
	/*Following 4 variables are for reporting purposes to output file, ET is for instanteous value for 2DSOIL*/
		//Pnet+=(photosynthesis_net/1.0e6*(60*initInfo.timeStep)); //cumulative mol CO2 m-2 ground time-1
		//Pgross+=(photosynthesis_gross/1.0e6*(60*initInfo.timeStep)); //cumulative mol CO2 m-2 ground time-1
		//ET+=(transpiration/1000*(60/initInfo.timeStep)); //cumulative mmol H2O m-2 ground time-1
	Pnet+=(photosynthesis_net*30/1.0e6*(60.*initInfo.timeStep))/initInfo.plantDensity; //cumulative g CH2O plant-1 time-1
	Pgross+=(photosynthesis_gross*30/1.0e6*(60.*initInfo.timeStep))/initInfo.plantDensity;// "
	cumulativeET+=(transpiration*18/1.0e3*(60.*initInfo.timeStep))/initInfo.plantDensity;//cumulative g h2o plant-1 inc-1
//	set_absPAR(sunlit_PFD*.8+shaded_PFD*.8);
	set_sunPFD(sunlit_PFD);
	set_shadePFD(shaded_PFD);
	set_PFD(weather.PFD);
	set_SRAD(weather.solRad);
	set_sunLAI(sunlit_LAI);
	set_shadeLAI(shaded_LAI);

	photosynthesis_netsunlitleaf = sunlit->A_net;
	photosynthesis_netshadedleaf = shaded->A_net;
	photosynthesis_grosssunlitleaf = sunlit->A_gross;
	photosynthesis_grossshadedleaf = shaded->A_gross;
	psistress_gs_factor = sunlit->get_psileaf_Pn_stress();

	transpiration_sunlitleaf = sunlit->ET;
	transpiration_shadedleaf = shaded->ET;
	if ((weather.PFD >= 350) && (weather.PFD <= 450))
	{
		if (photosynthesis_gross > photosynthesis_gross900) photosynthesis_gross900 = photosynthesis_gross;
		if (photosynthesis_net > photosynthesis_net900) photosynthesis_net900 = photosynthesis_net;
		if (sunlit->A_gross > photosynthesis_grossleaf900) photosynthesis_grossleaf900 = sunlit->A_gross;
	}
	delete sunlit;
	delete shaded;
	delete sun;
	delete light;
}

void CPlant::calcLeafArea(const TInitInfo info) 
{
	/***********************************************/
	//Calculate green leaf, all leaf, and dead leaf area at beginning of each time-step
	//Note: information on each senesced leaf is still maintained in the appropriate nodalunit array on the plant
	//Note: bigleaf = 0 is for individual organ model, = 1 is for bigleaf, SIMGUI based leaf model
	//		-for SIMGUI approach, all leaf area and mass information stored in NodalUnit[0] only
	/***********************************************/
	double garea = 0.0;
	double darea = 0.0;
	if (info.bigleaf == 0)
	{
		if (inodeNumber > 0)
		{
			for (int i = 0; i < inodeNumber; i++)
			{
				if (!nodalUnit[i].get_leaf()->isTerminated()) garea += nodalUnit[i].get_leaf()->get_greenArea();
				else
					darea += nodalUnit[i].get_leaf()->get_area();
			}
		}
	}
	if (info.bigleaf == 1 && !nodalUnit=='\0')
	{
		garea = nodalUnit[0].get_leaf()->get_greenArea();
		darea = nodalUnit[0].get_leaf()->get_area()-nodalUnit[0].get_leaf()->get_greenArea();
	}
	else if (info.bigleaf ==1) garea = darea = 0.;
	set_leafArea(garea+darea);
	set_LAI((garea+darea)*info.plantDensity/10000.);
	set_greenLeafArea(garea);
	set_greenLAI(garea*info.plantDensity/10000.);
	set_deadLeafArea(darea);
	set_deadLAI(darea*info.plantDensity/10000.);
	if (greenLAI > maxLAI) maxLAI = greenLAI;
	return;
}

double CPlant::calcPotentialLeafArea()
{
	double area = 0.0;
	for (int i = 0; i <= inodeNumber; i++)
	{
	//	nodalUnit[i].get_leaf->set_potentialArea();
	//	area += nodalUnit[i].get_leaf()->get_potentialArea();
	}
	potentialLeafArea = area;
	return area;
}

void CPlant::calcMaintRespiration(const TWeather & w, const TInitInfo info)
/********************************************************/
/* Simulate maintenance respiration costs for each organ type
/* //DHF: new routine is combination of SK (with q10fn) but modified to estimate Rm for each individual organ
//  after Ng and Loomis (1984_ routines and maintenance coefficients are under the leaf or organ classes
//   the leaf clas uses a higher respiration coefficient)
//	 also use Soo's 'agefn' that decreases respiration rate of canopy based on ratio of current green leaf to total leaf fixed over season
// SK:based on McCree's paradigm, See McCree(1988), Amthor (2000), Goudriaan and van Laar (1994)
// units very important here, be explicit whether dealing with gC, gCH2O, or gCO2
/* Note
/*	Big leaf model usage folows same approach, only single nodalunit is used*/

/*******************************************************/

{
	const double Q10 = 2.0; // SK: typical Q10 value for respiration, Loomis and Amthor (1999) Crop Sci 39:1584-1596
	double dt = info.timeStep/(24.*60.);
	double Rleaf=0., Rstem=0., Rroot=0., Rtuber=0., Rpool = 0.;
	//SK const double maintCoeff = 0.015; // gCH2O g-1DM day-1 at 20C for young plants, Goudriaan and van Laar (1994) Wageningen textbook p 54, 60-61	
	double q10fn = pow(Q10,(w.airT - 25.0)/10.);
	double agefn = (greenLeafArea+1) / (leafArea+1);
	if (initInfo.bigleaf == 0)
	{
		/*Rm fraction due to leaf and stem mass*/
		if (inodeNumber == 0) 
		{
			maintRespiration = 0.;
			instantResp = 0.;
			maintRespleaf = 0.;
			maintRespstem = 0.;
			maintResptuber = 0.;
			maintResproot = 0.;
			return;
		}
		for (int i=0; i < inodeNumber; i++)
		{
			if ((nodalUnit[i].isInitiated())&&(!nodalUnit[i].isTerminated()))
			{
				Rleaf += nodalUnit[i].get_leaf()->get_maintResp(info);
				Rstem += nodalUnit[i].get_stem()->get_maintResp();
				//nodalUnit[i].get_leaf()->import_CH2O(-nodalUnit[i].get_leaf()->get_maintResp()*dt*q10fn);
				//nodalUnit[i].get_leaf()->import_CH2O(-nodalUnit[i].get_stem()->get_maintResp()*dt*q10fn);
			}		
		}
		/*Rm fractino due to tubers and roots mass*/
		//Rleaf = Rleaf * __max(1,assimilate/reference); //adjust for metabolic activity
		if (this->get_tubers()->isInitiated()) Rtuber = this->get_tubers()->get_maintResp(this->age-develop->dTuberlinear());
		if (this->get_roots()->get_drymass() > 0) Rroot = this->get_roots()->get_maintResp();
		if (C_pool > 0) Rpool = nodalUnit[0].get_stem()->get_maintResp()/nodalUnit[0].get_stem()->get_drymass()*C_pool; //maintresp for soluble C pool, but don't penalize for seed reserve
		//maintRespiration = agefn*q10fn*maintCoeff*totalMass*dt;// gCH2O plant-1 dt-1
		//if(this->get_tubers()->isInitiated()) this->get_tubers()->import_CH2O(-Rtuber*dt*q10fn);
		//this->get_roots()->import_CH2O(-Rroot*dt*q10fn);
		//if (this->get_stemMass() < (C_reserve+C_pool)) Rstem = 0; // don't penalize for seed reserve supporting stem/leaf growth
		//if (this->get_leafMass() < (C_reserve+C_pool)) Rleaf = 0;
		//if (this->get_rootMass() < (C_reserve+C_pool)) Rroot = 0;
	}
	else
	{
		if (!nodalUnit=='\0'){
			Rleaf = nodalUnit[0].get_leaf()->get_maintResp(info);
			Rstem = nodalUnit[0].get_stem()->get_maintResp();
			if (this->get_tubers()->isInitiated()) Rtuber = this->get_tubers()->get_maintResp(this->age-develop->dTuberlinear());
			if (this->get_roots()->get_drymass() > 0) Rroot = this->get_roots()->get_maintResp();
			if (C_pool > 0) {
				Rpool = nodalUnit[0].get_stem()->get_maintResp()/nodalUnit[0].get_stem()->get_drymass()*C_pool; //maintresp for soluble C pool, but don't penalize for seed reserve
			}
		}
		else
		{
			maintRespiration = 0.;
			instantResp = 0.;
			maintRespleaf = 0.;
			maintRespstem = 0.;
			maintResptuber = 0.;
			maintResproot = 0.;
			return;
		}
	}
	maintRespiration = agefn*q10fn*(Rleaf+Rstem+Rroot+Rtuber+Rpool)*dt; //*gets subtracted from C fraction partitioned to each organ pool
	Rm += maintRespiration; //cumulative main respiration in g ch2o plant-1 dt-1
	instantResp = maintRespiration / 30. * 1000000. * info.plantDensity / dt  / (24.*3600.); //instant resp'n rate, umol co2 m-2 ground s-1
	maintRespleaf = agefn*q10fn*Rleaf*dt;
	maintRespstem = agefn*q10fn*Rstem*dt;
	maintResptuber = agefn*q10fn*Rtuber*dt;
	maintResproot = agefn*q10fn*Rroot*dt;
}
 
void CPlant::dailyave(int iCur, const TWeather& weather, const TInitInfo info, double dayinc)
{
	/********************************************/
	/* Routine provides 24 h averaged T, PAR, Daylength data
	/* Required by legacy SIMPOTATO phenology routines and leaf expansion routine
	/********************************************/

	//first reset variables at beginning of day
	this->PFD = weather.PFD;
	this->SRAD = weather.solRad;
	//this->temperature = weather.airT;
	if (fmod(iCur,((24.*60.)/initInfo.timeStep))==0) //first hour of day
	{
		Tdaylag = Tdayave;
		Tdaymax = Tdayave = Parave = Sradave = Photoave = temperature = 0;
		Tcanopyave = 0.;
		Pgross = 0.; Pnet = 0.; Rm= 0.; Rg = 0.;
		photosynthesis_gross900 = photosynthesis_net900 = photosynthesis_grossleaf900 = 0.;
		cumulativeET = 0.;
		Tdaymin = weather.airT; 
		Total_soil_water_uptake = 0.;
		average_LWP = 0.;
		average_gs = 0.;
		psistress_gs_factor_24h = 0.;
		Nstressfactorone_24h = 0.;
		Nstressfactortwo_24h = 0.;
		psi_leafexpansion_stress_24h = 0.;
		psi_leafexpansion_stress = 1.;
		if (iCur == 0) Tdaylag = weather.airT;
		if (develop->Emerged()==true)
		{
			age+=1;
		}
	}
	if (weather.time >= 0) { //then increment
		if (weather.airT > Tdaymax) Tdaymax = weather.airT;
		if (weather.airT < Tdaymin) Tdaymin = weather.airT;
		Tdayave += weather.airT;
		Tcanopyave += temperature;
		if (weather.PFD > 5) Photoave += info.timeStep;
		Parave += (weather.PFD * info.timeStep * 60.) / 1000000.; //mol m-2 timeinc-1
		Sradave += (weather.PFD * info.timeStep * 60. / 4.57 * 2. / 1000000.); // MJ srad m-2 timeinc-1
		average_LWP += weather.psil_;
		average_gs += this->get_conductance();
		psistress_gs_factor_24h += this->get_psistress_gs_factor();
		Nstressfactorone_24h += this->get_Nstressfactorone();
		Nstressfactortwo_24h += this->get_Nstressfactortwo();
		if (nodalUnit != NULL) psi_leafexpansion_stress_24h += this->get_nodalUnit()->get_leaf()->get_psi_leafexpansion_stress();
		if (nodalUnit != NULL) psi_leafexpansion_stress = this->get_nodalUnit()->get_leaf()->get_psi_leafexpansion_stress();
		//if (initInfo.bigleaf == 0) psi_leafexpansion_stress_24h += nodalUnit[i].get_leaf()->get_psi_leafexpansion_stress();
		//get_psi_leafexpansion_stress

		//plant->get_nodalUnit()->get_leaf()->get_psi_leafexpansion_stress()
	}
	if (fmod(iCur+1,((24.*60.)/initInfo.timeStep)) < 0.001)//i.e., last timestep prior to 24:00 or 0:00 next day
	//if (!(iCur==0)&&fmod(iCur,(24*60)/initInfo.timeStep)==0)
	{
		Parave += (weather.PFD * info.timeStep * 60.)/1000000.; //mol m-2 timeinc-1
		Sradave += (weather.PFD * info.timeStep * 60. / 4.57 * 2. / 1000000.); // MJ srad m-2 timeinc-1
		Tdayave = (Tdayave) / (24.*60./info.timeStep);
		Tcanopyave = (Tcanopyave) / (24.*60./info.timeStep);
		//Photoave = (Photoave+info.timeStep) / 60; // convert from minutes to hours
		Photoave = weather.dayLength; //override prior calc - more stable to use estimate from celestial calcs in 2DSOIL
		average_LWP = average_LWP / (24.*60./info.timeStep);
		average_gs = average_gs / (24.*60./info.timeStep);
		psistress_gs_factor_24h = psistress_gs_factor_24h / (24. * 60. / info.timeStep);
		Nstressfactorone_24h = Nstressfactorone_24h / (24. * 60. / info.timeStep);
		Nstressfactortwo_24h = Nstressfactortwo_24h / (24. * 60. / info.timeStep);


		psi_leafexpansion_stress_24h = psi_leafexpansion_stress_24h / (24. * 60. / info.timeStep);
		if (psi_leafexpansion_stress_24h < 0) psi_leafexpansion_stress_24h = 0.;
		if (psi_leafexpansion_stress < 0) psi_leafexpansion_stress = 1.; 
		this->set_TranspirationSupplyDemandRatio(this->get_totalSoilWaterUptake()/this->get_cumulativeET());
		swdf1 = this->get_transpirationSupplyDemandRatio();
		
		if (initInfo.Water_stress_off==1) swdf1 =1.;
	}
}

void CPlant::initiate(int iCur, const TWeather& weather)
{
	/*********************************************************/
	//Subroutine initiates new nodalunits according to leaf appearance rate estimated in develop method*/
	//Tubers will also initiate if called for from results in development->update()
	/* NodalUnit initiation dependent on C reserve status, N status, plus leaf appearance rate
	// Routine is primarily of importance to Bigleaf = 0 where individual nodal units are simulated

	// Current assumptions / limitations - Nov 2008:
	// -4 new nodalunits at plant emergence
	// -need 0.0525g of C_reserve to form new nodalunit and initial rootmass
	// - initial N also needed before leaf can form which is baed on percentage of initial nodal unit mass

	// -leaf appearance rate set to 0 after new leaf forms for basal and apical stems, or reduced by 1 for mainstem
	// -flowering occurs at 16 mainstem nodes at whichpoint 2 apicals along with nodalunits, form
	// -a single basal branch can form after 10 nodes exist on mainstem

	// If BigLeaf Model Used
	//	-much simplified
	//	-ignores initiation of multiple nodal units, branches, flowering etc.
	//	-Should add this functionality back and consider linking organs with overall canopy area estimates?
	/**********************************************************/

	int newnodes, stemnodes, mainstemnodes, basalbrno, flag, temp;
	newnodes = stemnodes = mainstemnodes = basalbrno = flag = temp =0;
	double nodemass, stemmass, leafmass;
	double Nreq, Nbuffer;
	double valve = 0.2 / (60./initInfo.timeStep); //limit CHO availability of reserve for new growth*/
	//double valve = 0.2 / 24; //assumes 20% of Creserve available for mobilization each day
	//double valve = 1; //no valve limitation
	nodemass = 0.0525/(60.*24./initInfo.timeStep); //initial node mass for leaf+stem at emergence - later in code, exact values are used
	if (nodalUnit != NULL) //establish N requirement for a new nodal unit to appear
	{
		Nreq = nodalUnit[0].get_leaf()->get_drymass_at_emergence()*NitrogenStatus.optimumleafNitrogenRatio + nodalUnit[0].get_stem()->get_drymass_at_emergence() * NitrogenStatus.optimumstemNitrogenRatio;
		Nbuffer = this->leafMass * (NitrogenStatus.actualleafNitrogenRatio - NitrogenStatus.minimumleafNitrogenRatio);
	}
	else
	{
		Nreq = 0.; 
		Nbuffer = 0.;
	}
	if (this->get_greenLeafArea() < 400) Nbuffer += N_reserve; //temporarily add N reserve to support new leaf growth when plant is small
	int Nredistribute = 0; // if new nodal unit, need to redistribute N to prevent mass balance issue
	
	
	//1 Is it tuber initiation?
	if ((develop->Tuberlinear()) && !(tubers->isInitiated()))
	{
		tubers->setInitiated();
	}

	//2 If plant starts out at emergence, establish initial root mass and N content
	if (((develop->Emerged()) && (this->get_rootMass() <= 0.)))	//  Note: to keep consistent with maizesim, root initial mass is now set in crop.cpp class at sowing date
	{
		//roots->set_EmergenceData(nodemass*4./2.); /*let initial rootmass equal to initial above ground stem mass*/ - see note above 12-15-10
		//C_seed = C_seed - this->get_rootMass();
		//N_seed -= this->get_roots()->get_currentNitrogenAmount();
	}

	if (initInfo.bigleaf == 1) //big leaf model
	{
		if ((develop->Vegetative()==true) && (nodalUnit == NULL))
		{
			nodalUnit = new CNodalUnit[1];
			inodeNumber = int(develop->get_LvsAppeared());
			nodalUnit[0].initialize(initInfo,0,0,0,0,develop);
			C_seed -= (nodalUnit[0].get_leaf()->get_drymass()+nodalUnit[0].get_stem()->get_drymass());
			N_reserve -= (nodalUnit[0].get_leaf()->get_currentNitrogenAmount()+nodalUnit[0].get_stem()->get_currentNitrogenAmount());
			this->set_N_seedused(nodalUnit[0].get_leaf()->get_currentNitrogenAmount()+nodalUnit[0].get_stem()->get_currentNitrogenAmount());
			NitrogenStatus.leafNitrogenAmount = nodalUnit[0].get_leaf()->get_currentNitrogenAmount();
			NitrogenStatus.stemNitrogenAmount = nodalUnit[0].get_stem()->get_currentNitrogenAmount();
			NitrogenStatus.rootNitrogenAmount = this->get_roots()->get_currentNitrogenAmount();
		}
	}

	if (initInfo.bigleaf == 0)
	{
		// case 1: //create initial leaves at emergence
		if ((develop->Vegetative()==true) && (nodalUnit == NULL)) 
		{
			nodalUnit = new CNodalUnit[2500]; // create enough leaf nodes for now, to be replaced by dynamic collection
			inodeNumber = int(develop->get_LvsAppeared());
			for (int i=0; i < inodeNumber; i++) //start with 4 new leaves on mainstem
			{
				nodalUnit[i].initialize(initInfo, i,0,0,i,develop);
				igreennodeNumber += 1;
				iyoungnodeNumber += 1;
				//need to subtract mass of leaf and stem from seedpiece C reserve
				C_seed = C_seed - nodalUnit[i].get_leaf()->get_drymass()-nodalUnit[i].get_stem()->get_drymass();
				C_seed_used += nodalUnit[i].get_leaf()->get_drymass()+nodalUnit[i].get_stem()->get_drymass();
				N_reserve -= (nodalUnit[i].get_leaf()->get_currentNitrogenAmount() + nodalUnit[i].get_stem()->get_currentNitrogenAmount());
				this->set_N_seedused(nodalUnit[i].get_leaf()->get_currentNitrogenAmount() + nodalUnit[i].get_stem()->get_currentNitrogenAmount());
				NitrogenStatus.leafNitrogenAmount += nodalUnit[i].get_leaf()->get_currentNitrogenAmount();// N for new stem, leaf is from Nreserve when those organs are instantiated at DOE 0 - no need to account for this!
				NitrogenStatus.stemNitrogenAmount += nodalUnit[i].get_stem()->get_currentNitrogenAmount();
			}
		}

		// case 2: //initialize new leaves if called for on appropriate stem at appropriate location
		temp = inodeNumber;
		if ((Bnewnode == true) && (nodalUnit != NULL))
		{
			leafmass = 0.0005; //these two values are initial leaf and stem mass requirements of a nodal unit at emergence - shouldn't be hardwired?
			stemmass = 0.0025;
			nodemass = leafmass + stemmass;
			Nreq = leafmass * NitrogenStatus.optimumleafNitrogenRatio + stemmass * NitrogenStatus.optimumstemNitrogenRatio;
			Nbuffer = this->leafMass * (NitrogenStatus.actualleafNitrogenRatio - NitrogenStatus.minimumleafNitrogenRatio);

			for (int i = 0; i < inodeNumber; i++)
			{
				stemnodes = 0;
				if (nodalUnit[i].get_LAR()>1)
				{
					if (nodalUnit[i].get_type()==0 && !develop->Flowered() && nodalUnit[i].get_node() == 0) 
					{//mainstem node addition - don't allow if mainstem has flowered
						if ((C_seed+C_pool)*valve < nodemass) return; //Check 1: not enough CHO reserve to support new leaf
						if (Nbuffer < Nreq || NitrogenStatus.stemNitrogenAmount < stemmass * NitrogenStatus.optimumstemNitrogenRatio) return; //Check 2: not enough N supply to support new leaf
						//exit routine if plant can't acess seed N reserve and there's not enough N in leaves for nodal unit
						for (int j = 0; j< inodeNumber; j++) //find highest node position
						{
							if (nodalUnit[j].get_type()==0 && nodalUnit[j].get_node() > stemnodes) stemnodes = nodalUnit[j].get_node();
						}
						flag = flag + 1;
						nodalUnit[temp].initialize(initInfo,temp,0,0,stemnodes+1,develop);
						
						if (C_seed > nodemass) {
							C_seed = C_seed - nodalUnit[temp].get_leaf()->get_drymass()-nodalUnit[temp].get_stem()->get_drymass();
							C_seed_used += nodalUnit[temp].get_leaf()->get_drymass()+nodalUnit[temp].get_stem()->get_drymass();
						}
						else if (C_pool > nodemass){
							C_pool = C_pool - nodalUnit[temp].get_leaf()->get_drymass()-nodalUnit[temp].get_stem()->get_drymass();
							C_pool_used += nodalUnit[temp].get_leaf()->get_drymass()-nodalUnit[temp].get_stem()->get_drymass();
						}
						else{
							C_seed += C_pool;
							C_seed = C_seed - nodalUnit[temp].get_leaf()->get_drymass()-nodalUnit[temp].get_stem()->get_drymass();
							C_pool_used += C_pool;
							C_seed_used += (nodalUnit[temp].get_leaf()->get_drymass()+nodalUnit[temp].get_stem()->get_drymass() - C_pool_used); 
							C_pool = 0;
						}
						Nredistribute = 1;
						newleafNitrogen += nodalUnit[temp].get_leaf()->get_currentNitrogenAmount();
						newstemNitrogen += nodalUnit[temp].get_stem()->get_currentNitrogenAmount();
						temp+=1;
						nodalUnit[i].set_LAR();//reset LAR of node to 0
					}
					if (nodalUnit[i].get_type()==1 && nodalUnit[i].get_node()==0)
					{//basal stem node addition
						stemnodes = 0;
						if ((C_pool*valve) < nodemass) return; //Check 1: Not enough CHO available 
						//if (C_new_organ < nodemass) return;					
						if (Nbuffer < Nreq || NitrogenStatus.stemNitrogenAmount < stemmass * NitrogenStatus.optimumstemNitrogenRatio) return; //not enough N to form new nodalunit
						for (int j = 0; j< inodeNumber; j++)
						{
							if (nodalUnit[j].get_type()==1 && nodalUnit[j].get_node() > stemnodes) stemnodes = nodalUnit[j].get_node();
						}			
						flag = flag + 1;
						nodalUnit[temp].initialize(initInfo, temp,1,nodalUnit[i].get_location(),stemnodes+1, develop);
						C_pool = C_pool - nodalUnit[temp].get_leaf()->get_drymass()-nodalUnit[temp].get_stem()->get_drymass();
						C_pool_used += 	C_seed_used += nodalUnit[temp].get_leaf()->get_drymass()+nodalUnit[temp].get_stem()->get_drymass();
						Nredistribute = 1;
						newleafNitrogen += nodalUnit[inodeNumber].get_leaf()->get_currentNitrogenAmount();
						newstemNitrogen += nodalUnit[inodeNumber].get_stem()->get_currentNitrogenAmount();
						temp+=1;
						nodalUnit[i].set_LAR();
					}
					if (nodalUnit[i].get_type()==2 && nodalUnit[i].get_node()==0)
					{//apical stem node addition
						stemnodes = 0;
						if ((C_pool*valve) < nodemass) return;
						//if (C_new_organ < nodemass) return;
						if (Nbuffer < Nreq || NitrogenStatus.stemNitrogenAmount < stemmass * NitrogenStatus.optimumstemNitrogenRatio) return; //not enough N to form new nodalunit
						for (int j = 0; j< inodeNumber; j++)
						{
							if (nodalUnit[j].get_type()==2 && nodalUnit[j].get_node() > stemnodes) stemnodes = nodalUnit[j].get_node();
						}
						flag = flag + 1;
						nodalUnit[temp].initialize(initInfo, temp,2,nodalUnit[i].get_location(),stemnodes+1,develop);
						C_pool = C_pool - nodalUnit[temp].get_leaf()->get_drymass()-nodalUnit[temp].get_stem()->get_drymass();
						C_pool_used += nodalUnit[temp].get_leaf()->get_drymass()-nodalUnit[temp].get_stem()->get_drymass();
						Nredistribute = 1;
						newleafNitrogen += nodalUnit[temp].get_leaf()->get_currentNitrogenAmount();
						newstemNitrogen += nodalUnit[temp].get_stem()->get_currentNitrogenAmount();
						temp+=1;
						nodalUnit[i].set_LAR();
					}
				}
			}
		}
		develop->set_LvsAppeared(flag);
		inodeNumber += flag; //update total number of nodes on the plant
		igreennodeNumber += flag;
		flag = 0;
		// case 3: determine if new lateral branches can form
		// for version 1.0 very simple rules:
		// - 2 apical branches form at mainstem flowering
		// - any number of basal branches can form at lag of 10 from mainstem nodes
		// - mainstems have 17 nodes prior to flowering
		// - unlimited number of nodes may appear on stems (except mainstem), but expansion and appearance
		//   should be limited by assimilate and nitrogen availability 
		mainstemnodes = 0;
		for (int i = 0; i < inodeNumber; i++)
		{
			//get # ms nodes and # basal branches
			if (nodalUnit[i].get_type()==0) mainstemnodes += 1;
			if ((nodalUnit[i].get_type()==1) && (nodalUnit[i].get_node() ==0)) basalbrno += 1;
		}

		//if ((mainstemnodes ==16 && !develop->Flowered())) develop->set_Flowered(weather.jday); //*2009 SPUDSIM paper
		/*2009 SPUDSIM Paper - Remove branch formation*/
		/* Two apical branches automatically form at flowering, 17 mainstem nodes */
		if ((mainstemnodes == 16 && !develop->Flowered())){ //assume plant has flowered at this stage and 2 apical branches can initiate nodes
			develop->set_Flowered(weather.jday);
			if ((C_seed+C_pool)*valve < 2.*nodemass) return; //Check 1 - not enough CHO to support new branch formation, factor of 2 due to two new apical branches
			//if (C_new_organ < (2*nodemass)) return;
			if (Nbuffer < 2.* Nreq || NitrogenStatus.stemNitrogenAmount < 2.*stemmass * NitrogenStatus.optimumstemNitrogenRatio) return; //Check 2 - not enough N to support new branch formation - factor of 2 required because 2 branches are being formed	
			nodalUnit[inodeNumber].initialize(initInfo,inodeNumber,2,0,0, develop);			
			C_pool = C_pool - nodalUnit[inodeNumber].get_leaf()->get_drymass()-nodalUnit[inodeNumber].get_stem()->get_drymass();
			C_pool_used += nodalUnit[inodeNumber].get_leaf()->get_drymass()-nodalUnit[inodeNumber].get_stem()->get_drymass();
			Nredistribute = 1;
			newleafNitrogen += nodalUnit[inodeNumber].get_leaf()->get_currentNitrogenAmount();
			newstemNitrogen += nodalUnit[inodeNumber].get_stem()->get_currentNitrogenAmount();
			nodalUnit[inodeNumber+1].initialize(initInfo,inodeNumber+1,2,1,0, develop);
			C_pool = C_pool - nodalUnit[inodeNumber+1].get_leaf()->get_drymass()-nodalUnit[inodeNumber+1].get_stem()->get_drymass();
			C_pool_used += nodalUnit[inodeNumber+1].get_leaf()->get_drymass()-nodalUnit[inodeNumber+1].get_stem()->get_drymass();
			Nredistribute = 1;
			newleafNitrogen += nodalUnit[inodeNumber].get_leaf()->get_currentNitrogenAmount();
			newstemNitrogen += nodalUnit[inodeNumber].get_stem()->get_currentNitrogenAmount();
			develop->set_LvsAppeared(2);
			inodeNumber +=2;
		}
		/* Basal lateral branches can form based on accumulated temperature, CHO, and N availability, but 1 at a time*/
		if ((nodalUnit[0].get_BAR()/(basalbrno+1)) > 1) //add new basalbranch, but increase time in between appearance of each one
		{
			if ((C_pool*valve) < nodemass) return; //Check 1 - not enough CHO to support new branch formation, factor of 2 due to two new apical branches
			//if (C_new_organ < nodemass) return;
			if (Nbuffer < Nreq || NitrogenStatus.stemNitrogenAmount < stemmass * NitrogenStatus.optimumstemNitrogenRatio) return; //Check 2 - not enough N to support new branch formation - factor of 2 required because 2 branches are being formed	
			if ((mainstemnodes - 9) <= 0) return; // Check 3 - induce an apical dominance type lag via the mainstem
			
			
			if ((basalbrno+3)>mainstemnodes) return;//restrict basal branch numbers to total number of mainstemnodes - 2 apical branches
			
			
			nodalUnit[inodeNumber].initialize(initInfo,inodeNumber,1,basalbrno+1,0,develop);
			C_pool = C_pool - nodalUnit[inodeNumber].get_leaf()->get_drymass()-nodalUnit[inodeNumber].get_stem()->get_drymass();
			C_pool_used += 	C_seed_used += nodalUnit[inodeNumber].get_leaf()->get_drymass()+nodalUnit[inodeNumber].get_stem()->get_drymass();
			Nredistribute = 1;
			newleafNitrogen += nodalUnit[inodeNumber].get_leaf()->get_currentNitrogenAmount();
			newstemNitrogen += nodalUnit[inodeNumber].get_stem()->get_currentNitrogenAmount();
			inodeNumber += 1;
			develop->set_LvsAppeared(1);
			nodalUnit[0].set_BAR();
		}
		
		if (Nredistribute == 1) 
			// When new nodal unit is formed, N must be taken either from existing nodalunits or from N_reserve
			// At this point in the time-stepping, N has already been allocated to each individual leaf and stem
			// This routine repartitions this N again among the organs.
			// Note that the whole plant leaf and stem indices, NitrogenStatus.leafNitrogenAmount and .stemNitrogenAmount 
			//	should not be affected unless N is being taken from the seed reserve
		{
			double tempNfromreserve = 0.;
			if (this->get_greenLeafArea()<400)
			{
				double tempNfromreserve = newleafNitrogen+newstemNitrogen;
				if (tempNfromreserve > N_reserve)
				{
					tempNfromreserve -= N_reserve; 
					this->set_N_seedused(N_reserve);
					N_reserve = 0.;
				}
				else
				{
					N_reserve -= tempNfromreserve;
					this->set_N_seedused(tempNfromreserve);
					NitrogenStatus.leafNitrogenAmount += newleafNitrogen; NitrogenStatus.stemNitrogenAmount += newstemNitrogen;
					return; // enough Nreserve available to fill new N demand, so no need to reallocate N among existing nodes
				}					
			}
			double tempNleafamount = 0.; double tempNstemamount = 0.; double minNreq = 0.; double deltaNreq = 0.; double luxNreq = 0.;
			for (int i = inodeNumber - 1; i >= 0; i--)
			{
				tempNleafamount += nodalUnit[i].get_leaf()->get_currentNitrogenAmount();
				tempNstemamount += nodalUnit[i].get_stem()->get_currentNitrogenAmount();
			}
			tempNleafamount -= newleafNitrogen; tempNstemamount -= newstemNitrogen; // don't double count N to recently initiated organs
			if ((this->get_leafArea() < 400) && (tempNfromreserve > 0)) //N reserve used, but not enough to satisfy new growth, so add portion used
			{
				double ratio = newleafNitrogen / newstemNitrogen;
				tempNleafamount += ratio * ((newleafNitrogen + newstemNitrogen) - tempNfromreserve);
				tempNstemamount += (1 - ratio) * ((newleafNitrogen + newstemNitrogen) - tempNfromreserve);
			}
			for (int i = inodeNumber - 1; i >= 0; i--) //First, each leaf gets minimum N content
			{	
				minNreq = nodalUnit[i].get_leaf()->get_drymass() * NitrogenStatus.minimumleafNitrogenRatio;
				if (minNreq > tempNleafamount) minNreq = tempNleafamount;
				nodalUnit[i].get_leaf()->set_currentNitrogenAmount(minNreq);
				tempNleafamount -= minNreq;
			}
			if (tempNleafamount > 0) //Second, starting with youngest leaf, N is added to optimum level
			{
				for (int i = inodeNumber - 1; i>=0; i--)
				{
					minNreq = nodalUnit[i].get_leaf()->get_drymass() * NitrogenStatus.minimumleafNitrogenRatio;
					deltaNreq = nodalUnit[i].get_leaf()->get_drymass() * (NitrogenStatus.optimumleafNitrogenRatio - NitrogenStatus.minimumleafNitrogenRatio);
					if (deltaNreq > tempNleafamount) deltaNreq = tempNleafamount;
					nodalUnit[i].get_leaf()->set_currentNitrogenAmount(deltaNreq+minNreq);
					tempNleafamount -= deltaNreq;
				}
			}
			if (tempNleafamount > 0) //Third, luxury N proportioned according to leaf N content
			{
				for (int i = inodeNumber - 1; i>-0; i--)
				{
					minNreq = nodalUnit[i].get_leaf()->get_drymass() * NitrogenStatus.minimumleafNitrogenRatio;
					deltaNreq = nodalUnit[i].get_leaf()->get_drymass() * (NitrogenStatus.optimumleafNitrogenRatio - NitrogenStatus.minimumleafNitrogenRatio);
					luxNreq = (minNreq + deltaNreq) * 0.9;
					if (luxNreq > tempNleafamount) luxNreq = tempNleafamount;
					nodalUnit[i].get_leaf()->set_currentNitrogenAmount(minNreq+deltaNreq+luxNreq);
					tempNleafamount -= luxNreq;
				}
				if (tempNleafamount > 0) //if still extra N sitting around, store it temporarily in newest leaf!
				{
					minNreq = nodalUnit[inodeNumber-1].get_leaf()->get_drymass() * NitrogenStatus.minimumleafNitrogenRatio;
					deltaNreq = nodalUnit[inodeNumber-1].get_leaf()->get_drymass() * (NitrogenStatus.optimumleafNitrogenRatio - NitrogenStatus.minimumleafNitrogenRatio);
					luxNreq = (minNreq + deltaNreq) * 0.9;
					nodalUnit[inodeNumber-1].get_leaf()->set_currentNitrogenAmount(minNreq + deltaNreq + luxNreq + tempNleafamount);
					tempNleafamount = 0.;
				}
			}
			//now stems
			for (int i=0; i< inodeNumber; i++)
			{
				if (nodalUnit[i].isInitiated() && !nodalUnit[i].isTerminated()) nodalUnit[i].get_stem()->set_currentNitrogenAmount(tempNstemamount/inodeNumber);
			}
		}


	}
	
}

void CPlant::leaf_output()
{
	int num = 31;//specify ending leaf number for output (leaf 0 is first leaf)
	int start = 0; //starting leaf one is interested in e.g output of leaf 0 to 5, 0 to i< num; output leaf 10 to 15, 11 to i < num, etc.
	if (nodalUnit==NULL) return;
	for (int i = start; i < num+1; i++)
	{
		if (nodalUnit[i].isInitiated()==true)
		{
			lmass[i-start] = nodalUnit[i].get_leaf()->get_drymass();
			page[i-start] = nodalUnit[i].get_leaf()->get_PAge();
			//page[i-start] = nodalUnit[i].get_leaf()->get_SLAleaf();
			cage[i-start] = nodalUnit[i].get_leaf()->get_Cage();
			//cage[i-start]=nodalUnit[i].get_leaf()->get_potGrowth24()/nodalUnit[i].get_leaf()->get_area0();
			area[i-start] = nodalUnit[i].get_leaf()->get_area();
		}
		if (lmass[i-start] < 0) lmass[i-start] = 0.;
		if (page[i-start] < 0) page[i-start] = 0.;
		if (cage[i-start] < 0) cage[i-start] = 0.;
		if (area[i-start] < 0) area[i-start] = 0.;
		if (!nodalUnit[i].isInitiated())
		{
			lmass[i-start] = 0.;
			page[i-start] = 0.;
			cage[i-start] = 0.;
			area[i-start] = 0.;
		}
	}
}

void CPlant::nitrogen_stress_substor(const TInitInfo info)
{
// Routine derived from SUBSTOR, although heavily modified
// Similar to SIMPOTATO except min, critical values for N vary
// Also do not have separate tuber or root buffer - too complicated
	if(nodalUnit==NULL) return;
	//  1) Set optimum, minimum leaf, stem, and tuber amounts
	if (get_develop()->get_xstage() < 2) 
	{
		NitrogenStatus.optimumleafNitrogenRatio = (4.5-0.5*(get_develop()->get_xstage() - 1))/100.;
		NitrogenStatus.minimumleafNitrogenRatio = NitrogenStatus.optimumleafNitrogenRatio - 0.02;
		NitrogenStatus.optimumstemNitrogenRatio = NitrogenStatus.optimumleafNitrogenRatio;
		NitrogenStatus.minimumstemNitrogenRatio = NitrogenStatus.minimumleafNitrogenRatio;
		NitrogenStatus.optimumtuberNitrogenRatio = 0.;
		NitrogenStatus.minimumtuberNitrogenRatio = 0.;
	} else
	{
		NitrogenStatus.optimumleafNitrogenRatio = (4.-1.0*(get_develop()->get_xstage() - 2))/100.;
		NitrogenStatus.minimumleafNitrogenRatio = NitrogenStatus.optimumleafNitrogenRatio - 0.02;
		NitrogenStatus.optimumstemNitrogenRatio = NitrogenStatus.optimumleafNitrogenRatio;
		NitrogenStatus.minimumstemNitrogenRatio = NitrogenStatus.minimumleafNitrogenRatio;
		NitrogenStatus.optimumtuberNitrogenRatio = 0.014;
		NitrogenStatus.minimumtuberNitrogenRatio = 0.014;
	}

	// 2) Set optimum, minimum stem and root amounts
	NitrogenStatus.optimumrootNitrogenRatio = (2.15 - 0.5*get_develop()->get_xstage())/100.;
	if (NitrogenStatus.optimumrootNitrogenRatio < 0.014) NitrogenStatus.optimumrootNitrogenRatio = 0.014;
	if (initInfo.bigleaf == 0 && nodalUnit) //assume each individual nodalunit has same N demand as whole plant index, DHF
	{
		for (int i = 0; i < inodeNumber; i++)
		{
			nodalUnit[i].get_leaf()->set_optimumNitrogenRatio(NitrogenStatus.optimumleafNitrogenRatio);
			nodalUnit[i].get_leaf()->set_minimumNitrogenRatio(NitrogenStatus.minimumleafNitrogenRatio);
			nodalUnit[i].get_stem()->set_optimumNitrogenRatio(NitrogenStatus.optimumstemNitrogenRatio);
			nodalUnit[i].get_stem()->set_minimumNitrogenRatio(NitrogenStatus.minimumstemNitrogenRatio);
		}	
	}
	else
	{
		nodalUnit[0].get_leaf()->set_optimumNitrogenRatio(NitrogenStatus.optimumleafNitrogenRatio);
		nodalUnit[0].get_leaf()->set_minimumNitrogenRatio(NitrogenStatus.minimumleafNitrogenRatio);
		nodalUnit[0].get_stem()->set_optimumNitrogenRatio(NitrogenStatus.optimumstemNitrogenRatio);
		nodalUnit[0].get_stem()->set_minimumNitrogenRatio(NitrogenStatus.minimumstemNitrogenRatio);
	}

// 3) Determine stress factors
	NitrogenStatus.actualleafNitrogenRatio = NitrogenStatus.leafNitrogenAmount / this->get_leafMass();
	if (info.Nitrogen_stress_off==1)
	{
		NitrogenStatus.Nitrogenstressfactor = 1.;
		NitrogenStatus.Nitrogendeficiencyone = 1.;
		NitrogenStatus.Nitrogendeficiencytwo = 1.;
		NitrogenStatus.Nitrogendeficiencythree = NitrogenStatus.Nitrogendeficiencyone;
		this->NitrogenStressFactor = NitrogenStatus.Nitrogenstressfactor;
		return;
	}
	if (get_develop()->get_istage() >= 2) NitrogenStatus.minimumvegetativeNitrogenRatio = NitrogenStatus.minimumleafNitrogenRatio;
	if (get_develop()->get_xstage() <= 1.1) // no N stress yet until several days after emergence
	{
		NitrogenStatus.Nitrogenstressfactor = 1.;
		NitrogenStatus.Nitrogendeficiencyone = 1.;
		NitrogenStatus.Nitrogendeficiencytwo = 1.;
	} else if (NitrogenStatus.actualleafNitrogenRatio >= NitrogenStatus.optimumleafNitrogenRatio)
	{
		NitrogenStatus.Nitrogenstressfactor = 1.;
		NitrogenStatus.Nitrogendeficiencyone = 1.;
		NitrogenStatus.Nitrogendeficiencytwo = 1.;
	} else
	{
		NitrogenStatus.Nitrogenstressfactor = 1 - pow(((NitrogenStatus.optimumleafNitrogenRatio - NitrogenStatus.actualleafNitrogenRatio)/(NitrogenStatus.optimumleafNitrogenRatio - NitrogenStatus.minimumleafNitrogenRatio)),3.);
		if (NitrogenStatus.Nitrogenstressfactor < 0.02) NitrogenStatus.Nitrogenstressfactor = 0.02;
		if (NitrogenStatus.Nitrogenstressfactor > 1.0) NitrogenStatus.Nitrogenstressfactor = 1.0;
		NitrogenStatus.Nitrogendeficiencyone = 1 - ((NitrogenStatus.optimumleafNitrogenRatio - NitrogenStatus.actualleafNitrogenRatio)/(NitrogenStatus.optimumleafNitrogenRatio - NitrogenStatus.minimumleafNitrogenRatio))*0.6;
		NitrogenStatus.Nitrogendeficiencytwo = 1 - ((NitrogenStatus.optimumleafNitrogenRatio - NitrogenStatus.actualleafNitrogenRatio)/(NitrogenStatus.optimumleafNitrogenRatio - NitrogenStatus.minimumleafNitrogenRatio))*0.75;
	}
	if (NitrogenStatus.Nitrogendeficiencyone < 0.1) NitrogenStatus.Nitrogendeficiencyone = 0.1;
	if (NitrogenStatus.Nitrogendeficiencyone > 1) NitrogenStatus.Nitrogendeficiencyone = 1;
	if (NitrogenStatus.Nitrogendeficiencytwo < 0.1) NitrogenStatus.Nitrogendeficiencytwo = 0.1;
	if (NitrogenStatus.Nitrogendeficiencytwo > 1) NitrogenStatus.Nitrogendeficiencytwo = 1;
	NitrogenStatus.Nitrogendeficiencythree = NitrogenStatus.Nitrogendeficiencyone;
	this->NitrogenStressFactor = NitrogenStatus.Nitrogenstressfactor;
}

void CPlant::nitrogen_stress(const TInitInfo info)
{
// ********************************
// Routine derived from SIMPOTATO NFACTO subroutine, N deficiency factor
// This routine sets optimum, minimum N percentages for leaf, stem, root, and tuber tissue based on developmental stage.
// These values decline with xstage.
// Currently assume each nodalUnit leaf and stem has same N opt and min requirements as whole plant index
// N stress factors are calculated based on ratios between differences of optimum, minimum, and actual N concentration in leaf tissue

// Additional Information:
// * N stresses depend on progression through developmental stage, xstage
//   , optimum and minimum leaf, stem, tuber and root contents for N
//   and actual N content as well.
// * N stress variables are used to modify leaf growth, photosynthesis,
//  and soil / root relationships (not used in SPUDSIM)
// * Nitrogenstressfactor runs from 0.02 to 1
// * Nitrogendeficiencyone runs from 0.1 to 1
// * Nitrogendeficiencytwo runs from 0.1 to 1
// * Optimumleafnitrogen runs from 0.08 at emergence to 0.0275 at harvest
// * Minimumleafnitrogen runs from 0.07 at emergence to 0.0175
// * Optimumtubernitrogen runs from 0.027 at TI to 0.015
// * Minimumtubernitrogen runs from 0.023 at TI to 0.011
// * Optimumstemnitrogen runs from 0.05 at emeregence to 0.026
// * Minimumstemnitrogen runs from 0.04 at emergence to 0.016
// * Optimumrootnitrogen runs from 0.024 to 0.014

// ***********************************
	NitrogenStatus.optimumrootNitrogenRatio = __max((2.4 - 0.25*get_develop()->get_xstage())/100.,0.014);
	//NitrogenStatus.optimumrootNitrogenRatio = __max((2.15 - 0.5*get_develop()->get_xstage())/100,0.014);

	if(nodalUnit==NULL) return;
//  1) Set optimum, minimum leaf and tuber amounts
	if (get_develop()->get_xstage() < 2) 
	{
		//NitrogenStatus.optimumleafNitrogenRatio = (8-(get_develop()->get_xstage() - 1))/100;
		//NitrogenStatus.minimumleafNitrogenRatio = NitrogenStatus.optimumleafNitrogenRatio - 0.01;
		NitrogenStatus.optimumleafNitrogenRatio = __max((4.5-0.5*(get_develop()->get_xstage() - 1))/100.,0.0275);
		NitrogenStatus.minimumleafNitrogenRatio = __min(NitrogenStatus.optimumleafNitrogenRatio - 0.02,0.0175);
		NitrogenStatus.optimumtuberNitrogenRatio = 0.;
		NitrogenStatus.minimumtuberNitrogenRatio = 0.;
	} else
	{
		//NitrogenStatus.optimumleafNitrogenRatio = (7-(get_develop()->get_xstage()-2))/100;
		//if (NitrogenStatus.optimumleafNitrogenRatio < 0.0275) NitrogenStatus.optimumleafNitrogenRatio = 0.0275;
		//NitrogenStatus.minimumleafNitrogenRatio = NitrogenStatus.optimumleafNitrogenRatio - 0.01;
		//NitrogenStatus.optimumtuberNitrogenRatio = (2.7 -1/3 * (get_develop()->get_xstage()-2))/100;
		//if (NitrogenStatus.optimumtuberNitrogenRatio < 0.015) NitrogenStatus.optimumtuberNitrogenRatio = 0.015;
		NitrogenStatus.optimumleafNitrogenRatio = __max((4.-1.0*(get_develop()->get_xstage() - 2.))/100.,0.0275);
		NitrogenStatus.minimumleafNitrogenRatio = __min(NitrogenStatus.optimumleafNitrogenRatio - 0.02,0.0175);
		NitrogenStatus.optimumtuberNitrogenRatio = 0.014;
		NitrogenStatus.minimumtuberNitrogenRatio = NitrogenStatus.optimumtuberNitrogenRatio - 0.004;
	}

// 2) Set optimum, minimum stem and root amounts
	//NitrogenStatus.optimumstemNitrogenRatio = (5 - 0.6*(get_develop()->get_xstage()-1))/100;
	//if (NitrogenStatus.optimumstemNitrogenRatio < 0.026) NitrogenStatus.optimumstemNitrogenRatio = 0.026;
	//NitrogenStatus.minimumstemNitrogenRatio = NitrogenStatus.optimumstemNitrogenRatio - 0.01;
	//NitrogenStatus.optimumrootNitrogenRatio = (2.4 - 0.25*get_develop()->get_xstage())/100;
	//if (NitrogenStatus.optimumrootNitrogenRatio < 0.014) NitrogenStatus.optimumrootNitrogenRatio = 0.014;
		NitrogenStatus.optimumstemNitrogenRatio = NitrogenStatus.optimumleafNitrogenRatio;
		NitrogenStatus.minimumstemNitrogenRatio = NitrogenStatus.minimumleafNitrogenRatio;
	//NitrogenStatus.optimumrootNitrogenRatio = __max((2.15 - 0.5*get_develop()->get_xstage())/100,0.014);

	if (initInfo.bigleaf == 0 && nodalUnit) //assume each individual nodalunit has same N demand as whole plant index, DHF
	{
		for (int i = 0; i < inodeNumber; i++)
		{
			nodalUnit[i].get_leaf()->set_optimumNitrogenRatio(NitrogenStatus.optimumleafNitrogenRatio);
			nodalUnit[i].get_leaf()->set_minimumNitrogenRatio(NitrogenStatus.minimumleafNitrogenRatio);
			nodalUnit[i].get_stem()->set_optimumNitrogenRatio(NitrogenStatus.optimumstemNitrogenRatio);
			nodalUnit[i].get_stem()->set_minimumNitrogenRatio(NitrogenStatus.minimumstemNitrogenRatio);
		}	
	}
	else
	{
		nodalUnit[0].get_leaf()->set_optimumNitrogenRatio(NitrogenStatus.optimumleafNitrogenRatio);
		nodalUnit[0].get_leaf()->set_minimumNitrogenRatio(NitrogenStatus.minimumleafNitrogenRatio);
		nodalUnit[0].get_stem()->set_optimumNitrogenRatio(NitrogenStatus.optimumstemNitrogenRatio);
		nodalUnit[0].get_stem()->set_minimumNitrogenRatio(NitrogenStatus.minimumstemNitrogenRatio);
	}
	NitrogenStatus.actualleafNitrogenRatio = NitrogenStatus.leafNitrogenAmount / this->get_leafMass();
	if (info.Nitrogen_stress_off==1)
	{
		NitrogenStatus.Nitrogenstressfactor = 1.;
		NitrogenStatus.Nitrogendeficiencyone = 1.;
		NitrogenStatus.Nitrogendeficiencytwo = 1.;
		NitrogenStatus.Nitrogendeficiencythree = NitrogenStatus.Nitrogendeficiencyone;
		this->NitrogenStressFactor = NitrogenStatus.Nitrogenstressfactor;
		return;
	}
// 3) Determine stress factors
	//base on weighted average of N content in canopy, not just stem or leaf due to N allocation difficulties, DHF
	double Noptimum = (NitrogenStatus.optimumleafNitrogenRatio + NitrogenStatus.optimumstemNitrogenRatio)*0.5;
	double Nminimum = (NitrogenStatus.minimumleafNitrogenRatio + NitrogenStatus.minimumstemNitrogenRatio)*0.5;
	double Nactual = (NitrogenStatus.leafNitrogenAmount+NitrogenStatus.stemNitrogenAmount)/(this->get_leafMass()+this->get_stemMass());
	if (get_develop()->get_istage() >= 2) NitrogenStatus.minimumvegetativeNitrogenRatio = NitrogenStatus.minimumleafNitrogenRatio;
	if (get_develop()->get_xstage() <= 1.1) // no N stress yet until several days after emergence
	{
		NitrogenStatus.Nitrogenstressfactor = 1.;
		NitrogenStatus.Nitrogendeficiencyone = 1.;
		NitrogenStatus.Nitrogendeficiencytwo = 1.;
	} else if (Nactual >= Noptimum)
	{
		NitrogenStatus.Nitrogenstressfactor = 1.;
		NitrogenStatus.Nitrogendeficiencyone = 1.;
		NitrogenStatus.Nitrogendeficiencytwo = 1.;
	} else
	{
		NitrogenStatus.Nitrogenstressfactor = 1 - pow((Noptimum - Nactual)/(Noptimum - Nminimum),3.);
		if (NitrogenStatus.Nitrogenstressfactor < 0.02) NitrogenStatus.Nitrogenstressfactor = 0.02;
		if (NitrogenStatus.Nitrogenstressfactor > 1.0) NitrogenStatus.Nitrogenstressfactor = 1.0;
		NitrogenStatus.Nitrogendeficiencyone = 1 - (Noptimum - Nactual)/(Noptimum - Nminimum)*0.6;
		NitrogenStatus.Nitrogendeficiencytwo = 1 - (Noptimum - Nactual)/(Noptimum - Nminimum)*0.75;

	/*			
		NitrogenStatus.Nitrogenstressfactor = 1 - pow(((NitrogenStatus.optimumleafNitrogenRatio - NitrogenStatus.actualleafNitrogenRatio)/(NitrogenStatus.optimumleafNitrogenRatio - NitrogenStatus.minimumleafNitrogenRatio)),3);
		if (NitrogenStatus.Nitrogenstressfactor < 0.02) NitrogenStatus.Nitrogenstressfactor = 0.02;
		if (NitrogenStatus.Nitrogenstressfactor > 1.0) NitrogenStatus.Nitrogenstressfactor = 1.0;
		NitrogenStatus.Nitrogendeficiencyone = 1 - ((NitrogenStatus.optimumleafNitrogenRatio - NitrogenStatus.actualleafNitrogenRatio)/(NitrogenStatus.optimumleafNitrogenRatio - NitrogenStatus.minimumleafNitrogenRatio))*0.6;
		NitrogenStatus.Nitrogendeficiencytwo = 1 - ((NitrogenStatus.optimumleafNitrogenRatio - NitrogenStatus.actualleafNitrogenRatio)/(NitrogenStatus.optimumleafNitrogenRatio - NitrogenStatus.minimumleafNitrogenRatio))*0.75;
	*/
	}
	if (NitrogenStatus.Nitrogendeficiencyone < 0.1) NitrogenStatus.Nitrogendeficiencyone = 0.1;
	if (NitrogenStatus.Nitrogendeficiencyone > 1) NitrogenStatus.Nitrogendeficiencyone = 1.;
	if (NitrogenStatus.Nitrogendeficiencytwo < 0.1) NitrogenStatus.Nitrogendeficiencytwo = 0.1;
	if (NitrogenStatus.Nitrogendeficiencytwo > 1) NitrogenStatus.Nitrogendeficiencytwo = 1.;
	NitrogenStatus.Nitrogendeficiencythree = NitrogenStatus.Nitrogendeficiencyone;
	this->NitrogenStressFactor = NitrogenStatus.Nitrogenstressfactor;
}

void CPlant::nitrogen_stress_off(const TInitInfo info)
{
/*Routine runs when no nitrogen stress is desired
/*Set all organs to optimum level of N and set availableNitrogen to equal total Nitrogen Demand at each time-step
*/
	NitrogenStatus.actualleafNitrogenRatio = NitrogenStatus.optimumleafNitrogenRatio;
	NitrogenStatus.actualstemNitrogenRatio = NitrogenStatus.optimumstemNitrogenRatio;
	NitrogenStatus.actualrootNitrogenRatio = NitrogenStatus.optimumrootNitrogenRatio;
	NitrogenStatus.leafNitrogendemand = this->dblLeafgro * NitrogenStatus.optimumleafNitrogenRatio;
	NitrogenStatus.stemNitrogendemand = this->dblStemgro * NitrogenStatus.optimumstemNitrogenRatio;
	NitrogenStatus.rootNitrogendemand = this->dblRootgro * NitrogenStatus.optimumrootNitrogenRatio;
	if (this->get_tubers()->isInitiated() == true ) 
	{
		NitrogenStatus.actualtuberNitrogenRatio = NitrogenStatus.optimumtuberNitrogenRatio;
		NitrogenStatus.tuberNitrogendemand = this->dblTubgro * NitrogenStatus.optimumtuberNitrogenRatio;
		NitrogenStatus.totalNitrogendemand = NitrogenStatus.leafNitrogendemand + NitrogenStatus.stemNitrogendemand + NitrogenStatus.rootNitrogendemand + NitrogenStatus.tuberNitrogendemand;
	}
	else NitrogenStatus.totalNitrogendemand = NitrogenStatus.leafNitrogendemand + NitrogenStatus.stemNitrogendemand + NitrogenStatus.rootNitrogendemand;
	if (this->get_tubers()->isInitiated() == true)
	{
		NitrogenStatus.tuberNitrogenAmount = NitrogenStatus.optimumtuberNitrogenRatio * this->get_tuberMass() + NitrogenStatus.tuberNitrogendemand;
	}
	else NitrogenStatus.tuberNitrogenAmount = 0.;
	NitrogenStatus.leafNitrogenAmount = NitrogenStatus.optimumleafNitrogenRatio * this->get_leafMass() + NitrogenStatus.leafNitrogendemand;
	NitrogenStatus.stemNitrogenAmount = NitrogenStatus.optimumstemNitrogenRatio * this->get_stemMass() + NitrogenStatus.stemNitrogendemand;
	NitrogenStatus.rootNitrogenAmount = NitrogenStatus.optimumrootNitrogenRatio * this->get_rootMass() + NitrogenStatus.rootNitrogendemand;
	this->set_TotalPlantNitrogen(NitrogenStatus.leafNitrogenAmount + NitrogenStatus.stemNitrogenAmount + NitrogenStatus.rootNitrogenAmount + NitrogenStatus.tuberNitrogenAmount);
}

void CPlant::nitrogen_balance_pretubers_substor(const TInitInfo info)
{
//***********************************
//* Coefficients come from SUBSTOR - more reasonable, allocation routine similar to SUBSTOR but modified here
//* Main idea is that stem and leaf N demand is identical and both can serve as N buffer

	double leafNdemand, stemNdemand, rootNdemand, surplusN;
	double leafnewgrowthNdemand, stemnewgrowthNdemand, rootnewgrowthNdemand, newgrowthNdemand;
	double leafoldgrowthNdemand, stemoldgrowthNdemand, rootoldgrowthNdemand, oldgrowthNdemand;
	if (this->get_tubers()->isInitiated() == true) return;
	NitrogenStatus.actualleafNitrogenRatio = 0.;
	NitrogenStatus.actualstemNitrogenRatio = 0.;
	NitrogenStatus.actualrootNitrogenRatio = 0.;
	if ((this->get_leafMass()) > 0) NitrogenStatus.actualleafNitrogenRatio = NitrogenStatus.leafNitrogenAmount / (this->get_leafMass()+dblLeafgro);
	if ((this->get_stemMass()) > 0) NitrogenStatus.actualstemNitrogenRatio = NitrogenStatus.stemNitrogenAmount / (this->get_stemMass()+dblStemgro);
	if ((this->get_rootMass()) > 0) NitrogenStatus.actualrootNitrogenRatio = NitrogenStatus.rootNitrogenAmount/ (this->get_rootMass()+dblRootgro);
	NitrogenStatus.leafNitrogendemand = (this->get_leafMass()+this->dblLeafgro)*(NitrogenStatus.optimumleafNitrogenRatio - NitrogenStatus.actualleafNitrogenRatio);
	NitrogenStatus.stemNitrogendemand = (this->get_stemMass()+this->dblStemgro)*(NitrogenStatus.optimumstemNitrogenRatio - NitrogenStatus.actualstemNitrogenRatio);
	NitrogenStatus.rootNitrogendemand = (this->get_rootMass()+this->dblRootgro)*(NitrogenStatus.optimumrootNitrogenRatio - NitrogenStatus.actualrootNitrogenRatio);
	NitrogenStatus.totalNitrogendemand = NitrogenStatus.leafNitrogendemand + NitrogenStatus.stemNitrogendemand + NitrogenStatus.rootNitrogendemand;	
	
	this->set_N_hourlygrowthdemand(__max(0.,NitrogenStatus.totalNitrogendemand));

	NitrogenStatus.leafNitrogenbuffer = __max(((NitrogenStatus.actualleafNitrogenRatio - NitrogenStatus.minimumleafNitrogenRatio)*this->get_leafMass()),0.);
	NitrogenStatus.stemNitrogenbuffer = __max(((NitrogenStatus.actualstemNitrogenRatio - NitrogenStatus.minimumstemNitrogenRatio)*this->get_stemMass()),0.);
	
	leafnewgrowthNdemand = this->dblLeafgro * NitrogenStatus.optimumleafNitrogenRatio;
	stemnewgrowthNdemand = this->dblStemgro * NitrogenStatus.optimumstemNitrogenRatio;
	rootnewgrowthNdemand = this->dblRootgro * NitrogenStatus.optimumrootNitrogenRatio;
	newgrowthNdemand = leafnewgrowthNdemand + stemnewgrowthNdemand + rootnewgrowthNdemand; 	
	
	leafoldgrowthNdemand = NitrogenStatus.leafNitrogendemand - leafnewgrowthNdemand;
	stemoldgrowthNdemand = NitrogenStatus.stemNitrogendemand - stemnewgrowthNdemand;
	rootoldgrowthNdemand = NitrogenStatus.rootNitrogendemand - rootnewgrowthNdemand;
	oldgrowthNdemand = leafoldgrowthNdemand + stemoldgrowthNdemand + rootoldgrowthNdemand;

	// scenarios
	// 0) N in excess, store in haulm
	if (NitrogenStatus.availableNitrogen >= NitrogenStatus.totalNitrogendemand)
	{
		surplusN = (NitrogenStatus.availableNitrogen-NitrogenStatus.totalNitrogendemand);
		leafNdemand = NitrogenStatus.leafNitrogendemand + 0.75*surplusN;
		stemNdemand = NitrogenStatus.stemNitrogendemand + 0.25*surplusN;
		rootNdemand = NitrogenStatus.rootNitrogendemand;
		surplusN = 0.;
	}
	else
	{
	// 1) If N demand not met, add available seedpiece reserve to plant N supply as long as plant is still dependent on seedpiece
		if (this->get_greenLeafArea() < 10000. && N_reserve > 0. && ((NitrogenStatus.actualleafNitrogenRatio-NitrogenStatus.optimumleafNitrogenRatio)* this->get_leafMass() + NitrogenStatus.availableNitrogen < NitrogenStatus.totalNitrogendemand)) // 400 cm2 leaf area same restriction as seedpiece CHO in C_allocation1
		// but N availabe from seed only if N in leaf is below the optimal level
		{
			double diff = NitrogenStatus.totalNitrogendemand - NitrogenStatus.availableNitrogen;
			if (diff >= N_reserve)
			{
				NitrogenStatus.availableNitrogen += N_reserve;
				this->set_N_seedused(N_reserve);
				this->set_N_hourlyseedused(N_reserve);
				N_reserve = 0.;
			}
			else 
			{
				NitrogenStatus.availableNitrogen += diff;
				N_reserve -= diff;
				this->set_N_seedused(diff);
				this->set_N_hourlyseedused(diff);
			}
		}
	// 1a) Demand can now be satisfied
		if (NitrogenStatus.availableNitrogen >= NitrogenStatus.totalNitrogendemand)
		{
			leafNdemand = NitrogenStatus.leafNitrogendemand;
			stemNdemand = NitrogenStatus.stemNitrogendemand;
			rootNdemand = NitrogenStatus.rootNitrogendemand;
		}
	// 2) Still not enough N from seedpiece to satisfy demand, determine if enough N is available to satisfy just new growth N
		else
		{
	// 2a) Enough to satisfy new growth N demand and some older growth, allocated to meet new growth first, then replenish leaf and stem buffers
			if (NitrogenStatus.availableNitrogen >= newgrowthNdemand)
			{
				surplusN = NitrogenStatus.availableNitrogen - newgrowthNdemand;
				leafNdemand = leafnewgrowthNdemand;
				stemNdemand = stemnewgrowthNdemand;
				rootNdemand = rootnewgrowthNdemand;
				if (surplusN > 0) //allocate remaining N to either restore leaf N content to optimum, and if still in excess, do the same for stem
				{
					if (leafoldgrowthNdemand >= surplusN) 
					{
						leafNdemand += surplusN;
						surplusN = 0.;
					}
					else
					{
						leafNdemand += leafoldgrowthNdemand;
						surplusN -= leafoldgrowthNdemand;
						stemNdemand += surplusN;
						surplusN = 0.;
					}
				}
			}
	// 2b) New growth N demand can't be satisfied.  Will need to utilize N from leaf and stem buffer to satisfy new growth, limited to 5% each hour
			else
			{
				/*
				surplusN = __min(0.05*(NitrogenStatus.leafNitrogenbuffer  + NitrogenStatus.stemNitrogenbuffer),(newgrowthNdemand - NitrogenStatus.availableNitrogen));
				NitrogenStatus.stemNitrogenAmount -= 0.75*surplusN;
				NitrogenStatus.leafNitrogenAmount -= 0.25*surplusN;
				*/
				surplusN = NitrogenStatus.leafNitrogenbuffer + NitrogenStatus.stemNitrogenbuffer;
				NitrogenStatus.stemNitrogenAmount -= NitrogenStatus.stemNitrogenbuffer;
				NitrogenStatus.leafNitrogenAmount -= NitrogenStatus.leafNitrogenbuffer;
				NitrogenStatus.availableNitrogen += surplusN;
				this->set_N_hourlytranslocated(surplusN);
	// 2c) Buffer can used to meet N demand
				if ((NitrogenStatus.availableNitrogen) >= (newgrowthNdemand))
				{
					leafNdemand = leafnewgrowthNdemand;
					stemNdemand = stemnewgrowthNdemand;
					rootNdemand = rootnewgrowthNdemand;
					surplusN -= (leafnewgrowthNdemand + stemnewgrowthNdemand + rootnewgrowthNdemand);				}
	// 2d) Still not enough buffer to meet N demand for new growth, reduce new growth accordingly
				else
				{
					double Nfactor = (NitrogenStatus.availableNitrogen+surplusN)/newgrowthNdemand;
					dblLeafgro *= Nfactor;
					dblStemgro *= Nfactor;
					dblRootgro *= Nfactor;
					leafNdemand = this->dblLeafgro * NitrogenStatus.optimumleafNitrogenRatio;
					stemNdemand = this->dblStemgro * NitrogenStatus.optimumstemNitrogenRatio;
					rootNdemand = this->dblRootgro * NitrogenStatus.optimumrootNitrogenRatio;
					surplusN = 0.;
				}
			}
		}
	}
	NitrogenStatus.leafNitrogendemand = leafNdemand;
	NitrogenStatus.stemNitrogendemand = stemNdemand;
	NitrogenStatus.rootNitrogendemand = rootNdemand;
	NitrogenStatus.totalNitrogendemand = leafNdemand + stemNdemand + rootNdemand;
	NitrogenStatus.availableNitrogen -= (leafNdemand + stemNdemand + rootNdemand);
 	NitrogenStatus.leafNitrogenAmount += NitrogenStatus.leafNitrogendemand;
	NitrogenStatus.stemNitrogenAmount += NitrogenStatus.stemNitrogendemand;
	NitrogenStatus.rootNitrogenAmount += NitrogenStatus.rootNitrogendemand;
	// Recalculate N content of leaves, stems, roots, now that new growth has conceptually been added to existing
	if (this->get_leafMass()+dblLeafgro > 0) NitrogenStatus.actualleafNitrogenRatio = NitrogenStatus.leafNitrogenAmount / (this->get_leafMass()+dblLeafgro);
	if (this->get_stemMass()+dblStemgro > 0) NitrogenStatus.actualstemNitrogenRatio = NitrogenStatus.stemNitrogenAmount / (this->get_stemMass()+dblStemgro);
	if (this->get_rootMass()+dblRootgro > 0) NitrogenStatus.actualrootNitrogenRatio = NitrogenStatus.rootNitrogenAmount / (this->get_rootMass()+dblRootgro);
	if (NitrogenStatus.availableNitrogen < 0) NitrogenStatus.availableNitrogen = 0.;
	if (NitrogenStatus.availableNitrogen > 0)
	{
		if (this->get_leafMass() > 0 && this->get_stemMass() > 0)
		{
			if (NitrogenStatus.availableNitrogen >= (leafoldgrowthNdemand + stemoldgrowthNdemand))
			{
				NitrogenStatus.leafNitrogenAmount += leafoldgrowthNdemand;
				NitrogenStatus.stemNitrogenAmount += stemoldgrowthNdemand;
				NitrogenStatus.availableNitrogen -= (leafoldgrowthNdemand + stemoldgrowthNdemand);
				if (NitrogenStatus.availableNitrogen > 0.)
				{
					NitrogenStatus.leafNitrogenAmount += 0.75*NitrogenStatus.availableNitrogen;
					NitrogenStatus.stemNitrogenAmount += 0.25*NitrogenStatus.availableNitrogen;
					NitrogenStatus.availableNitrogen = 0.;
				}
			}
			else
			{
				NitrogenStatus.stemNitrogenAmount += NitrogenStatus.availableNitrogen*0.25;
				NitrogenStatus.leafNitrogenAmount += NitrogenStatus.availableNitrogen*0.75;
				NitrogenStatus.availableNitrogen = 0.;
			}

		}

		NitrogenStatus.actualleafNitrogenRatio = NitrogenStatus.leafNitrogenAmount / (this->get_leafMass()+dblLeafgro);
		NitrogenStatus.actualstemNitrogenRatio = NitrogenStatus.stemNitrogenAmount / (this->get_stemMass()+dblStemgro);
	}
	this->set_TotalPlantNitrogen(NitrogenStatus.leafNitrogenAmount+NitrogenStatus.stemNitrogenAmount+NitrogenStatus.rootNitrogenAmount);


}
void CPlant::nitrogen_balance_pretubers_simplified(const TInitInfo info)
{
// *********************************************
// Routine derived from SIMPOTATO PTTNVG and PARTTN subroutines, but heavily modified
// Major change is to base N demand on new growth only
// C allocation is not directly affected by N demand except through reduction in photosynthetic rate based on pooled leaf N content sitting below minimum
//  threshold value
	double ratio = 0., temp;
	temp = 0.;
	if (this->get_tubers()->isInitiated() == true) return;
	NitrogenStatus.actualleafNitrogenRatio = 0.;
	NitrogenStatus.actualstemNitrogenRatio = 0.;
	NitrogenStatus.actualrootNitrogenRatio = 0.;
	if ((this->get_leafMass()) > 0) NitrogenStatus.actualleafNitrogenRatio = NitrogenStatus.leafNitrogenAmount / (this->get_leafMass());
	if ((this->get_stemMass()) > 0) NitrogenStatus.actualstemNitrogenRatio = NitrogenStatus.stemNitrogenAmount / (this->get_stemMass());
	if ((this->get_rootMass()) > 0) NitrogenStatus.actualrootNitrogenRatio = NitrogenStatus.rootNitrogenAmount/ (this->get_rootMass());
	NitrogenStatus.leafNitrogendemand = (this->dblLeafgro*NitrogenStatus.optimumleafNitrogenRatio);
	NitrogenStatus.stemNitrogendemand = (this->dblStemgro*NitrogenStatus.optimumstemNitrogenRatio);
	NitrogenStatus.rootNitrogendemand = (this->dblRootgro*NitrogenStatus.optimumrootNitrogenRatio);
	NitrogenStatus.totalNitrogendemand = NitrogenStatus.leafNitrogendemand + NitrogenStatus.stemNitrogendemand + NitrogenStatus.rootNitrogendemand;
	// scenarios
	// 0) If N demand not met, add available seedpiece reserve to plant N supply as long as plant is still dependent on seedpiece
	if (NitrogenStatus.availableNitrogen < NitrogenStatus.totalNitrogendemand)
	{
		if (this->get_greenLeafArea() < 1000. && N_reserve > 0)
		// 400 cm2 leaf area same restriction as seedpiece CHO in C_allocation1
		// but N availabe from seed only if N in leaf is below the optimal level
		{
			double diff = NitrogenStatus.totalNitrogendemand - NitrogenStatus.availableNitrogen;
			if (diff >= N_reserve)
			{
				NitrogenStatus.availableNitrogen += N_reserve;
				this->set_N_seedused(N_reserve);
				N_reserve = 0.;
			}
			else 
			{
				NitrogenStatus.availableNitrogen += diff;
				N_reserve -= diff;
				this->set_N_seedused(diff);
			}
		}
	}
	// 1) Luxury consumption - more N available than needed, increase leaf demand to store surplus
	if (NitrogenStatus.availableNitrogen > NitrogenStatus.totalNitrogendemand)
	{
		if (NitrogenStatus.totalNitrogendemand > 0) // compute surplus after removing demand for new growth
		{
			temp = NitrogenStatus.availableNitrogen - NitrogenStatus.totalNitrogendemand;
		}
		else
		{
			if (NitrogenStatus.leafNitrogendemand > 0) temp -= NitrogenStatus.leafNitrogendemand;
			if (NitrogenStatus.stemNitrogendemand > 0) temp -= NitrogenStatus.stemNitrogendemand;
			if (NitrogenStatus.rootNitrogendemand > 0) temp -= NitrogenStatus.rootNitrogendemand;
			temp += NitrogenStatus.availableNitrogen;
		}
		NitrogenStatus.leafNitrogendemand += temp;
		NitrogenStatus.totalNitrogendemand = NitrogenStatus.leafNitrogendemand + NitrogenStatus.rootNitrogendemand + NitrogenStatus.stemNitrogendemand;
	}
	// 2) Not enough N - check for surplus in leaves and redistribute
	else
	{
		if (NitrogenStatus.actualleafNitrogenRatio > NitrogenStatus.optimumleafNitrogenRatio) // surplus condition in existing leaf biomass
		{
			NitrogenStatus.leafNitrogenbuffer = (this->get_leafMass())*(NitrogenStatus.actualleafNitrogenRatio-NitrogenStatus.optimumleafNitrogenRatio);
			NitrogenStatus.availableNitrogen += NitrogenStatus.leafNitrogenbuffer;
			NitrogenStatus.leafNitrogenAmount -= NitrogenStatus.leafNitrogenbuffer;
			NitrogenStatus.leafNitrogenbuffer = 0.;
			// 2a) Leaf Surplus is enough to meet demand, move leaf N into stems and roots to meet demand
			if (NitrogenStatus.totalNitrogendemand <= (NitrogenStatus.availableNitrogen + NitrogenStatus.leafNitrogenbuffer))
			{
				//NitrogenStatus.stemNitrogenAmount += NitrogenStatus.stemNitrogendemand;
				//NitrogenStatus.rootNitrogenAmount += NitrogenStatus.rootNitrogendemand;
				//NitrogenStatus.leafNitrogenAmount += NitrogenStatus.leafNitrogendemand;
				//NitrogenStatus.availableNitrogen -= (NitrogenStatus.stemNitrogendemand + NitrogenStatus.rootNitrogendemand + NitrogenStatus.leafNitrogendemand);
			}
			// 2b) Leaf surplus not enough, so take N away from minimum N contents in leaf and stem to meet demand
			else
			{
				if (this->get_leafMass() > 0) NitrogenStatus.actualleafNitrogenRatio = NitrogenStatus.leafNitrogenAmount / this->get_leafMass();
				if (this->get_stemMass() > 0) NitrogenStatus.actualstemNitrogenRatio = NitrogenStatus.stemNitrogenAmount / this->get_stemMass();
				NitrogenStatus.leafNitrogenbuffer = max(this->get_leafMass()*(NitrogenStatus.actualleafNitrogenRatio - NitrogenStatus.minimumleafNitrogenRatio),0.);
				NitrogenStatus.leafNitrogenAmount -= NitrogenStatus.leafNitrogenbuffer;
				NitrogenStatus.stemNitrogenbuffer = max(this->get_stemMass()*(NitrogenStatus.actualstemNitrogenRatio - NitrogenStatus.minimumstemNitrogenRatio),0.);
				NitrogenStatus.stemNitrogenAmount -= NitrogenStatus.stemNitrogenbuffer;
				NitrogenStatus.availableNitrogen += NitrogenStatus.leafNitrogenbuffer + NitrogenStatus.stemNitrogenbuffer;
				NitrogenStatus.leafNitrogenbuffer = 0.;
				NitrogenStatus.stemNitrogenbuffer = 0.;
				// recheck demand now that leaf and stem buffer were added
				// N more than enough, increase leaf and stem demand:
				if (NitrogenStatus.totalNitrogendemand < NitrogenStatus.availableNitrogen) 
				{
					if ((NitrogenStatus.leafNitrogendemand > 0)||(NitrogenStatus.stemNitrogendemand > 0))
					{
						//ratio will be greater than 1 here, indicating that N will provide new mass at optimum N level plus will return some N to older mass
						ratio = (NitrogenStatus.availableNitrogen - NitrogenStatus.rootNitrogendemand) / (NitrogenStatus.leafNitrogendemand+NitrogenStatus.stemNitrogendemand);
						NitrogenStatus.leafNitrogendemand *= ratio;
						NitrogenStatus.stemNitrogendemand *= ratio;
					}
					else
					{
						NitrogenStatus.leafNitrogendemand = (NitrogenStatus.availableNitrogen - NitrogenStatus.totalNitrogendemand)/2.;
						NitrogenStatus.stemNitrogendemand = (NitrogenStatus.availableNitrogen - NitrogenStatus.totalNitrogendemand)/2.;
					}
				}
				//If N not enough, must reduce growth of new organs
				else
				{
					ratio = NitrogenStatus.availableNitrogen / NitrogenStatus.totalNitrogendemand;
				
					dblLeafgro *= ratio;
					dblStemgro *= ratio;
					dblRootgro *= ratio;
				
					NitrogenStatus.leafNitrogendemand *= ratio;
					NitrogenStatus.stemNitrogendemand *= ratio;
					NitrogenStatus.rootNitrogendemand *= ratio;
				}
			}
		}
	}
	NitrogenStatus.totalNitrogendemand = NitrogenStatus.leafNitrogendemand + NitrogenStatus.stemNitrogendemand + NitrogenStatus.rootNitrogendemand;
	// Check for imbalances in N demand versus available N and correct
	if (NitrogenStatus.totalNitrogendemand > 1.01*NitrogenStatus.availableNitrogen)
	{
		ratio = NitrogenStatus.availableNitrogen/NitrogenStatus.totalNitrogendemand;
		NitrogenStatus.leafNitrogendemand *= ratio;
		NitrogenStatus.stemNitrogendemand *= ratio;
		NitrogenStatus.rootNitrogendemand *= ratio;
	}
	// Add available N to plant parts
	//note, this section also makes sure that any residual stem and leaf N demand not satisfied from leaf N buffer is satisfied with available N
 	NitrogenStatus.leafNitrogenAmount += NitrogenStatus.leafNitrogendemand;
	NitrogenStatus.stemNitrogenAmount += NitrogenStatus.stemNitrogendemand;
	NitrogenStatus.rootNitrogenAmount += NitrogenStatus.rootNitrogendemand;
	// Recalculate N content of leaves, stems, roots, now that new growth has conceptually been added to existing
	if (this->get_leafMass()+dblLeafgro > 0) NitrogenStatus.actualleafNitrogenRatio = NitrogenStatus.leafNitrogenAmount / (this->get_leafMass()+dblLeafgro);
	if (this->get_stemMass()+dblStemgro > 0) NitrogenStatus.actualstemNitrogenRatio = NitrogenStatus.stemNitrogenAmount / (this->get_stemMass()+dblStemgro);
	if (this->get_rootMass()+dblRootgro > 0) NitrogenStatus.actualrootNitrogenRatio = NitrogenStatus.rootNitrogenAmount / (this->get_rootMass()+dblRootgro);
	NitrogenStatus.availableNitrogen -=(NitrogenStatus.leafNitrogendemand+NitrogenStatus.stemNitrogendemand+NitrogenStatus.rootNitrogendemand);
	if (NitrogenStatus.availableNitrogen < 0) NitrogenStatus.availableNitrogen = 0;
	if (this->get_leafMass() > 0) //DHF - Add surplus N into leaf pool or stem pool
	{
		NitrogenStatus.leafNitrogenAmount += NitrogenStatus.availableNitrogen; //DHF - Adding available N into leaf class as surplus
		NitrogenStatus.availableNitrogen = 0.;
		NitrogenStatus.actualleafNitrogenRatio = NitrogenStatus.leafNitrogenAmount / (this->get_leafMass()+dblLeafgro);
	}
	else if (this->get_stemMass() > 0)
	{
		NitrogenStatus.stemNitrogenAmount += NitrogenStatus.availableNitrogen;
		NitrogenStatus.availableNitrogen = 0.;
		NitrogenStatus.actualstemNitrogenRatio = NitrogenStatus.stemNitrogenAmount / (this->get_stemMass()+dblStemgro);
	}
	this->set_TotalPlantNitrogen(NitrogenStatus.leafNitrogenAmount+NitrogenStatus.stemNitrogenAmount+NitrogenStatus.rootNitrogenAmount+NitrogenStatus.tuberNitrogenAmount);
}	
			
void CPlant::nitrogen_balance_pretubers(const TInitInfo info)
{
// ********************************
// Routine derived from SIMPOTATO PRTTNVG and PARTTN subroutines
// Difference between subroutines is based on developmental stage of the plant (pre vs post TI)
// Routine adjusts actual C growth based on N deficiency - therefore, must run after C_Allocation1 method.
// Simulates whole N demand by haulm and roots based on potential growth
// act__Nitrogen variables uses units of g N g organ-1
// ____Nitrogen and ___Nitrogendemand variables use units of g N plant-1
// NOTE: major difference from SIMPOTATO is that Navailable is NOT potential, but actual growth
// NOTE: although logic appears convoluted at first, trace-through indicates algorithm works as intended - I added additional comments

// First calculate N demand for existing biomass and potential growth
// Compare N demand to total N uptake
// If N is limiting, take the following steps in order to reduce N
//    demand until availn meets demand:
// (1) Use N from seed reserve if any is left (this is done in carbon_allocation 1)
// (2) If N content is above the critical N level (i.e. prior luxury consumption of N) then redistribute N.
// (3) Reduce growth to fit within available N
// (4) N content for root growth cannot be remobilized, only stems and leaves
// Right now, subroutine is independent of bigleaf or individual leaf model as N demand based on total leaf or stem mass
// ********************************

	double ratio = 0., temp;
	temp = 0.;
	if (this->get_tubers()->isInitiated() == true) return; //only run routine if pre-TI
	
	// Ndemand equals difference between Nopt - Nactual for all (existing + new) growth of each organ class
	//	Recall that Noptimal is the N content needed to avoid any stress

	NitrogenStatus.actualleafNitrogenRatio = 0.;
	NitrogenStatus.actualstemNitrogenRatio = 0.;
	NitrogenStatus.actualrootNitrogenRatio = 0.;
	if ((this->get_leafMass() + this->dblLeafgro) > 0) NitrogenStatus.actualleafNitrogenRatio = NitrogenStatus.leafNitrogenAmount / (this->get_leafMass()+this->dblLeafgro);
	if ((this->get_stemMass() + this->dblStemgro) > 0) NitrogenStatus.actualstemNitrogenRatio = NitrogenStatus.stemNitrogenAmount / (this->get_stemMass()+this->dblStemgro);
	if ((this->get_rootMass() + this->dblRootgro) > 0) NitrogenStatus.actualrootNitrogenRatio = NitrogenStatus.rootNitrogenAmount/ (this->get_rootMass()+this->dblRootgro);
	NitrogenStatus.leafNitrogendemand = (this->get_leafMass()+this->dblLeafgro)*(NitrogenStatus.optimumleafNitrogenRatio - NitrogenStatus.actualleafNitrogenRatio);
	NitrogenStatus.stemNitrogendemand = (this->get_stemMass()+this->dblStemgro)*(NitrogenStatus.optimumstemNitrogenRatio - NitrogenStatus.actualstemNitrogenRatio);
	NitrogenStatus.rootNitrogendemand = (this->get_rootMass()+this->dblRootgro)*(NitrogenStatus.optimumrootNitrogenRatio - NitrogenStatus.actualrootNitrogenRatio);
	NitrogenStatus.totalNitrogendemand = NitrogenStatus.leafNitrogendemand + NitrogenStatus.stemNitrogendemand + NitrogenStatus.rootNitrogendemand;
 	
	// Scenarios
	// 0) First, if N demand not met, simply add any available seedpiece reserve to plant N supply as long as plant is still dependent on seedpiece = DHF
	if (NitrogenStatus.availableNitrogen < NitrogenStatus.totalNitrogendemand)
	{
		if (this->get_greenLeafArea() < 10000 && N_reserve > 0 && (NitrogenStatus.optimumleafNitrogenRatio > NitrogenStatus.actualleafNitrogenRatio)) 
		// 400 cm2 leaf area same restriction as seedpiece CHO in C_allocation1
		// but N availabe from seed only if N in leaf is below the optimal level
		{
			double diff = NitrogenStatus.totalNitrogendemand - NitrogenStatus.availableNitrogen;
			if (diff >= N_reserve)
			{
				NitrogenStatus.availableNitrogen += N_reserve;
				this->set_N_seedused(N_reserve);
				N_reserve = 0.;
			}
			else 
			{
				NitrogenStatus.availableNitrogen += diff;
				N_reserve -= diff;
				this->set_N_seedused(diff);
			}
		}
	}
	// 1) Luxury consumption - more N available than needed, increase leaf demand by up to 90% of current value
	if (NitrogenStatus.availableNitrogen >= NitrogenStatus.totalNitrogendemand)
	{
		if (NitrogenStatus.totalNitrogendemand > 0) { //following if-then needed to make sure only current surplus N is added to temp variable
			temp = NitrogenStatus.availableNitrogen - NitrogenStatus.totalNitrogendemand; 
		}
		else
		{
			if (NitrogenStatus.leafNitrogendemand > 0) temp -= NitrogenStatus.leafNitrogendemand;
			if (NitrogenStatus.stemNitrogendemand > 0) temp -= NitrogenStatus.stemNitrogendemand;
			if (NitrogenStatus.rootNitrogendemand > 0) temp -= NitrogenStatus.rootNitrogendemand;
			temp += NitrogenStatus.availableNitrogen;
		}
		//if (temp > NitrogenStatus.leafNitrogendemand*0.9) temp = NitrogenStatus.leafNitrogendemand*0.9; 
		NitrogenStatus.leafNitrogendemand += temp; 
		NitrogenStatus.totalNitrogendemand = NitrogenStatus.leafNitrogendemand + NitrogenStatus.rootNitrogendemand + NitrogenStatus.stemNitrogendemand;
	}
	// 2) Not enough N - check for surplus in leaves and redistribute
	else 
	{
		if (NitrogenStatus.actualleafNitrogenRatio > NitrogenStatus.optimumleafNitrogenRatio) //surplus condition
		{
			NitrogenStatus.leafNitrogenbuffer = (this->get_leafMass()+this->dblLeafgro)*(NitrogenStatus.actualleafNitrogenRatio-NitrogenStatus.optimumleafNitrogenRatio);
	// 2a) Leaf Surplus is enough to meet demand, move leaf N into stems and roots to meet demand
			if (NitrogenStatus.totalNitrogendemand <= (NitrogenStatus.availableNitrogen + NitrogenStatus.leafNitrogenbuffer))
			{
				//use up leaf N buffer, giving stems priority over roots
				//at end of this section, stem and/or root N demand still may not be satisfied and leaf N buffer may or may not be zero
				//however, availN still not used yet
				if (NitrogenStatus.stemNitrogendemand > 0) //remobilize N as needed
				{
					if (NitrogenStatus.stemNitrogendemand > NitrogenStatus.leafNitrogenbuffer)
					{
						NitrogenStatus.stemNitrogenAmount += NitrogenStatus.leafNitrogenbuffer;
						NitrogenStatus.stemNitrogendemand -= NitrogenStatus.leafNitrogenbuffer;
						NitrogenStatus.leafNitrogenAmount -= NitrogenStatus.leafNitrogenbuffer;
						NitrogenStatus.leafNitrogenbuffer = 0.0;
					}
					else
					{
						NitrogenStatus.stemNitrogenAmount += NitrogenStatus.stemNitrogendemand;
						NitrogenStatus.leafNitrogenAmount -= NitrogenStatus.stemNitrogendemand;
						NitrogenStatus.leafNitrogenbuffer -= NitrogenStatus.stemNitrogendemand;
						NitrogenStatus.stemNitrogendemand = 0.0;
					}
				}
				if ((NitrogenStatus.rootNitrogendemand > 0) && (NitrogenStatus.leafNitrogenbuffer > 0))
				{
					if (NitrogenStatus.rootNitrogendemand > NitrogenStatus.leafNitrogenbuffer)
					{
						NitrogenStatus.rootNitrogenAmount += NitrogenStatus.leafNitrogenbuffer;
						NitrogenStatus.rootNitrogendemand -= NitrogenStatus.leafNitrogenbuffer;
						NitrogenStatus.leafNitrogenAmount -= NitrogenStatus.leafNitrogenbuffer;
						NitrogenStatus.leafNitrogenbuffer = 0.;
					}
					else
					{
						NitrogenStatus.rootNitrogenAmount += NitrogenStatus.rootNitrogendemand;
						NitrogenStatus.leafNitrogenAmount -= NitrogenStatus.rootNitrogendemand;
						NitrogenStatus.leafNitrogenbuffer -= NitrogenStatus.rootNitrogendemand;
						NitrogenStatus.rootNitrogendemand = 0.;
					}
				}
			}
	// 2b) Surplus not enough to meet N demand, so give new growth priority for N and recalculate demands
			// check leaf and stem new growth N requirements first:
			else
			{
				NitrogenStatus.leafNitrogenAmount -= NitrogenStatus.leafNitrogenbuffer;
				NitrogenStatus.availableNitrogen += NitrogenStatus.leafNitrogenbuffer;
				NitrogenStatus.leafNitrogendemand = this->dblLeafgro*NitrogenStatus.optimumleafNitrogenRatio;
				NitrogenStatus.stemNitrogendemand = this->dblStemgro*NitrogenStatus.optimumstemNitrogenRatio;
				NitrogenStatus.totalNitrogendemand = NitrogenStatus.leafNitrogendemand + NitrogenStatus.stemNitrogendemand + NitrogenStatus.rootNitrogendemand;
				if (NitrogenStatus.availableNitrogen > NitrogenStatus.totalNitrogendemand) //now there's enough N, determine additional fraction to move into older biomass
				{
					//in this case, ratio will exceed a value of 1, so leaf and stem N will go towards both new growth and older dry mass
					ratio = (NitrogenStatus.availableNitrogen - NitrogenStatus.rootNitrogendemand) / (NitrogenStatus.leafNitrogendemand + NitrogenStatus.stemNitrogendemand);
					NitrogenStatus.leafNitrogendemand = NitrogenStatus.leafNitrogendemand * ratio;
					NitrogenStatus.stemNitrogendemand = NitrogenStatus.stemNitrogendemand * ratio;
				}
				else // N still not enough, so restrict N to all organ classes in equal proportions
				{
					//in this case, ratio will be less than 1, so all organ N demand will be decreased in same ratio
					ratio = NitrogenStatus.availableNitrogen / NitrogenStatus.totalNitrogendemand;
 					
					dblLeafgro *= ratio;
					dblStemgro *= ratio;
					dblRootgro *= ratio;
					
					NitrogenStatus.leafNitrogendemand *= ratio;
					NitrogenStatus.stemNitrogendemand *= ratio;
					NitrogenStatus.rootNitrogendemand *= ratio;
				}
			}
			NitrogenStatus.totalNitrogendemand = NitrogenStatus.leafNitrogendemand + NitrogenStatus.stemNitrogendemand + NitrogenStatus.rootNitrogendemand;
			ratio = 0.;
		}
	// 2c) Still not enough N after using up leaf N buffer, take N in leaves and stems above minimum concentration and add to N available for new growth
		if (NitrogenStatus.totalNitrogendemand > NitrogenStatus.availableNitrogen)
		{
			// 2c-i) find leaf N above the minimum level
			NitrogenStatus.actualleafNitrogenRatio = 0.;
			if (this->get_leafMass() > 0) NitrogenStatus.actualleafNitrogenRatio = NitrogenStatus.leafNitrogenAmount / this->get_leafMass();
			NitrogenStatus.leafNitrogenbuffer = max(this->get_leafMass()*(NitrogenStatus.actualleafNitrogenRatio - NitrogenStatus.minimumleafNitrogenRatio),0.);
			if (NitrogenStatus.leafNitrogenbuffer < NitrogenStatus.leafNitrogenAmount){
				NitrogenStatus.leafNitrogenAmount -= NitrogenStatus.leafNitrogenbuffer;
			}
			/*else {
				//indicates that leaf N is already below minimum possible level
				//original SIMPOTATO code uses error flag here to indicate this to user
			}
			*/
			// 2c-ii) find stem N avove the minimum level
			NitrogenStatus.actualstemNitrogenRatio = 0.;
			if (this->get_stemMass() > 0) NitrogenStatus.actualstemNitrogenRatio = NitrogenStatus.stemNitrogenAmount / this->get_stemMass();
			NitrogenStatus.stemNitrogenbuffer = this->get_stemMass() * (NitrogenStatus.actualstemNitrogenRatio-NitrogenStatus.minimumstemNitrogenRatio);
			if (NitrogenStatus.stemNitrogenbuffer < 0) NitrogenStatus.stemNitrogenbuffer = 0.;
			if (NitrogenStatus.stemNitrogenbuffer < NitrogenStatus.stemNitrogenAmount) {
				NitrogenStatus.stemNitrogenAmount -= NitrogenStatus.stemNitrogenbuffer;
			}
			/* else {
				//indicates that stem N is already below minimum possible level
				//original SIMPOTATO code uses error flag here to indicate this to user
				}
			*/
			NitrogenStatus.availableNitrogen = NitrogenStatus.availableNitrogen + NitrogenStatus.leafNitrogenbuffer + NitrogenStatus.stemNitrogenbuffer;
			NitrogenStatus.leafNitrogendemand = this->dblLeafgro*NitrogenStatus.optimumleafNitrogenRatio;
			NitrogenStatus.stemNitrogendemand = this->dblStemgro*NitrogenStatus.optimumstemNitrogenRatio;
			NitrogenStatus.totalNitrogendemand = NitrogenStatus.leafNitrogendemand + NitrogenStatus.stemNitrogendemand + NitrogenStatus.rootNitrogendemand;
			// recheck demand now that leaf and stem buffer were added
			// N more than enough, increase leaf and stem demand:
			if (NitrogenStatus.totalNitrogendemand < NitrogenStatus.availableNitrogen) 
			{
				if ((NitrogenStatus.leafNitrogendemand > 0)||(NitrogenStatus.stemNitrogendemand > 0))
				{
					//ratio will be greater than 1 here, indicating that N will provide new mass at optimum N level plus will return some N to older mass
					ratio = (NitrogenStatus.availableNitrogen - NitrogenStatus.rootNitrogendemand) / (NitrogenStatus.leafNitrogendemand+NitrogenStatus.stemNitrogendemand);
					NitrogenStatus.leafNitrogendemand *= ratio;
					NitrogenStatus.stemNitrogendemand *= ratio;
				}
				else
				{
					NitrogenStatus.leafNitrogendemand = (NitrogenStatus.availableNitrogen - NitrogenStatus.totalNitrogendemand)/2.;
					NitrogenStatus.stemNitrogendemand = (NitrogenStatus.availableNitrogen - NitrogenStatus.totalNitrogendemand)/2.;
				}
			}
			//If N not enough, reduce growth
			else
			{
				ratio = NitrogenStatus.availableNitrogen / NitrogenStatus.totalNitrogendemand;
				
				dblLeafgro *= ratio;
				dblStemgro *= ratio;
				dblRootgro *= ratio;
				
				NitrogenStatus.leafNitrogendemand *= ratio;
				NitrogenStatus.stemNitrogendemand *= ratio;
				NitrogenStatus.rootNitrogendemand *= ratio;
			}
		}
	}
	NitrogenStatus.totalNitrogendemand = NitrogenStatus.leafNitrogendemand + NitrogenStatus.stemNitrogendemand + NitrogenStatus.rootNitrogendemand;
	// Check for imbalances in N demand versus available N and correct
	if (NitrogenStatus.totalNitrogendemand > 1.01*NitrogenStatus.availableNitrogen)
	{
		ratio = NitrogenStatus.availableNitrogen/NitrogenStatus.totalNitrogendemand;
		NitrogenStatus.leafNitrogendemand *= ratio;
		NitrogenStatus.stemNitrogendemand *= ratio;
		NitrogenStatus.rootNitrogendemand *= ratio;
	}
	// Add available N to plant parts
	//note, this section also makes sure that any residual stem and leaf N demand not satisfied from leaf N buffer is satisfied with available N
 	NitrogenStatus.leafNitrogenAmount += NitrogenStatus.leafNitrogendemand;
	NitrogenStatus.stemNitrogenAmount += NitrogenStatus.stemNitrogendemand;
	NitrogenStatus.rootNitrogenAmount += NitrogenStatus.rootNitrogendemand;
	// Recalculate N content of leaves, stems, roots
	if (this->get_leafMass()+dblLeafgro > 0) NitrogenStatus.actualleafNitrogenRatio = NitrogenStatus.leafNitrogenAmount / (this->get_leafMass()+dblLeafgro);
	if (this->get_stemMass()+dblStemgro > 0) NitrogenStatus.actualstemNitrogenRatio = NitrogenStatus.stemNitrogenAmount / (this->get_stemMass()+dblStemgro);
	if (this->get_rootMass()+dblRootgro > 0) NitrogenStatus.actualrootNitrogenRatio = NitrogenStatus.rootNitrogenAmount / (this->get_rootMass()+dblRootgro);
	NitrogenStatus.availableNitrogen -=(NitrogenStatus.leafNitrogendemand+NitrogenStatus.stemNitrogendemand+NitrogenStatus.rootNitrogendemand);
	if (NitrogenStatus.availableNitrogen < 0) NitrogenStatus.availableNitrogen = 0.;
	if (this->get_leafMass() > 0) //DHF - Add surplus N into leaf pool or stem pool
	{
		NitrogenStatus.leafNitrogenAmount += NitrogenStatus.availableNitrogen; //DHF - Adding available N into leaf class as surplus
		NitrogenStatus.availableNitrogen = 0.;
		NitrogenStatus.actualleafNitrogenRatio = NitrogenStatus.leafNitrogenAmount / (this->get_leafMass()+dblLeafgro);
	}
	else if (this->get_stemMass() > 0)
	{
		NitrogenStatus.stemNitrogenAmount += NitrogenStatus.availableNitrogen;
		NitrogenStatus.availableNitrogen = 0.;
		NitrogenStatus.actualstemNitrogenRatio = NitrogenStatus.stemNitrogenAmount / (this->get_stemMass()+dblStemgro);
	}
	this->set_TotalPlantNitrogen(NitrogenStatus.leafNitrogenAmount+NitrogenStatus.stemNitrogenAmount+NitrogenStatus.rootNitrogenAmount+NitrogenStatus.tuberNitrogenAmount);
}	

void CPlant::nitrogen_balance_posttubers2(const TInitInfo info)
{
	//Modified from the previous version to simply reallocate N.  Will not adjust partitioning coefficients*/
	// If N is limting, take the following steps
	// (1) Use N from seed reserve (no restriction on plant age or leaf area)
	// (2) Check if leaf N above the optimum level, see if this is enough to meet demand, and redistribute accordingly
	// (3) If not (2), add any N for leaf, stem or tuber above the minimum amount to an N pool for new growth
	// (4) Restrict new growth for each organ based on trying to meet minimum N amount for growth

	// Is N limiting?
	if (this->get_tubers()->isInitiated()==false) return; // skip if tubers not intiated
	double temp = 0., temp1 = 0., temp2 = 0., ratio = 0.;
	NitrogenStatus.leafNitrogenbuffer = NitrogenStatus.stemNitrogenbuffer = NitrogenStatus.tuberNitrogenbuffer = 0;
	// Calculate N demand for existing biomass and new growth
	NitrogenStatus.actualleafNitrogenRatio =0.0;
    if ((this->get_leafMass()+ this->dblLeafgro) > 0) NitrogenStatus.actualleafNitrogenRatio=NitrogenStatus.leafNitrogenAmount/(this->get_leafMass()+this->dblLeafgro);
	NitrogenStatus.actualstemNitrogenRatio =0.0;
    if ((this->get_stemMass()+this->dblStemgro) > 0) NitrogenStatus.actualstemNitrogenRatio=NitrogenStatus.stemNitrogenAmount/(this->get_stemMass()+this->dblStemgro);
    NitrogenStatus.actualrootNitrogenRatio = 0.0;
    if ((this->get_rootMass()+this->dblRootgro) > 0) NitrogenStatus.actualrootNitrogenRatio=NitrogenStatus.rootNitrogenAmount/(this->get_rootMass()+this->dblRootgro);
	NitrogenStatus.actualtuberNitrogenRatio=0.0;
	if ((this->get_tuberMass()+this->dblTubgro)>0.0) NitrogenStatus.actualtuberNitrogenRatio=NitrogenStatus.tuberNitrogenAmount/(this->get_tuberMass()+this->dblTubgro);
    NitrogenStatus.leafNitrogendemand=__max(0,(this->get_leafMass()+this->dblLeafgro)*(NitrogenStatus.optimumleafNitrogenRatio-NitrogenStatus.actualleafNitrogenRatio));
    NitrogenStatus.stemNitrogendemand=__max(0,(this->get_stemMass()+this->dblStemgro)*(NitrogenStatus.optimumstemNitrogenRatio-NitrogenStatus.actualstemNitrogenRatio));
    NitrogenStatus.rootNitrogendemand=__max(0,(this->get_rootMass()+this->dblRootgro)*(NitrogenStatus.optimumrootNitrogenRatio-NitrogenStatus.actualrootNitrogenRatio));
	NitrogenStatus.tuberNitrogendemand=__max(0,(this->get_tuberMass()+this->dblTubgro)*(NitrogenStatus.optimumtuberNitrogenRatio-NitrogenStatus.actualtuberNitrogenRatio));
	NitrogenStatus.totalNitrogendemand=NitrogenStatus.leafNitrogendemand+NitrogenStatus.stemNitrogendemand+NitrogenStatus.rootNitrogendemand+NitrogenStatus.tuberNitrogendemand;
	this->set_N_hourlygrowthdemand(NitrogenStatus.totalNitrogendemand);
	
	// Is N limiting?
	if (NitrogenStatus.availableNitrogen < NitrogenStatus.totalNitrogendemand) //(1) Add N_reserve from seed piece
	{
		if (N_reserve > 0)
		{
			double diff = NitrogenStatus.totalNitrogendemand - NitrogenStatus.availableNitrogen;
			if (diff >= N_reserve)
			{
				NitrogenStatus.availableNitrogen += N_reserve;
				this->set_N_seedused(N_reserve);
				this->set_N_hourlyseedused(N_reserve);
				N_reserve = 0.;
			}
			else
			{
				NitrogenStatus.availableNitrogen += diff;
				N_reserve -= diff;
				this->set_N_seedused(diff);
				this->set_N_hourlyseedused(diff);
			}
		}
	}
	// Check (2) Was seed N enough to meet demand?  If not look into excess leaf N buffer
	if (NitrogenStatus.availableNitrogen < NitrogenStatus.totalNitrogendemand)
	{
		if (NitrogenStatus.actualleafNitrogenRatio > NitrogenStatus.optimumleafNitrogenRatio)
		{
			NitrogenStatus.leafNitrogenbuffer = NitrogenStatus.leafNitrogenAmount-(this->get_leafMass()+this->dblLeafgro)*NitrogenStatus.optimumleafNitrogenRatio; //add any N in excess of optimum amount to leafbuffer
			if ((NitrogenStatus.availableNitrogen + NitrogenStatus.leafNitrogenbuffer) > NitrogenStatus.totalNitrogendemand) //more N taken than needed, return the rest
			{
				NitrogenStatus.leafNitrogenbuffer = NitrogenStatus.totalNitrogendemand - NitrogenStatus.availableNitrogen;
				NitrogenStatus.leafNitrogenAmount -= NitrogenStatus.leafNitrogenbuffer;
				NitrogenStatus.availableNitrogen += NitrogenStatus.leafNitrogenbuffer;
			}
			else
			{
				NitrogenStatus.availableNitrogen += NitrogenStatus.leafNitrogenbuffer;
				NitrogenStatus.leafNitrogenAmount -= NitrogenStatus.leafNitrogenbuffer;
			}
			NitrogenStatus.actualleafNitrogenRatio=NitrogenStatus.leafNitrogenAmount/(this->get_leafMass()+this->dblLeafgro);
		}
	}
	// Check (3) If Seed N and excess leaf N still not enough, start pooling N from each organ class above the minimum
	if (NitrogenStatus.availableNitrogen < NitrogenStatus.totalNitrogendemand)
	{
		//undo leafbuffer calc from prior time-step if needed
		if (NitrogenStatus.leafNitrogenbuffer > 0)
		{
			NitrogenStatus.availableNitrogen -= NitrogenStatus.leafNitrogenbuffer;
			NitrogenStatus.leafNitrogenAmount += NitrogenStatus.leafNitrogenbuffer;
			NitrogenStatus.leafNitrogenbuffer = 0.;
			NitrogenStatus.actualleafNitrogenRatio=NitrogenStatus.leafNitrogenAmount/(this->get_leafMass()+this->dblLeafgro);
		}
		if (NitrogenStatus.actualleafNitrogenRatio > NitrogenStatus.minimumleafNitrogenRatio)
		{
			NitrogenStatus.leafNitrogenbuffer = NitrogenStatus.leafNitrogenAmount - this->get_leafMass() * NitrogenStatus.minimumleafNitrogenRatio;
			NitrogenStatus.leafNitrogenAmount -= NitrogenStatus.leafNitrogenbuffer;
		}		
		if (NitrogenStatus.actualstemNitrogenRatio > NitrogenStatus.minimumstemNitrogenRatio)
		{
			NitrogenStatus.stemNitrogenbuffer = NitrogenStatus.stemNitrogenAmount - this->get_stemMass() * NitrogenStatus.minimumstemNitrogenRatio;
			NitrogenStatus.stemNitrogenAmount -= NitrogenStatus.stemNitrogenbuffer;
		}
		if (NitrogenStatus.actualtuberNitrogenRatio > NitrogenStatus.minimumtuberNitrogenRatio)
		{
			NitrogenStatus.tuberNitrogenbuffer = NitrogenStatus.tuberNitrogenAmount - this->get_tuberMass() * NitrogenStatus.minimumtuberNitrogenRatio;
			NitrogenStatus.tuberNitrogenAmount -= NitrogenStatus.tuberNitrogenbuffer;
		}
		NitrogenStatus.leafNitrogendemand = dblLeafgro * NitrogenStatus.optimumleafNitrogenRatio + __max(0,(this->get_leafMass() * (NitrogenStatus.minimumleafNitrogenRatio - NitrogenStatus.actualleafNitrogenRatio)));
		NitrogenStatus.stemNitrogendemand = dblStemgro * NitrogenStatus.optimumstemNitrogenRatio + __max(0,(this->get_stemMass() * (NitrogenStatus.minimumstemNitrogenRatio - NitrogenStatus.actualstemNitrogenRatio)));
		NitrogenStatus.tuberNitrogendemand = dblTubgro * NitrogenStatus.optimumtuberNitrogenRatio + __max(0,(this->get_tuberMass()*(NitrogenStatus.minimumtuberNitrogenRatio- NitrogenStatus.actualtuberNitrogenRatio)));
		NitrogenStatus.availableNitrogen += NitrogenStatus.leafNitrogenbuffer + NitrogenStatus.stemNitrogenbuffer + NitrogenStatus.tuberNitrogenbuffer;
		NitrogenStatus.totalNitrogendemand = NitrogenStatus.leafNitrogendemand + NitrogenStatus.stemNitrogendemand + NitrogenStatus.rootNitrogendemand + NitrogenStatus.tuberNitrogendemand;
		if (NitrogenStatus.availableNitrogen >= NitrogenStatus.totalNitrogendemand) //now have enough demand for new growth and more - redistribute N back to pools
		{
			double redistribute = NitrogenStatus.availableNitrogen - NitrogenStatus.totalNitrogendemand;

			//redistribute N to leaf first, then stem, then tuber
			if (redistribute > (NitrogenStatus.leafNitrogenbuffer + NitrogenStatus.stemNitrogenbuffer))
			{
				NitrogenStatus.leafNitrogendemand += NitrogenStatus.leafNitrogenbuffer;
				NitrogenStatus.stemNitrogendemand += NitrogenStatus.stemNitrogenbuffer;
				if ((this->get_tuberMass()+dblTubgro) > 0)
				{ //check if tubers have actually grown since initiation
					NitrogenStatus.tuberNitrogendemand += (redistribute - NitrogenStatus.leafNitrogenbuffer - NitrogenStatus.stemNitrogenbuffer);
				}
				else //otherwise, put back into canopy using weighted average of N demand for leaves and stems
				{
					temp1 = (this->get_leafMass()+this->dblLeafgro) *(NitrogenStatus.optimumleafNitrogenRatio-NitrogenStatus.actualleafNitrogenRatio);
					temp2 = (this->get_stemMass()+this->dblStemgro) *(NitrogenStatus.optimumstemNitrogenRatio-NitrogenStatus.actualstemNitrogenRatio);
					temp = redistribute - NitrogenStatus.leafNitrogenbuffer - NitrogenStatus.stemNitrogenbuffer;
					NitrogenStatus.leafNitrogendemand += (temp1 / (temp1+temp2)*temp);
					NitrogenStatus.stemNitrogendemand += (temp2 / (temp1+temp2)*temp);
				}

				this->set_N_hourlytranslocated(redistribute - NitrogenStatus.totalNitrogendemand - NitrogenStatus.leafNitrogenbuffer - NitrogenStatus.stemNitrogenbuffer);
			}
			else if (redistribute > NitrogenStatus.leafNitrogenbuffer)
			{
				NitrogenStatus.leafNitrogendemand += NitrogenStatus.leafNitrogenbuffer;
				NitrogenStatus.stemNitrogendemand += (redistribute - NitrogenStatus.leafNitrogenbuffer);
				this->set_N_hourlytranslocated(redistribute - NitrogenStatus.totalNitrogendemand - NitrogenStatus.leafNitrogenbuffer);
			}
			else
			{
				NitrogenStatus.leafNitrogendemand += redistribute;
				this->set_N_hourlytranslocated(redistribute - NitrogenStatus.totalNitrogendemand);
			}
			NitrogenStatus.totalNitrogendemand = NitrogenStatus.leafNitrogendemand+NitrogenStatus.stemNitrogendemand+NitrogenStatus.rootNitrogendemand+NitrogenStatus.tuberNitrogendemand;
		}
		else //will need to reduce new growth accordingly based on N available for each organ since there's still not enough available to meet new demand
		{
			ratio = (NitrogenStatus.availableNitrogen - NitrogenStatus.rootNitrogendemand)/ (NitrogenStatus.totalNitrogendemand - NitrogenStatus.rootNitrogendemand); 
			if (ratio > 0)
			{
				NitrogenStatus.leafNitrogendemand *= ratio;
				NitrogenStatus.stemNitrogendemand *= ratio;
				NitrogenStatus.tuberNitrogendemand *= ratio;
			}
			else //* assume root growth has 100% priority for now unless leaf N content dropped below minimum due to senescence
			{
				if (NitrogenStatus.availableNitrogen > __max(0.,this->get_leafMass()*(NitrogenStatus.minimumleafNitrogenRatio-NitrogenStatus.actualleafNitrogenRatio)))
				{
					NitrogenStatus.rootNitrogendemand = NitrogenStatus.availableNitrogen - __max(0.,this->get_leafMass()*(NitrogenStatus.minimumleafNitrogenRatio-NitrogenStatus.actualleafNitrogenRatio));
					NitrogenStatus.leafNitrogendemand = 0.;
					NitrogenStatus.stemNitrogendemand = 0.;
					NitrogenStatus.tuberNitrogendemand = 0.;
				}
				else  //this should force 0 root grown and only allocate N to existing leaf mass at end of routine
				{
					NitrogenStatus.rootNitrogendemand = 0.;
					NitrogenStatus.leafNitrogendemand = 0.;
					NitrogenStatus.stemNitrogendemand = 0.;
					NitrogenStatus.tuberNitrogendemand = 0.;
				}
				//NitrogenStatus.rootNitrogendemand = NitrogenStatus.availableNitrogen;
			}
			dblLeafgro = NitrogenStatus.leafNitrogendemand / NitrogenStatus.optimumleafNitrogenRatio;
			dblStemgro = NitrogenStatus.stemNitrogendemand / NitrogenStatus.optimumstemNitrogenRatio;
			dblTubgro = NitrogenStatus.tuberNitrogendemand / NitrogenStatus.optimumtuberNitrogenRatio;
			dblRootgro = NitrogenStatus.rootNitrogendemand / NitrogenStatus.optimumrootNitrogenRatio;
			this->set_N_hourlytranslocated(NitrogenStatus.leafNitrogenbuffer + NitrogenStatus.stemNitrogenbuffer + NitrogenStatus.tuberNitrogenbuffer);

			/*
			dblLeafgro = ratio * dblLeafgro / (dblLeafgro + dblStemgro + dblTubgro + dblRootgro); //will this work?
			dblStemgro = ratio * dblStemgro / (dblLeafgro + dblStemgro + dblTubgro + dblRootgro);
			dblTubgro = ratio * dblStemgro / (dblLeafgro + dblStemgro + dblTubgro + dblRootgro);
			NitrogenStatus.leafNitrogendemand = dblLeafgro * NitrogenStatus.optimumleafNitrogenRatio;
			NitrogenStatus.stemNitrogendemand = dblStemgro * NitrogenStatus.optimumstemNitrogenRatio;
			NitrogenStatus.tuberNitrogendemand = dblTubgro * NitrogenStatus.optimumtuberNitrogenRatio;
			NitrogenStatus.totalNitrogendemand = NitrogenStatus.leafNitrogendemand+NitrogenStatus.stemNitrogendemand+NitrogenStatus.rootNitrogendemand+NitrogenStatus.tuberNitrogendemand;
			*/
		}
	}
	NitrogenStatus.leafNitrogenAmount += NitrogenStatus.leafNitrogendemand;
	NitrogenStatus.stemNitrogenAmount += NitrogenStatus.stemNitrogendemand;
	NitrogenStatus.rootNitrogenAmount += NitrogenStatus.rootNitrogendemand;
	NitrogenStatus.tuberNitrogenAmount += NitrogenStatus.tuberNitrogendemand;
	NitrogenStatus.totalNitrogendemand=NitrogenStatus.leafNitrogendemand+NitrogenStatus.stemNitrogendemand+NitrogenStatus.rootNitrogendemand+NitrogenStatus.tuberNitrogendemand;
	NitrogenStatus.availableNitrogen -= (NitrogenStatus.leafNitrogendemand + NitrogenStatus.stemNitrogendemand + NitrogenStatus.rootNitrogendemand + NitrogenStatus.tuberNitrogendemand);
	if ((this->get_leafMass() + dblLeafgro) > 0)
		NitrogenStatus.actualleafNitrogenRatio = NitrogenStatus.leafNitrogenAmount/(this->get_leafMass() + dblLeafgro);
	else
		NitrogenStatus.actualleafNitrogenRatio = 0.; 
	if ((this->get_stemMass() + dblStemgro) > 0)
		NitrogenStatus.actualstemNitrogenRatio = NitrogenStatus.stemNitrogenAmount/(this->get_stemMass() + dblStemgro);
	if ((this->get_tuberMass() + dblTubgro) > 0)
		NitrogenStatus.actualtuberNitrogenRatio = NitrogenStatus.tuberNitrogenAmount/(this->get_tuberMass() + dblTubgro);
	NitrogenStatus.actualrootNitrogenRatio = NitrogenStatus.rootNitrogenAmount / (this->get_rootMass() + dblRootgro);
	if (NitrogenStatus.availableNitrogen < 0) NitrogenStatus.availableNitrogen = 0.;
	if (this->get_leafMass() > 0 && NitrogenStatus.availableNitrogen > 0) //DHF - Add surplus N into leaf pool or stem pool, up to 8% N in leaf, stem
	{
		double surplusNstorage = (0.08 - NitrogenStatus.actualleafNitrogenRatio)*(this->dblLeafgro + this->get_leafMass());
		if (surplusNstorage >= NitrogenStatus.availableNitrogen) // add full amount to leaves
		{
			NitrogenStatus.leafNitrogendemand += NitrogenStatus.availableNitrogen;
			NitrogenStatus.leafNitrogenAmount += NitrogenStatus.availableNitrogen;
			NitrogenStatus.availableNitrogen = 0.;
		}
		else if (surplusNstorage > 0 && surplusNstorage < NitrogenStatus.availableNitrogen) //add partial amount to leaves
		{
			NitrogenStatus.leafNitrogendemand += surplusNstorage;
			NitrogenStatus.leafNitrogenAmount += surplusNstorage;
			NitrogenStatus.availableNitrogen -= surplusNstorage;
		}
		NitrogenStatus.actualleafNitrogenRatio = NitrogenStatus.leafNitrogenAmount / (this->dblLeafgro + this->get_leafMass());
		if (NitrogenStatus.actualleafNitrogenRatio > 0.09)
		{
			double excessN = (this->dblLeafgro + this->get_leafMass())*(NitrogenStatus.actualleafNitrogenRatio - 0.09);
			NitrogenStatus.leafNitrogenAmount -= excessN;
			NitrogenStatus.leafNitrogendemand -= excessN;
			NitrogenStatus.availableNitrogen += excessN;
			NitrogenStatus.actualleafNitrogenRatio = NitrogenStatus.leafNitrogenAmount / (this->dblLeafgro + this->get_leafMass());
		}
		/*
		NitrogenStatus.leafNitrogenAmount += NitrogenStatus.availableNitrogen;
		NitrogenStatus.availableNitrogen = 0;
		NitrogenStatus.actualleafNitrogenRatio = NitrogenStatus.leafNitrogenAmount / (this->dblLeafgro + this->get_leafMass());
		*/
	}
	if (this->get_stemMass() > 0 && NitrogenStatus.availableNitrogen > 0)
	{
		double surplusNstorage = (0.08 - NitrogenStatus.actualstemNitrogenRatio)*(this->dblStemgro + this->get_stemMass());
		if (surplusNstorage >= NitrogenStatus.availableNitrogen) // add full amount to leaves
		{
			NitrogenStatus.stemNitrogendemand += NitrogenStatus.availableNitrogen;
			NitrogenStatus.stemNitrogenAmount += NitrogenStatus.availableNitrogen;
			NitrogenStatus.availableNitrogen = 0.;
		}
		else if (surplusNstorage > 0 && surplusNstorage < NitrogenStatus.availableNitrogen) //add partial amount to leaves
		{
			NitrogenStatus.stemNitrogendemand += surplusNstorage;
			NitrogenStatus.stemNitrogenAmount += surplusNstorage;
			NitrogenStatus.availableNitrogen -= surplusNstorage;
		}
		NitrogenStatus.actualstemNitrogenRatio = NitrogenStatus.stemNitrogenAmount / (this->dblStemgro + this->get_stemMass());
	}
	if (this->get_tuberMass() > 0 && NitrogenStatus.availableNitrogen >0)
	{
		double surplusNstorage = (0.05 - NitrogenStatus.actualtuberNitrogenRatio)*(this->dblTubgro + this->get_tuberMass());
		if (surplusNstorage >= NitrogenStatus.availableNitrogen) // add full amount to leaves
		{
			NitrogenStatus.tuberNitrogendemand += NitrogenStatus.availableNitrogen;
			NitrogenStatus.tuberNitrogenAmount += NitrogenStatus.availableNitrogen;
			NitrogenStatus.availableNitrogen = 0.;
		}
		else if (surplusNstorage > 0 && surplusNstorage < NitrogenStatus.availableNitrogen) //add partial amount to leaves
		{
			NitrogenStatus.tuberNitrogendemand += surplusNstorage;
			NitrogenStatus.tuberNitrogenAmount += surplusNstorage;
			NitrogenStatus.availableNitrogen -= surplusNstorage;
		}
		NitrogenStatus.actualtuberNitrogenRatio = NitrogenStatus.tuberNitrogenAmount / (this->dblTubgro + this->get_tuberMass());
	}
	if (this->get_rootMass() > 0 && NitrogenStatus.availableNitrogen >0)
	{
		double surplusNstorage = (0.1 - NitrogenStatus.actualrootNitrogenRatio)*(this->dblRootgro + this->get_rootMass());
		if (surplusNstorage >= NitrogenStatus.availableNitrogen) // add full amount to leaves
		{
			NitrogenStatus.rootNitrogendemand += NitrogenStatus.availableNitrogen;
			NitrogenStatus.rootNitrogenAmount += NitrogenStatus.availableNitrogen;
			NitrogenStatus.availableNitrogen = 0.;
		}
		else if (surplusNstorage > 0 && surplusNstorage < NitrogenStatus.availableNitrogen) //add partial amount to leaves
		{
			NitrogenStatus.rootNitrogendemand += surplusNstorage;
			NitrogenStatus.rootNitrogenAmount += surplusNstorage;
			NitrogenStatus.availableNitrogen -= surplusNstorage;
		}
		NitrogenStatus.actualrootNitrogenRatio = NitrogenStatus.rootNitrogenAmount / (this->dblRootgro + this->get_rootMass());
	}
	this->set_TotalPlantNitrogen(NitrogenStatus.leafNitrogenAmount+NitrogenStatus.stemNitrogenAmount+NitrogenStatus.rootNitrogenAmount+NitrogenStatus.tuberNitrogenAmount);
}

void CPlant::nitrogen_balance_posttubers_substor(const TInitInfo info)
{
//***********************************
//* Coefficients come from SUBSTOR - more reasonable, allocation routine similar to SUBSTOR but modified here
//* Main idea is that stem and leaf N demand is identical and both can serve as N buffer

	double leafNdemand = 0., stemNdemand = 0., rootNdemand = 0., tuberNdemand = 0., surplusN = 0.;
	double leafnewgrowthNdemand = 0., stemnewgrowthNdemand = 0., rootnewgrowthNdemand = 0., tubernewgrowthNdemand = 0., newgrowthNdemand = 0.;
	double leafoldgrowthNdemand = 0., stemoldgrowthNdemand = 0., rootoldgrowthNdemand = 0., tuberoldgrowthNdemand = 0., oldgrowthNdemand = 0.;
	if (!this->get_tubers()->isInitiated()) return;
	NitrogenStatus.actualleafNitrogenRatio = 0.;
	NitrogenStatus.actualstemNitrogenRatio = 0.;
	NitrogenStatus.actualrootNitrogenRatio = 0.;
	NitrogenStatus.actualtuberNitrogenRatio = 0.;
	if ((this->get_leafMass()) > 0) NitrogenStatus.actualleafNitrogenRatio = NitrogenStatus.leafNitrogenAmount / (this->get_leafMass()+dblLeafgro);
	if ((this->get_stemMass()) > 0) NitrogenStatus.actualstemNitrogenRatio = NitrogenStatus.stemNitrogenAmount / (this->get_stemMass()+dblStemgro);
	if ((this->get_rootMass()) > 0) NitrogenStatus.actualrootNitrogenRatio = NitrogenStatus.rootNitrogenAmount / (this->get_rootMass()+dblRootgro);
	if ((this->get_tuberMass()) > 0) NitrogenStatus.actualtuberNitrogenRatio = NitrogenStatus.tuberNitrogenAmount / (this->get_tuberMass()+dblTubgro);
	NitrogenStatus.leafNitrogendemand = (this->get_leafMass()+this->dblLeafgro)*(NitrogenStatus.optimumleafNitrogenRatio - NitrogenStatus.actualleafNitrogenRatio);
	NitrogenStatus.stemNitrogendemand = (this->get_stemMass()+this->dblStemgro)*(NitrogenStatus.optimumstemNitrogenRatio - NitrogenStatus.actualstemNitrogenRatio);
	NitrogenStatus.rootNitrogendemand = (this->get_rootMass()+this->dblRootgro)*(NitrogenStatus.optimumrootNitrogenRatio - NitrogenStatus.actualrootNitrogenRatio);
	NitrogenStatus.tuberNitrogendemand = (this->get_tuberMass()+this->dblTubgro)*(NitrogenStatus.optimumtuberNitrogenRatio - NitrogenStatus.actualtuberNitrogenRatio);
	NitrogenStatus.totalNitrogendemand = NitrogenStatus.leafNitrogendemand + NitrogenStatus.stemNitrogendemand + NitrogenStatus.rootNitrogendemand + NitrogenStatus.tuberNitrogendemand;
	
	NitrogenStatus.leafNitrogenbuffer = __max(((NitrogenStatus.actualleafNitrogenRatio - NitrogenStatus.minimumleafNitrogenRatio)*this->get_leafMass()),0.);
	NitrogenStatus.stemNitrogenbuffer = __max(((NitrogenStatus.actualstemNitrogenRatio - NitrogenStatus.minimumstemNitrogenRatio)*this->get_stemMass()),0.);
	NitrogenStatus.tuberNitrogenbuffer = __max(((NitrogenStatus.actualtuberNitrogenRatio - NitrogenStatus.minimumtuberNitrogenRatio)*this->get_tuberMass()),0.);
	leafnewgrowthNdemand = this->dblLeafgro * NitrogenStatus.optimumleafNitrogenRatio;
	stemnewgrowthNdemand = this->dblStemgro * NitrogenStatus.optimumstemNitrogenRatio;
	rootnewgrowthNdemand = this->dblRootgro * NitrogenStatus.optimumrootNitrogenRatio;
	tubernewgrowthNdemand = this->dblTubgro * NitrogenStatus.optimumtuberNitrogenRatio;
	newgrowthNdemand = leafnewgrowthNdemand + stemnewgrowthNdemand + rootnewgrowthNdemand + tubernewgrowthNdemand; 	
	
	leafoldgrowthNdemand = NitrogenStatus.leafNitrogendemand - leafnewgrowthNdemand;
	stemoldgrowthNdemand = NitrogenStatus.stemNitrogendemand - stemnewgrowthNdemand;
	rootoldgrowthNdemand = NitrogenStatus.rootNitrogendemand - rootnewgrowthNdemand;
	tuberoldgrowthNdemand = NitrogenStatus.tuberNitrogendemand - tubernewgrowthNdemand;
	oldgrowthNdemand = leafoldgrowthNdemand + stemoldgrowthNdemand + rootoldgrowthNdemand + tuberoldgrowthNdemand;

	// scenarios
	// 0) N in excess, store in haulm and tubers - will need to develop better scheme for luxury consumption
	if (NitrogenStatus.availableNitrogen >= NitrogenStatus.totalNitrogendemand)
	{
		//store excess in tops and tubers
		surplusN = (NitrogenStatus.availableNitrogen-NitrogenStatus.totalNitrogendemand);
		leafNdemand = __max(0.,NitrogenStatus.leafNitrogendemand) + 0.4*surplusN;
		stemNdemand = __max(0.,NitrogenStatus.stemNitrogendemand) + 0.4*surplusN;
		rootNdemand = __max(0.,NitrogenStatus.rootNitrogendemand);
		tuberNdemand = __max(0.,NitrogenStatus.tuberNitrogendemand) + 0.2*surplusN;
		surplusN = 0.;
	}
	else
// 2) Is there enough N to satisfy new growth?
	{
// 2a) Yes, satisfy new growth, then allocate any extra N to leaves, stems, then tubers
		if (NitrogenStatus.availableNitrogen >= newgrowthNdemand)
		{
			surplusN = NitrogenStatus.availableNitrogen - newgrowthNdemand;
			leafNdemand = leafnewgrowthNdemand;
			stemNdemand = stemnewgrowthNdemand;
			rootNdemand = rootnewgrowthNdemand;
			tuberNdemand = tubernewgrowthNdemand;
			if (surplusN > 0) // allocate remaining N to either resotre leaf N content, stem N content, or tuber N content to optimum in that order
			{
				if ((stemoldgrowthNdemand + leafoldgrowthNdemand) < surplusN)
				{
					leafNdemand += leafoldgrowthNdemand;
					surplusN -= leafoldgrowthNdemand;
					stemNdemand += stemoldgrowthNdemand;
					surplusN -= stemoldgrowthNdemand;
					tuberNdemand += surplusN;
					surplusN = 0.;
				}
				else if ((stemoldgrowthNdemand + leafoldgrowthNdemand) >= surplusN)
				{
					if (leafoldgrowthNdemand > surplusN)
					{
						leafNdemand += surplusN;
						surplusN = 0.;
					}
					else
					{
						leafNdemand += leafoldgrowthNdemand;
						surplusN -= leafoldgrowthNdemand;
						if (stemoldgrowthNdemand > surplusN)
						{
							stemNdemand += surplusN;
							surplusN = 0.;
						}
						else
						{
							stemNdemand += stemoldgrowthNdemand;
							surplusN -= stemoldgrowthNdemand;
							tuberNdemand += surplusN;
							surplusN = 0.;
						}
					}
				}
			}
		}
	// 2b) No, so give tuber N demand first priority, then adjust top growth based on availability of top N buffer
		else if (NitrogenStatus.availableNitrogen >= tubernewgrowthNdemand)
		{
			double extraN = NitrogenStatus.availableNitrogen - tubernewgrowthNdemand;
			surplusN = __min(0.05*(NitrogenStatus.leafNitrogenbuffer + NitrogenStatus.stemNitrogenbuffer),(leafnewgrowthNdemand + stemnewgrowthNdemand + rootnewgrowthNdemand - extraN));
			NitrogenStatus.stemNitrogenAmount -= 0.75*surplusN;
			NitrogenStatus.leafNitrogenAmount -= 0.25*surplusN;
			NitrogenStatus.availableNitrogen += surplusN;
	// 2c) Enough available N and haulm N to satisfy all new vegetative growth
			if ((extraN + surplusN) >= (leafnewgrowthNdemand+stemnewgrowthNdemand+rootnewgrowthNdemand))
			{
				leafNdemand = leafnewgrowthNdemand;
				stemNdemand = stemnewgrowthNdemand;
				rootNdemand = rootnewgrowthNdemand;
				tuberNdemand = tubernewgrowthNdemand;
				surplusN = 0.;
			}
	// 2d) Not enough available N and haulm N to satisfy all new vegetative growth, so growth tubers accordingly, reduce growth of other organs
			else
			{
				if ((leafnewgrowthNdemand + stemnewgrowthNdemand + rootnewgrowthNdemand) > 0)
				{
					double factor = (extraN + surplusN) / (leafnewgrowthNdemand + stemnewgrowthNdemand + rootnewgrowthNdemand);
					this->dblLeafgro *= factor;
					this->dblStemgro *= factor;
					this->dblRootgro *= factor;
					leafNdemand = this->dblLeafgro * NitrogenStatus.optimumleafNitrogenRatio;
					stemNdemand = this->dblStemgro * NitrogenStatus.optimumstemNitrogenRatio;
					rootNdemand = this->dblRootgro * NitrogenStatus.optimumrootNitrogenRatio;
					tuberNdemand = tubernewgrowthNdemand;
					surplusN = 0.;
				}
			}
		}
	// 3) Not enough AvailableN to satisfy new tuber growth, so just use haulm N for hits purpose, and limit growth to other organs
		else if (NitrogenStatus.availableNitrogen < tubernewgrowthNdemand)
		{
			this->dblLeafgro = 0.;
			this->dblStemgro = 0.;
			this->dblRootgro = 0.;
			leafNdemand = 0.;
			stemNdemand = 0.;
			rootNdemand = 0.;
			surplusN = 0.05 * (NitrogenStatus.leafNitrogenbuffer + NitrogenStatus.stemNitrogenbuffer);
			NitrogenStatus.leafNitrogenAmount -= 0.25 * surplusN;
			NitrogenStatus.stemNitrogenAmount -= 0.75 * surplusN;
			NitrogenStatus.availableNitrogen += surplusN;
	// 3a) Still can't meet new tuber growth even with haulm N
			if ((NitrogenStatus.availableNitrogen) < tubernewgrowthNdemand)
			{
				this->dblTubgro *= NitrogenStatus.availableNitrogen / tubernewgrowthNdemand;
				tuberNdemand = dblTubgro * NitrogenStatus.optimumtuberNitrogenRatio;
			}
			else
			{
				tuberNdemand = tubernewgrowthNdemand + __max(0.,(NitrogenStatus.availableNitrogen - tubernewgrowthNdemand));
				surplusN = 0.;
			}
		}
	}
	NitrogenStatus.leafNitrogendemand = leafNdemand;
	NitrogenStatus.stemNitrogendemand = stemNdemand;
	NitrogenStatus.rootNitrogendemand = rootNdemand;
	NitrogenStatus.tuberNitrogendemand = tuberNdemand;
	NitrogenStatus.leafNitrogenAmount += NitrogenStatus.leafNitrogendemand;
	NitrogenStatus.stemNitrogenAmount += NitrogenStatus.stemNitrogendemand;
	NitrogenStatus.rootNitrogenAmount += NitrogenStatus.rootNitrogendemand;
	NitrogenStatus.tuberNitrogenAmount += NitrogenStatus.tuberNitrogendemand;
	NitrogenStatus.totalNitrogendemand=NitrogenStatus.leafNitrogendemand+NitrogenStatus.stemNitrogendemand+NitrogenStatus.rootNitrogendemand+NitrogenStatus.tuberNitrogendemand;
	NitrogenStatus.availableNitrogen -= (NitrogenStatus.leafNitrogendemand + NitrogenStatus.stemNitrogendemand + NitrogenStatus.rootNitrogendemand + NitrogenStatus.tuberNitrogendemand);
	if ((this->get_leafMass() + dblLeafgro) > 0)
		NitrogenStatus.actualleafNitrogenRatio = NitrogenStatus.leafNitrogenAmount/(this->get_leafMass() + dblLeafgro);
	else
		NitrogenStatus.actualleafNitrogenRatio = 0.; 
	if ((this->get_stemMass() + dblStemgro) > 0)
		NitrogenStatus.actualstemNitrogenRatio = NitrogenStatus.stemNitrogenAmount/(this->get_stemMass() + dblStemgro);
	if ((this->get_tuberMass() + dblTubgro) > 0)
		NitrogenStatus.actualtuberNitrogenRatio = NitrogenStatus.tuberNitrogenAmount/(this->get_tuberMass() + dblTubgro);
	NitrogenStatus.actualrootNitrogenRatio = NitrogenStatus.rootNitrogenAmount / (this->get_rootMass() + dblRootgro);
	if (NitrogenStatus.availableNitrogen < 0) NitrogenStatus.availableNitrogen = 0.;
	if (this->get_leafMass() > 0) //DHF - Add surplus N into leaf pool or stem pool
	{
		NitrogenStatus.leafNitrogenAmount += NitrogenStatus.availableNitrogen; //DHF - Adding available N into leaf class as surplus
		NitrogenStatus.availableNitrogen = 0.;
		NitrogenStatus.actualleafNitrogenRatio = NitrogenStatus.leafNitrogenAmount / (this->get_leafMass()+dblLeafgro);
	}
	else if (this->get_stemMass() > 0)
	{
		NitrogenStatus.stemNitrogenAmount += NitrogenStatus.availableNitrogen;
		NitrogenStatus.availableNitrogen = 0.;
		NitrogenStatus.actualstemNitrogenRatio = NitrogenStatus.stemNitrogenAmount / (this->get_stemMass()+dblStemgro);
	}
	this->set_TotalPlantNitrogen(NitrogenStatus.leafNitrogenAmount+NitrogenStatus.stemNitrogenAmount+NitrogenStatus.rootNitrogenAmount+NitrogenStatus.tuberNitrogenAmount);

}


void CPlant::nitrogen_balance_posttubers(const TInitInfo info)
{
// Routine derived from SIMPOTATO PARTTN subroutines
// Difference between subroutines is based on developmental stage of the plant (post TI)
// Routine adjusts actual C growth based on N deficiency - therefore, must run after C_Allocation1 method.
// Simulates whole N demand by haulm and roots based on potential growth

// First calculate N demand for existing biomass and potential growth
// Compare N demand to total N uptake
// If N is limiting, take the following steps in order to reduce N
//    demand until availn meets demand:
// (1) Use N from seed reserve if any is left
// (2) If leaf N content is above the critical or optimum N level (i.e. prior luxury consumption of N) then redistribute N.
// (3) if leaf, stem or tuber N concentration is more than the minimum amount, reduce to minimum and add freed up N to available N for new growth
// (4) Shift growth from leaves to stems and tubers according to an N:C ratio from tabular form
// (5) Remove C from leaves to free more N for tuber growth
// Right now, subroutine is independent of bigleaf or individual leaf model as N demand based on total leaf or stem mass
// ********************************

	if (this->get_tubers()->isInitiated()==false) return; //only after TI
	bool tablelocation = false;
	double temp = 0., ratio = 0., ncratio = 0.;
	double availC = 0., availN = 0., xc = 0., xn = 0.;
	temp = 0.;
	double table[23][4] = {{0.0275,0,0,100},{0.028,0.5,0.5,99},{0.03,1.5,1.5,97},{0.032,3,3,94},{0.034,5,5,90},{0.036,7.5,7.5,85},
	{0.038,10.5,10.5,79},{0.04,13.5,13.5,73},{0.042,17.5,15.5,67},{0.044,22,17,61},{0.046,27,18,55},{0.048,32,19,49},{0.05,37,20,43},{0.052,42,21,37},{0.054,47,22,31},
	{0.056,50,22.5,27.5},{0.058,53,23,24},{0.06,56,23,21},{0.062,60,23,17},{0.064,65,23,12},{0.066,69,23,8},{0.068,74,23,3},{0.06929,77,23,0}}; //will use elements 0 to 22, 0 to 3 vs SIMPOTATO 1 thruogh 23, 1 through 4

	//   CALCULATE N DEMAND FOR EXISTING BIOMASS AND NEW GROWTH
    NitrogenStatus.actualleafNitrogenRatio =0.0;
    if ((this->get_leafMass()+ this->dblLeafgro) > 0) NitrogenStatus.actualleafNitrogenRatio=NitrogenStatus.leafNitrogenAmount/(this->get_leafMass()+this->dblLeafgro);
	NitrogenStatus.actualstemNitrogenRatio =0.0;
    if ((this->get_stemMass()+this->dblStemgro) > 0) NitrogenStatus.actualstemNitrogenRatio=NitrogenStatus.stemNitrogenAmount/(this->get_stemMass()+this->dblStemgro);
    NitrogenStatus.actualrootNitrogenRatio = 0.0;
    if ((this->get_rootMass()+this->dblRootgro) > 0) NitrogenStatus.actualrootNitrogenRatio=NitrogenStatus.rootNitrogenAmount/(this->get_rootMass()+this->dblRootgro);
	NitrogenStatus.actualtuberNitrogenRatio=0.0;
	if ((this->get_tuberMass()+this->dblTubgro)>0.0) NitrogenStatus.actualtuberNitrogenRatio=NitrogenStatus.tuberNitrogenAmount/(this->get_tuberMass()+this->dblTubgro);
    NitrogenStatus.leafNitrogendemand=(this->get_leafMass()+this->dblLeafgro)*(NitrogenStatus.optimumleafNitrogenRatio-NitrogenStatus.actualleafNitrogenRatio);
    NitrogenStatus.stemNitrogendemand=(this->get_stemMass()+this->dblStemgro)*(NitrogenStatus.optimumstemNitrogenRatio-NitrogenStatus.actualstemNitrogenRatio);
    NitrogenStatus.rootNitrogendemand=(this->get_rootMass()+this->dblRootgro)*(NitrogenStatus.optimumrootNitrogenRatio-NitrogenStatus.actualrootNitrogenRatio);
	NitrogenStatus.tuberNitrogendemand=(this->get_tuberMass()+this->dblTubgro)*(NitrogenStatus.optimumtuberNitrogenRatio-NitrogenStatus.actualtuberNitrogenRatio);
//	if (NitrogenStatus.leafNitrogendemand < 0) NitrogenStatus.leafNitrogendemand =0.0;
//    if (NitrogenStatus.stemNitrogendemand < 0) NitrogenStatus.stemNitrogendemand =0.0;
//	if (NitrogenStatus.rootNitrogendemand < 0) NitrogenStatus.rootNitrogendemand =0.0;
//    if (NitrogenStatus.tuberNitrogendemand < 0) NitrogenStatus.tuberNitrogendemand =0.0;
    NitrogenStatus.totalNitrogendemand=NitrogenStatus.leafNitrogendemand+NitrogenStatus.stemNitrogendemand+NitrogenStatus.rootNitrogendemand+NitrogenStatus.tuberNitrogendemand;
    //  IF N STATUS IS ADEQUATE, RETURN.  NO N WILL BE TAKEN FROM THE SOIL
//    if (NitrogenStatus.totalNitrogendemand <= 0) return;
	//  if excess nitrogen is available, put some of it in leaves; the rest will not be taken from soil
  		// 0) First, if N demand not met, simply add any available seedpiece reserve to plant N supply as long as plant is still dependent on seedpiece = DHF
	if (NitrogenStatus.availableNitrogen < NitrogenStatus.totalNitrogendemand)
	{
		if (this->get_greenLeafArea() < 10000 && N_reserve > 0 && (NitrogenStatus.optimumleafNitrogenRatio > NitrogenStatus.actualleafNitrogenRatio)) 
		// 400 cm2 leaf area same restriction as seedpiece CHO in C_allocation1
		// but N availabe from seed only if N in leaf is below the optimal level
		{
			double diff = NitrogenStatus.totalNitrogendemand - NitrogenStatus.availableNitrogen;
			if (diff >= N_reserve)
			{
				NitrogenStatus.availableNitrogen += N_reserve;
				this->set_N_seedused(N_reserve);
				N_reserve = 0.;
			}
			else 
			{
				NitrogenStatus.availableNitrogen += diff;
				N_reserve -= diff;
				this->set_N_seedused(diff);
			}
		}
	}
	
	if (NitrogenStatus.availableNitrogen >= NitrogenStatus.totalNitrogendemand)
	{
		if (NitrogenStatus.totalNitrogendemand > 0) { //following if-then needed to make sure only current surplus N is added to temp variable
			temp = NitrogenStatus.availableNitrogen - NitrogenStatus.totalNitrogendemand; 
		}
		else
		{
			if (NitrogenStatus.leafNitrogendemand > 0) temp -= NitrogenStatus.leafNitrogendemand;
			if (NitrogenStatus.stemNitrogendemand > 0) temp -= NitrogenStatus.stemNitrogendemand;
			if (NitrogenStatus.rootNitrogendemand > 0) temp -= NitrogenStatus.rootNitrogendemand;
			if (NitrogenStatus.tuberNitrogendemand > 0) temp -= NitrogenStatus.tuberNitrogendemand;
			temp += NitrogenStatus.availableNitrogen;
		}
        // if (temp > NitrogenStatus.leafNitrogendemand*0.2) temp=NitrogenStatus.leafNitrogendemand*0.2;//Originaly from SIMPOTATO, idea was to limit uptake from soil if leaf N > 1.2 times demand; N uptake controlled my m-m kinetics in spudsim, DHF
         NitrogenStatus.leafNitrogendemand += temp;
         NitrogenStatus.totalNitrogendemand=NitrogenStatus.leafNitrogendemand+NitrogenStatus.stemNitrogendemand+NitrogenStatus.rootNitrogendemand+NitrogenStatus.tuberNitrogendemand;
	}
	else
	{

	//   IF SURPLUS N IN LEAVES ABOVE TCNP, MOVE TO AVAILN
         if (NitrogenStatus.actualleafNitrogenRatio > NitrogenStatus.optimumleafNitrogenRatio)
		 {
            NitrogenStatus.leafNitrogenbuffer=NitrogenStatus.leafNitrogenAmount-(this->get_leafMass()+this->dblLeafgro)*NitrogenStatus.optimumleafNitrogenRatio;
            if ((NitrogenStatus.availableNitrogen+NitrogenStatus.leafNitrogenbuffer)>=NitrogenStatus.totalNitrogendemand)
			{
				NitrogenStatus.leafNitrogenAmount -= (NitrogenStatus.totalNitrogendemand - NitrogenStatus.availableNitrogen);
				NitrogenStatus.availableNitrogen = NitrogenStatus.totalNitrogendemand;
			}
			else
			{
				NitrogenStatus.leafNitrogenAmount -= NitrogenStatus.leafNitrogenbuffer;
				NitrogenStatus.availableNitrogen += NitrogenStatus.leafNitrogenbuffer;
			}
			NitrogenStatus.actualleafNitrogenRatio = 0.;
			if ((this->get_leafMass()+this->dblLeafgro) > 0) NitrogenStatus.actualleafNitrogenRatio=NitrogenStatus.leafNitrogenAmount/(this->get_leafMass()+this->dblLeafgro);
		 }
	//   IF AVAILN STILL < ANDEM AND IF SURPLUS N IN LEAVES, STEMS, OR TUBERS 
	//     IS ABOVE MINIMUM, MOVE UP TO 4% TO AVAILN.  IT WILL BE DISTRIBUTED 
	//     AMONG THEM ACCORDING TO GROWTH.
		 if (NitrogenStatus.availableNitrogen < NitrogenStatus.totalNitrogendemand) 
		 {
			NitrogenStatus.leafNitrogenbuffer = 0.;
			NitrogenStatus.stemNitrogenbuffer = 0.;
			NitrogenStatus.tuberNitrogenbuffer = 0.;
			if (NitrogenStatus.actualleafNitrogenRatio > NitrogenStatus.minimumleafNitrogenRatio)
			{
				NitrogenStatus.leafNitrogenbuffer = NitrogenStatus.leafNitrogenAmount - this->get_leafMass() * NitrogenStatus.minimumleafNitrogenRatio;
				if (NitrogenStatus.leafNitrogenbuffer>0.04*NitrogenStatus.leafNitrogenAmount) NitrogenStatus.leafNitrogenbuffer = 0.04*NitrogenStatus.leafNitrogenAmount;
			}
			if (NitrogenStatus.actualstemNitrogenRatio > NitrogenStatus.minimumstemNitrogenRatio)
			{
				NitrogenStatus.stemNitrogenbuffer = NitrogenStatus.stemNitrogenAmount - this->get_stemMass() * NitrogenStatus.minimumstemNitrogenRatio;
				if (NitrogenStatus.stemNitrogenbuffer > 0.04*NitrogenStatus.stemNitrogenAmount) NitrogenStatus.stemNitrogenbuffer = 0.04*NitrogenStatus.stemNitrogenAmount;
			}
			if (NitrogenStatus.actualtuberNitrogenRatio > NitrogenStatus.minimumtuberNitrogenRatio)
			{
				NitrogenStatus.tuberNitrogenbuffer = NitrogenStatus.tuberNitrogenAmount - this->get_tuberMass() * NitrogenStatus.minimumtuberNitrogenRatio;
				if (NitrogenStatus.tuberNitrogenbuffer > 0.04*NitrogenStatus.tuberNitrogenAmount) NitrogenStatus.tuberNitrogenbuffer = 0.04* NitrogenStatus.tuberNitrogenAmount;
			}
			NitrogenStatus.availableNitrogen += NitrogenStatus.leafNitrogenbuffer+NitrogenStatus.stemNitrogenbuffer+NitrogenStatus.tuberNitrogenbuffer;
			NitrogenStatus.leafNitrogenAmount -= NitrogenStatus.leafNitrogenbuffer;
			NitrogenStatus.stemNitrogenAmount -= NitrogenStatus.stemNitrogenbuffer;
			NitrogenStatus.tuberNitrogenAmount -= NitrogenStatus.tuberNitrogenbuffer;
			NitrogenStatus.leafNitrogendemand = dblLeafgro * NitrogenStatus.optimumleafNitrogenRatio;
			NitrogenStatus.stemNitrogendemand = dblStemgro * NitrogenStatus.optimumstemNitrogenRatio;
			NitrogenStatus.tuberNitrogendemand = dblTubgro * NitrogenStatus.optimumtuberNitrogenRatio;
			NitrogenStatus.totalNitrogendemand = NitrogenStatus.leafNitrogendemand+NitrogenStatus.stemNitrogendemand+NitrogenStatus.rootNitrogendemand+NitrogenStatus.tuberNitrogendemand;
			if (NitrogenStatus.availableNitrogen > NitrogenStatus.totalNitrogendemand)
			{
				if (NitrogenStatus.totalNitrogendemand > 0) 
				{
					ratio = NitrogenStatus.availableNitrogen / NitrogenStatus.totalNitrogendemand;
					NitrogenStatus.leafNitrogendemand *= ratio;
					NitrogenStatus.stemNitrogendemand *= ratio;
					NitrogenStatus.tuberNitrogendemand *= ratio;
					NitrogenStatus.totalNitrogendemand = NitrogenStatus.leafNitrogendemand+NitrogenStatus.stemNitrogendemand+NitrogenStatus.rootNitrogendemand+NitrogenStatus.tuberNitrogendemand;
					//   ADJUST TNDEM FOR ROUNDING ERRORS
					if (NitrogenStatus.totalNitrogendemand > NitrogenStatus.availableNitrogen) NitrogenStatus.leafNitrogendemand -= (NitrogenStatus.totalNitrogendemand-NitrogenStatus.availableNitrogen);
				}
			}
			//   RECALCULATE N CONCENTRATIONS AFTER MOVING N OUT OF ORGANS
			NitrogenStatus.actualleafNitrogenRatio = 0.;
			if ((this->get_leafMass()+dblLeafgro) > 0) NitrogenStatus.actualleafNitrogenRatio = NitrogenStatus.leafNitrogenAmount / (this->get_leafMass() + dblLeafgro);
			NitrogenStatus.actualstemNitrogenRatio = 0.;
			if ((this->get_stemMass()+dblStemgro) > 0) NitrogenStatus.actualstemNitrogenRatio = NitrogenStatus.stemNitrogenAmount / (this->get_stemMass() + dblStemgro);
			NitrogenStatus.actualtuberNitrogenRatio = 0.;
			if ((this->get_tuberMass()+dblTubgro) > 0) NitrogenStatus.actualtuberNitrogenRatio = NitrogenStatus.tuberNitrogenAmount / (this->get_tuberMass() + dblTubgro);
		 }
		//   IF N IS LESS THAN OPTIMUM, SHIFT GROWTH FROM LEAVES TO STEMS 
		//      AND TUBERS ACCORDING N:C RATIO IN ARRAY TABLE
		//         IF (ANDEM.GT.AVAILN .OR. TANC.LT.TCNP) THEN
		 if (NitrogenStatus.totalNitrogendemand > NitrogenStatus.availableNitrogen)
		 {
			 ncratio = 0.;
			//C   CALCULATE CARBOHYDRATE AND NITROGEN AVAILABLE FOR TOP & TUBER GROWTH
			availC = (this->get_C_supply() - dblRootgro);//DHF - need to readjust this - added 24hour since we are comparing cumulative N availability to CHO fixed at current time-step.  THere's a disconnect with this thinking here...
			availN = (NitrogenStatus.availableNitrogen - NitrogenStatus.rootNitrogendemand);
			if (availC > 0) ncratio = availN/availC;
			ncratio = (NitrogenStatus.leafNitrogenAmount+NitrogenStatus.stemNitrogenAmount+NitrogenStatus.rootNitrogenAmount+NitrogenStatus.tuberNitrogenAmount)/(this->get_leafMass()+this->get_rootMass()+this->get_stemMass()+this->get_tuberMass());
			if (ncratio < 0.0275) ncratio = 0.0275;

			//   FIND ROW IN TABLE FOR CURRENT N:C RATIO
			int i = -1;
			while (tablelocation == false)
			{
				i = i + 1;
				if (ncratio < table[i][0]) tablelocation = true;
			}
			
			//dblLeafgro = availC*table[i-1][1]/100.;
			//dblStemgro = availC*table[i-1][2]/100.;
			//dblTubgro = availC*table[i-1][3]/100.;
			
            if (i < 22)
			{
				double interp = 0.;
				if (ncratio > table[i-1][0]/100.) interp = (ncratio-table[i-1][0])/(table[i][0]-table[i-1][0]);
				
				//dblLeafgro += availC*interp*(table[i][1]-table[i-1][1])/100.;
				//dblStemgro += availC*interp*(table[i][2]-table[i-1][2])/100.;
				//dblTubgro += availC*interp*(table[i][3]-table[i-1][3])/100.;
				
			}
			NitrogenStatus.leafNitrogendemand=dblLeafgro*NitrogenStatus.optimumleafNitrogenRatio;
			NitrogenStatus.stemNitrogendemand=dblStemgro*NitrogenStatus.optimumstemNitrogenRatio;
			NitrogenStatus.tuberNitrogendemand=dblTubgro*NitrogenStatus.optimumtuberNitrogenRatio;
			
			/*if (dblTubgro > this->get_tubers()->get_tubmax())
			{
				dblTubgro = this->get_tubers()->get_tubmax();
				NitrogenStatus.tuberNitrogendemand = dblTubgro*NitrogenStatus.optimumtuberNitrogenRatio;
				if (NitrogenStatus.availableNitrogen > (NitrogenStatus.leafNitrogendemand+NitrogenStatus.stemNitrogendemand+NitrogenStatus.rootNitrogendemand+NitrogenStatus.tuberNitrogendemand))
				{
		 			//  RATIO OF LEAF TO LEAF+STEM GROWTH
					ratio = (availC-dblTubgro)/availC;
					if (ratio < 0.5) ratio = 0.5;
					if (ratio > 0.75) ratio = 0.75;
					availN -= NitrogenStatus.totalNitrogendemand;
					//   GROLF/(GROLF+GROSTM) = RATIO : THUS AS PORTION OF TUBER GROWTH DECREASES,
					//      LEAF GROWTH IS FAVORED OVER STEM GROWTH
					//   GROLF*TCNP + GROSTM*SCNP = AVLN
					//   SOLVED FOR GROLF:
					//   GROLF = (AVLN-GROSTM*SCNP)/TCNP
					//   RATIO/(1-RATIO)*GROSTM*TCNP + GROSTM*SCNP = AVLN
					//   GROSTM*[TCNP*RATIO/(1-RATIO) + SCNP] = AVLN
					//   SOLVED FOR GROSTM:
					//   GROSTM = AVLN/[SCNP + TCNP*RATIO/(1-RATIO)]
					
					dblStemgro = availN/(NitrogenStatus.optimumstemNitrogenRatio+NitrogenStatus.optimumleafNitrogenRatio*ratio/(1.0-ratio));
					dblLeafgro = (availN-dblStemgro*NitrogenStatus.optimumstemNitrogenRatio)/NitrogenStatus.optimumleafNitrogenRatio;
					
					NitrogenStatus.leafNitrogendemand=dblLeafgro*NitrogenStatus.optimumleafNitrogenRatio;
					NitrogenStatus.stemNitrogendemand=dblStemgro*NitrogenStatus.optimumstemNitrogenRatio;
				}
			}
			*/
			NitrogenStatus.totalNitrogendemand = NitrogenStatus.leafNitrogendemand+NitrogenStatus.stemNitrogendemand+NitrogenStatus.rootNitrogendemand+NitrogenStatus.tuberNitrogendemand;
			if ((NitrogenStatus.totalNitrogendemand < NitrogenStatus.availableNitrogen)&&(NitrogenStatus.totalNitrogendemand>0))
			{
				ratio = (NitrogenStatus.availableNitrogen-NitrogenStatus.rootNitrogendemand)/(NitrogenStatus.totalNitrogendemand-NitrogenStatus.rootNitrogendemand);
				NitrogenStatus.leafNitrogendemand *= ratio;
				NitrogenStatus.stemNitrogendemand *= ratio;
				NitrogenStatus.tuberNitrogendemand *= ratio;
				NitrogenStatus.totalNitrogendemand = NitrogenStatus.leafNitrogendemand+NitrogenStatus.stemNitrogendemand+NitrogenStatus.rootNitrogendemand+NitrogenStatus.tuberNitrogendemand;
				//   CHECK FOR SMALL ROUNDING ERRORS AND ADJUST: 
				//   ANDEM SHOULD EQUAL AVAILN HERE!
				if (NitrogenStatus.totalNitrogendemand > NitrogenStatus.availableNitrogen) NitrogenStatus.leafNitrogendemand -= (NitrogenStatus.totalNitrogendemand-NitrogenStatus.availableNitrogen);
				if (NitrogenStatus.totalNitrogendemand > NitrogenStatus.availableNitrogen) NitrogenStatus.totalNitrogendemand = NitrogenStatus.availableNitrogen;
			}
			if (NitrogenStatus.totalNitrogendemand > NitrogenStatus.availableNitrogen)
			{
				//  At this point, if tndem & GROLF >0 
				//     then we have a small error from the TABLE values
				if (dblLeafgro > 0)
				{
					ratio = (NitrogenStatus.availableNitrogen - NitrogenStatus.rootNitrogendemand)/(NitrogenStatus.totalNitrogendemand-NitrogenStatus.rootNitrogendemand);
					dblLeafgro *= ratio;
					NitrogenStatus.leafNitrogendemand *= ratio;
					dblStemgro *= ratio;
					NitrogenStatus.stemNitrogendemand *= ratio;
					dblTubgro *= ratio;
					NitrogenStatus.tuberNitrogendemand *= ratio;
					NitrogenStatus.totalNitrogendemand = NitrogenStatus.leafNitrogendemand+NitrogenStatus.stemNitrogendemand+NitrogenStatus.rootNitrogendemand+NitrogenStatus.tuberNitrogendemand;
					//  If tndem<=0 then low N will cause N to be extracted from the leaves
					//   NITROGEN ABOVE TMNC CAN BE MOVED OUT OF LEAVES AT 4% /DAY (SEE ABOVE)
					//     SO IF STILL TANC > TMNC THEN JUST REDUCE TUBER GROWTH TO MATCH ALVN
				}
				else if (NitrogenStatus.actualleafNitrogenRatio > NitrogenStatus.minimumleafNitrogenRatio)
				{
					NitrogenStatus.tuberNitrogendemand=NitrogenStatus.availableNitrogen;
					if (NitrogenStatus.tuberNitrogendemand > 0)
					{
						dblTubgro = NitrogenStatus.tuberNitrogendemand / NitrogenStatus.optimumtuberNitrogenRatio;
					}
					else 
					{
						NitrogenStatus.tuberNitrogendemand = 0.;
						dblTubgro = 0.;
					}
				}
				//   TOO MUCH CARBON FOR TUBER GROWTH COMPARED TO AVAIL N, SO PULL 
				//       N & C OUT OF LEAVES UP TO 20% OF LEAF WEIGHT
				//    XC IS C LOSS FROM LEAVES; XN IS N LOSS FROM LEAVES
				//    1/2 OF XC WILL GO TO TUBERS, REST IS SENESCED
				//    AVLC IS C ALREADY AVAILABLE FOR TUBER GROWTH
				//    AVLN IS N ALREADY AVAILABLE FOR TUBER GROWTH
				//
				//   { XC/2 + AVLC } * TUBCNP   =   AVLN  +  XC * TMNC
				//    C FOR TUBER     CONVERTED               N FREED FROM
				//     GROWTH        TO N BASIS              LEAVES BY LOSS OF XC
				//
				//    SOLVE FOR XC:
				//      XC * TUBCNP/2 - XC * TMNC = AVLN - AVLC * TUBCNP
				//      XC = (AVLN - AVLC*TUBCNP)/(0.5*TUBCNP - TMNC)
				//    XN WILL BE XC*TMNC
				//
				else
				{
					xc = (availN - availC*NitrogenStatus.optimumtuberNitrogenRatio)/(0.5*NitrogenStatus.optimumtuberNitrogenRatio - NitrogenStatus.minimumleafNitrogenRatio);
					if ((availC + xc/2.2) > this->get_tubers()->get_tubmax()) xc = 2. * (this->get_tubers()->get_tubmax()-availC);
					if (xc > 0)
					{
						if (this->get_leafMass()>25)
						{
							if (xc > (0.2*this->get_leafMass())) xc = 0.2*this->get_leafMass();
						}
						else
						{
							if (xc > this->get_leafMass()) xc = this->get_leafMass();
						}
						xn = xc * NitrogenStatus.minimumleafNitrogenRatio;
						//   THIS IS GOING DOWN TOO FAST WITH .4  ->> TO .2
						dblLeafgro = -xc;
						NitrogenStatus.leafNitrogendemand = - xn;
						dblTubgro = xc/2. + availC;
						NitrogenStatus.tuberNitrogendemand = availN + xn;
						//    CHECK FOR NEGATIVE NITROGEN CONTENT OF TOP
						if ((NitrogenStatus.leafNitrogenAmount+NitrogenStatus.leafNitrogendemand) < 0)
						{
							NitrogenStatus.tuberNitrogendemand = NitrogenStatus.availableNitrogen - NitrogenStatus.rootNitrogendemand + NitrogenStatus.leafNitrogenAmount;
							NitrogenStatus.tuberNitrogendemand = -1. * NitrogenStatus.leafNitrogenAmount;
							//dblTubgro = NitrogenStatus.tuberNitrogendemand / NitrogenStatus.optimumtuberNitrogenRatio;
						}
					}
					//    IF XC=0, i.e. ONLY N IS NEEDED FOR MAXIMUM TUBER GROWTH, TAKE NEEDED
					//       N FROM LEAVES AND PUT INTO TUBDEM.  FREED C WILL GO INTO DEADLF
					else
					{
						xn = NitrogenStatus.totalNitrogendemand-NitrogenStatus.availableNitrogen;
						if (this->get_leafMass() > 25)
						{
							if (xn > (0.2*NitrogenStatus.leafNitrogenAmount)) xn = 0.2*NitrogenStatus.leafNitrogenAmount;
						}
						else
						{
							if (xn > NitrogenStatus.leafNitrogenAmount) xn = NitrogenStatus.leafNitrogenAmount;
						}
						xc = xn / NitrogenStatus.minimumleafNitrogenRatio;
						dblLeafgro = -xc;
						NitrogenStatus.leafNitrogendemand = - xn;
						NitrogenStatus.tuberNitrogendemand = availN + xn;
					}
				}
			}
		}
	} 
	NitrogenStatus.leafNitrogenAmount += NitrogenStatus.leafNitrogendemand;
	NitrogenStatus.stemNitrogenAmount += NitrogenStatus.stemNitrogendemand;
	NitrogenStatus.rootNitrogenAmount += NitrogenStatus.rootNitrogendemand;
	NitrogenStatus.tuberNitrogenAmount += NitrogenStatus.tuberNitrogendemand;
	NitrogenStatus.totalNitrogendemand=NitrogenStatus.leafNitrogendemand+NitrogenStatus.stemNitrogendemand+NitrogenStatus.rootNitrogendemand+NitrogenStatus.tuberNitrogendemand;
	NitrogenStatus.availableNitrogen -= (NitrogenStatus.leafNitrogendemand + NitrogenStatus.stemNitrogendemand + NitrogenStatus.rootNitrogendemand + NitrogenStatus.tuberNitrogendemand);
	if ((this->get_leafMass() + dblLeafgro) > 0)
		NitrogenStatus.actualleafNitrogenRatio = NitrogenStatus.leafNitrogenAmount/(this->get_leafMass() + dblLeafgro);
	else
		NitrogenStatus.actualleafNitrogenRatio = 0.; 
	if (NitrogenStatus.availableNitrogen < 0) NitrogenStatus.availableNitrogen = 0.;
	if (this->get_leafMass() > 0) //DHF - Add surplus N into leaf pool or stem pool
	{
		NitrogenStatus.leafNitrogenAmount += NitrogenStatus.availableNitrogen; //DHF - Adding available N into leaf class as surplus
		NitrogenStatus.availableNitrogen = 0.;
		NitrogenStatus.actualleafNitrogenRatio = NitrogenStatus.leafNitrogenAmount / (this->get_leafMass()+dblLeafgro);
	}
	else if (this->get_stemMass() > 0)
	{
		NitrogenStatus.stemNitrogenAmount += NitrogenStatus.availableNitrogen;
		NitrogenStatus.availableNitrogen = 0.;
		NitrogenStatus.actualstemNitrogenRatio = NitrogenStatus.stemNitrogenAmount / (this->get_stemMass()+dblStemgro);
	}
	this->set_TotalPlantNitrogen(NitrogenStatus.leafNitrogenAmount+NitrogenStatus.stemNitrogenAmount+NitrogenStatus.rootNitrogenAmount+NitrogenStatus.tuberNitrogenAmount);
}

void CPlant::potential_growth(const TInitInfo info, const TWeather & weather)
{
	//********************************************************************/
    // Calculate potential dry mass gain during time increment for each organ class
	// Based on potential leaf expansion rates determined from update_nodal method
	// At end of this routine, all C demand variables incorporate growth respiration cost
	//	-this cost ultimately gets accounted for when carbon is allocated in the C_allocation2, C_allocation3, and SetMass routines
	// Similar % of C allocation to stem, root, tubers as in SIMPOTATO
	// Root growth is adjusted based on water (to add) or nitrogen stress as in SIMPOTATO
	/* Notes
		-demand to root and shoots further modified in C_allocation1 as done in MAIZSIM based on alteration of C to roots in response to soil water potential
		-C partitioning to leaf,stem,tuber also modified based on N:C ratio if remobilized N not from leaf pool not sufficient in N_balance routines
		-C partitioning also modified if Cdemand exceeds Cpool 
	// * Note
		-need to add water component commented out 
	// *******************************************************************/
	double a,b,dblRstress = 1.;
	// Potential leaf growth
	dblPleafgro = dblPmainleafgro+dblPbasalleafgro+dblPapicalleafgro;
	dblnostressPleafgro = dblnostressPmainleafgro; 
	// Potential stem growth
	dblPstemgro = dblPleafgro*Rgleaf/Rgorgan / 3.; // demand as per SIMPOTATO
	dblPstemgro = dblPleafgro * Rgleaf/Rgorgan/2.; //modified as per experimental observations
	dblnostressPstemgro = dblnostressPleafgro * Rgleaf / Rgorgan / 2.;
	// Potential tuber growth
	if (tubers->isInitiated()) dblPtubgro = tubers->Growth(develop,Tdaylag,info)/Rgorgan;// *1)  calculate potential tuber growth and dry solids content as per SIMPOTATO
	else
		dblPtubgro = 0.;
	// Potential root growth
	//dblProotgro = 3 * develop->get_DTT17()/ info.plantDensity * (1-develop->get_CDD17()/35); // from SIMGUI - not realistic at all!
	if (tubers->isInitiated()) {
		dblProotgro = dblPleafgro*Rgleaf/Rgorgan*0.2 + dblPstemgro*0.2; // from SUBSTOR - SIMGUI code not robust with wide range of Temps
		dblnostressProotgro = dblnostressPleafgro*Rgleaf/Rgorgan*0.2 + dblnostressPstemgro*0.2; 
	}
	else {
		dblProotgro = dblPleafgro*Rgleaf/Rgorgan*0.5 + dblPstemgro*0.5; // tried to bump up root growth earlier in season
		dblnostressProotgro = dblnostressPleafgro*Rgleaf/Rgorgan*0.5 + dblnostressPstemgro*0.5; 
		//dblProotgro = dblPstemgro*0.4; //was 0.4, not enough early root growth
		//dblnostressProotgro = dblnostressPstemgro*0.4;
	}
	if (develop->get_istage() < 3 && NitrogenStatus.Nitrogendeficiencytwo < 1) { // water or N stress affects ala SIMPOTATO
		// Root growth can increase up to 25% under water or N stress, - W.G. Burton (1989), The Potato, p. 124.
		//DHF - increase in root growth due to water stress handled by 2dsoil
		//DHF - increase due to N stress handled here below but only consider if this is less than water stress
		//  That is, use the minimum / most limiting stress so we don't double penalize
		if (weather.pcrs < dblProotgro) {//pcrs is the amount of cho sent to roots at prior time-step to account for extra root growth needed for increased W uptake
			dblProotgro = dblProotgro / max(NitrogenStatus.Nitrogendeficiencytwo, 0.8);
			dblnostressProotgro = dblnostressProotgro / max(NitrogenStatus.Nitrogendeficiencytwo, 0.8);
		}

		if (info.Water_stress_simulation_type ==1){
			a = min(NitrogenStatus.Nitrogendeficiencytwo, weather.swdf2);
			b = min(a,1.);
			dblRstress = max(b,0.8);
			dblProotgro = dblProotgro / dblRstress;
		}	
		
	}
	if ((this->get_develop()->get_istage() == 3) && (dblProotgro > (assimilate*CH2O_MW/CO2_MW + C_deadpool)*0.1)) 
	{
		dblProotgro = (assimilate*CH2O_MW/CO2_MW + C_deadpool) * 0.1; //limit root growth to 10% of assimilate supply after dominant tuber growth begins
	}
	if ((this->get_develop()->get_istage() == 3) && (dblnostressProotgro > (assimilate*CH2O_MW/CO2_MW + C_deadpool)*0.1)) 
	{
		dblnostressProotgro = (assimilate*CH2O_MW/CO2_MW + C_deadpool) * 0.1; //limit root growth to 10% of assimilate supply after dominant tuber growth begins
	}
	if (dblProotgro < 0) dblProotgro = 0.;
	if (dblnostressProotgro < 0) dblnostressProotgro = 0.;
	dblnostressPtubgro = dblPtubgro;
	dblLeafgro = dblPleafgro;
	dblStemgro = dblPstemgro;
	dblRootgro = dblProotgro;
	dblTubgro = dblPtubgro;
}

void CPlant::calcMaxPoolSize()
{
	/*****************************************/
	/* Idea is to estimate max size for Cpool based on ability of leaves to store excess C
	/* As Cpool increases, sometype of down-regulation of photosynthetic activity occurs
	/* No penalty, however, until leaf area exceeds 400 cm2
	/* Only used for individual leaf approach
	/******************************************/
	double C_PoolRoom1 = 0.; //how much more C can leaves hold before being saturated - based on sla
	double C_PoolRoom2 = 0.; //as above, based on haulm mass
	double C_PoolRatio = 0.; //ratio tracks how much room is left in leaves to hold more C - at 1.0, no more reserve can be stored
	double temp = 0.;
	if ((initInfo.bigleaf == 0) && (inodeNumber > 0))
	{
		for (int i=0; i<inodeNumber; i++) //get total room left in leaf for C based on current area, SLA, and minimum possible SLA
		{
			C_pool_room += nodalUnit[i].get_leaf()->get_greenArea()*((1./nodalUnit[i].get_leaf()->get_SLAmin())-(1./nodalUnit[i].get_leaf()->get_SLAleaf()));
		}
	}
	if ((initInfo.bigleaf == 1) && (inodeNumber > 0))
	{
		C_PoolRoom1 = nodalUnit[0].get_leaf()->get_greenArea()*((1./nodalUnit[0].get_leaf()->get_SLAmin())-(1./nodalUnit[0].get_leaf()->get_SLAleaf()));
		//alternative approach
		C_PoolRoom2 = (this->get_leafMass()+this->get_stemMass())*0.1;
		C_pool_room = min(C_PoolRoom1, C_PoolRoom2);

	}
	C_pool_room = C_PoolRoom1;
	C_PoolRatio = C_pool / C_pool_room;
	//temp = exp(-1.0*C_PoolRatio);
	temp = (1. - C_PoolRatio);
	//if (assimilate_old > 0){
	//	C_PoolRatio = C_pool / (assimilate_old*CH2O_MW/CO2_MW);

	PhotosyntheticFeedback = min(1.,temp);
	if (PhotosyntheticFeedback <0.1) PhotosyntheticFeedback = 0.1;
	//	temp = (4-C_PoolRatio)/4.;
//		if (temp <0) temp = 0;
//		if (temp > 1) temp = 1;
//		PhotosyntheticFeedback = temp;
//	}
//	else
//		PhotosyntheticFeedback = 1;//only at night?

	EffectiveCanopySLA = this->get_greenLeafArea()/(this->get_leafMass()+C_pool);
	if (this->get_leafArea() < 400) PhotosyntheticFeedback = 1.; //no penalty until plant is actively growing in exp. phase	
	//PhotosyntheticFeedback = 1;


}


void CPlant::reset()
{
	/******************************/
	/* All temporary variables are reset to 0 at the beginning of each time-step
	/* Variables include those that track potential and actual C and area growth demand and gas exchagne rates
	/* Also include architectural information
	/* ***************************/
	dblPleafgro = dblPmainleafgro = dblPbasalleafgro = dblPapicalleafgro = dblPtubgro = dblPstemgro = dblProotgro=0.;
	dblTubgro = dblLeafgro = dblMainleafgro = dblBasalleafgro = dblApicalleafgro = dblStemgro = dblRootgro = 0.;
	dblPyoungleafgro = dblPoldleafgro = dblYoungleafgro = dblOldleafgro = 0.;
	dblnostressPleafgro = dblnostressPstemgro = dblnostressProotgro = dblnostressPtubgro = 0.;
	dblnostressLeafgro = dblnostressStemgro = dblnostressRootgro = dblnostressTubgro = 0.;
	iyoungnodeNumber = ioldnodeNumber = igreennodeNumber = 0;
	mainleafArea = basalleafArea = apicalleafArea = 0.;
	mainstemMass = basalstemMass = apicalstemMass = 0.;
	mainleafMass = basalleafMass = apicalleafMass = deadMass = 0.;
	mainleafNo = basalleafNo = apicalleafNo = basalstemNo = apicalstemNo = 0;
	assimilate = 0.; maintRespiration = 0.; reference = 0.;
	photosynthesis_gross = 0.; photosynthesis_net = 0.; transpiration = 0.; instantResp = 0.;
	shootPart_old = shootPart;
	rootPart_old = rootPart;
	tubPart_old = tubPart;
	C_demand_old = C_demand;
	shootPart = rootPart = tubPart = 0.;
	swdf1 = psistress_gs_factor = 1.;

	this->leafageEffect = 0.;
	this->set_N_hourlyseedused(0.);
	this->set_N_hourlytranslocated(0.);
	this->set_N_hourlygrowthdemand(0.);
	N_dead = 0.; newleafNitrogen = 0.; newstemNitrogen = 0.;
	C_pool_used = C_seed_used = 0.;
	C_deadpool = 0.;
}


void CPlant::set_mass()
//Method sums up mass for whole organ classes and plant N status
/* These are running totals over course of season */
{
	if (initInfo.bigleaf == 0)
	{
		for (int i = 0; i<inodeNumber; i++)
		{
			if (nodalUnit[i].get_leaf()->isTerminated()) deadMass += nodalUnit[i].get_leaf()->get_deadmass();
			// obtain values for stem and leaf branch types
			if (nodalUnit[i].get_type() == 0){
				mainstemMass += nodalUnit[i].get_stem()->get_drymass();
				mainleafMass += nodalUnit[i].get_leaf()->get_drymass();
				mainleafArea += nodalUnit[i].get_leaf()->get_greenArea();
				mainleafNo += 1;
			}
			else if (nodalUnit[i].get_type() == 1){
				basalstemMass += nodalUnit[i].get_stem()->get_drymass();
				basalleafMass += nodalUnit[i].get_leaf()->get_drymass();
				basalleafArea += nodalUnit[i].get_leaf()->get_greenArea();
				basalleafNo += 1;
				if (nodalUnit[i].get_node()==0 && nodalUnit[i].isInitiated()) basalstemNo  +=1; 
			}
			else if (nodalUnit[i].get_type() == 2){
				apicalstemMass += nodalUnit[i].get_stem()->get_drymass();
				apicalleafMass += nodalUnit[i].get_leaf()->get_drymass();
				apicalleafArea += nodalUnit[i].get_leaf()->get_greenArea();
				apicalleafNo += 1;
				if (nodalUnit[i].get_node()==0) apicalstemNo  +=1; 
			}
		}
	}
	else
	{
		if (nodalUnit!=NULL) 
		{
			mainstemMass = nodalUnit[0].get_stem()->get_drymass();
			deadMass = nodalUnit[0].get_leaf()->get_deadmass();
			mainleafMass = nodalUnit[0].get_leaf()->get_drymass();
			mainleafArea = nodalUnit[0].get_leaf()->get_greenArea();
			mainleafNo = basalleafNo = apicalleafNo = 0;
			basalstemMass = basalleafMass = basalleafArea = 0.;
			apicalstemMass = apicalleafMass = apicalleafArea = 0.;
		}
		else
		{
			NitrogenStatus.HourlyNitrogenDemand = 0.;
			NitrogenStatus.CumulativeNitrogenDemand = 0.;
		}
	}

	stemMass = mainstemMass+basalstemMass+apicalstemMass; 
	leafMass = mainleafMass+basalleafMass+apicalleafMass;// + C_pool*(C_MW/CH2O_MW)/C_content;
	rootMass = this->get_roots()->get_drymass();
	tuberMass = this->get_tubers()->get_drymass();
	C_deadpool_season += C_deadpool;
	//totalMass = stemMass + leafMass + deadMass + rootMass + tuberMass + C_reserve + C_pool + C_deadpool - C_seed;
	totalMass = stemMass + leafMass + deadMass + rootMass + tuberMass + C_pool + C_deadpool + C_pool_root;
	// N status
	//NitrogenStatus.leafNitrogenAmount -= newleafNitrogen; //remove N cost of any new nodal units from appropriate N pool
	//NitrogenStatus.stemNitrogenAmount -= newstemNitrogen;
	if (leafMass > 0) NitrogenStatus.actualleafNitrogenRatio = NitrogenStatus.leafNitrogenAmount / leafMass;
	else
		NitrogenStatus.actualleafNitrogenRatio = 0.;
	if (stemMass > 0) NitrogenStatus.actualstemNitrogenRatio = NitrogenStatus.stemNitrogenAmount / stemMass;
	else
		NitrogenStatus.actualstemNitrogenRatio = 0.;
	if (rootMass > 0) NitrogenStatus.actualrootNitrogenRatio = NitrogenStatus.rootNitrogenAmount / rootMass;
	else
		NitrogenStatus.actualrootNitrogenRatio = 0.;
	if (tuberMass > 0) NitrogenStatus.actualtuberNitrogenRatio = NitrogenStatus.tuberNitrogenAmount / tuberMass;
	else
		NitrogenStatus.actualtuberNitrogenRatio = 0.;
	NCRatioLeaves = NitrogenStatus.actualleafNitrogenRatio;
	NCRatioStems = NitrogenStatus.actualstemNitrogenRatio;
	NCRatioRoots = NitrogenStatus.actualrootNitrogenRatio;
	NCRatioTubers = NitrogenStatus.actualtuberNitrogenRatio;

	double leafNdemand = (NitrogenStatus.optimumleafNitrogenRatio - NitrogenStatus.actualleafNitrogenRatio) * leafMass;
	double stemNdemand = (NitrogenStatus.optimumstemNitrogenRatio - NitrogenStatus.actualstemNitrogenRatio) * stemMass;
	double rootNdemand = (NitrogenStatus.optimumrootNitrogenRatio - NitrogenStatus.actualrootNitrogenRatio) * rootMass;
	double tuberNdemand = (NitrogenStatus.optimumtuberNitrogenRatio - NitrogenStatus.actualtuberNitrogenRatio) * tuberMass;
	TotalNitrogen = NitrogenStatus.leafNitrogenAmount + NitrogenStatus.stemNitrogenAmount + NitrogenStatus.rootNitrogenAmount + NitrogenStatus.tuberNitrogenAmount;
	NCRatioTotalPlant = TotalNitrogen / totalMass;
	//totalDeadNitrogen = N_dead;	//already set in update_nodal method
	totalReserveNitrogen = NitrogenStatus.availableNitrogen;
	//NitrogenStatus.HourlyNitrogenDemand = __max(0,(leafNdemand + stemNdemand + rootNdemand + tuberNdemand)) + this->get_HourlyNitrogenSeedUsed() + this->get_HourlyNitrogenTranslocated();
	//NitrogenStatus.HourlyNitrogenDemand = leafNdemand + stemNdemand + rootNdemand + tuberNdemand;
	NitrogenStatus.HourlyNitrogenDemand = __max(0.,this->get_HourlyNitrogenGrowthDemand());
	NitrogenStatus.HourlyNitrogenDemand = __max(0.,(leafNdemand + stemNdemand + rootNdemand + tuberNdemand));
	NitrogenStatus.CumulativeNitrogenDemand += NitrogenStatus.HourlyNitrogenDemand;
	
	return;
}

void CPlant::update(int iCur, const TWeather & weather, double lwpd, const TInitInfo info, double dayinc)
{
	/***********************/
	/* This subroutine iterates through all aspects of plant growth and development
	/* Called by Controller:: Run method during each crop.dll call from 2dSOIL
	/**************************************/
	reset(); //reset variables at start of each time step
	/******************************************************/
	/* FIX ROUTINES SO N N Tracking or Stress - for PAPER*/
	/*  To be replaced with routine that will set all NitrogenStatus components to non-stressed level and keep them there*/

	NitrogenStatus.Nitrogendeficiencyone = 1.;
	NitrogenStatus.Nitrogendeficiencytwo = 1.;
	NitrogenStatus.Nitrogendeficiencythree = 1.;
	NitrogenStatus.Nitrogenstressfactor = 1.;
	this->NitrogenStressFactor = 1.;
	/******************************************************/
	dailyave(iCur, weather, info, dayinc); //track daily temperature values
	calcLeafArea(info);
	develop->update(iCur, weather, info, NitrogenStatus, get_greenLeafArea(), greenLAI, maxLAI, dayinc, Tdaymax, Tdaymin, Tdayave, Tdaylag, Sradave, Photoave, 0, C_pool); //update for particular stage C_reserve+C_poool
	update_nodes(iCur, weather, info, lwpd); //update nodal unit status, get senescence, estimate potential leaf expansion
	update_tubers(Tubfraction); //use lookup table to estimate tuber sink strength on photosynthesis
	update_roots();
	if (!develop->Matured() && develop->Emerged())
	{	
		
		nitrogen_stress(info); //calculate N stress factors on growth and development
		calcMaxPoolSize();//calculate max permittable size of CHO reserve and photosynthetic feedback factors, 2009 Paper
		calcGasExchange(weather, info); //calculate CHO fixation prior to potential expansion and limited by N stress
		potential_growth(info, weather); //calculate potential dry mass gain for current time increment		
		calcMaintRespiration(weather,info); //adjust for respiration of existing biomass
		C_allocation1(iCur, weather, info, lwpd); //allocation CHO pool and compute actual expansion
		C_allocation2(iCur, weather, info, lwpd); //allocate CHO and N within leaves in the canopy
		C_allocation3(iCur, weather, info); //allocate CHO and N to remaining organs
		initiate(iCur, weather); //add new organs (right now, just nodal units since tuber and roots already present) where necessary
	}
	set_mass();

}


void CPlant::update_nodes(int iCur, const TWeather& weather, const TInitInfo info, double lwpd)
{
	/**************************************/
	/* Updates canopy leaf development by simulating aging of existing leaves, potential growth and appearance of new nodalunits
	//Note: bigleaf = 0 is for individual organ leaf model, = 1 is for bigleaf, SIMGUI based leaf model
	//		-for SIMGUI approach, all leaf area and mass information stored in NodalUnit[0] only
	/* Accesses methods from nodalunit.cpp and leaf.cpp classes for aging, expansion, senescence
	/* For individual organ model (bigleaf = 0):
			Sort through each individual nodalunit and get aging and potential growth estimates
	/*		Determines if flowering occurs; checks leaf appearance and branch appearance rates along with Creserve to determine
	/*		if new nodalunits and branches can form
	/*		A new leaf can form on a given lateral branch at any time as long as sufficient N and C are available
	/*		Lateral branch appearance rate only influenced by CHO, N, and thermal time - no apical dominance as in POTATO
	/*		If leaves have been senesced, adds deadpool of that leaf to C_reserve - deadpool then reset to 0 to prevent C balancing issues
	/*		Also computes overall canopy level leafageeffect that weights individual leaf age on p.s. response
	/*		Also includes whole canopy shading effect to increase senescence rate if canopy too large as in SIMPOTATO
	/* NOTES:
			-Need to incorporate N effect on new growth of leaves / stems (make sure not dupliced in development-UPDATE method
				when computing leaf and branch appearance rates

	/* For big leaf model:
			All  leaf information is stored in nodalunit element 0, thus no sorting required
			Flowering, branch initiation not simulated
			Method simply estimates potential expansion of entire canopy
			Sensence assumes 50% CHO recovered, 75% of N as per SIMPOTATO
	/**************************************/

	double templeafdead = 0.;
	double tempstemdead = 0.;
	static double dSenescedarea; // SIMPOTATO variable - if > 0, indicates shading effect forcing additional senescence
	if (initInfo.bigleaf == 0)
	{
		if (develop->Vegetative()==true && inodeNumber > 0) //update leaf status
		{
			NitrogenStatus.leafNitrogenAmount = NitrogenStatus.stemNitrogenAmount = 0;
			Bnewnode = false;
			for (int i=0; i < inodeNumber; i++)
			{
				if (nodalUnit[i].isInitiated())
				{
					igreennodeNumber += 1;
					nodalUnit[i].update(iCur, develop, info, weather, NitrogenStatus, Tdaylag, Sradave, greenLAI, 0, 0., lwpd, swdf1, C_pool_room); //use organ update? - 0 indicates compute potential expansion
					//want to estimate potential expansion and and CHO demand here
					if (nodalUnit[i].get_type()==0) dblPmainleafgro += nodalUnit[i].get_leaf()->get_potentialDWinc()/Rgleaf;
					if (nodalUnit[i].get_type()==1) dblPbasalleafgro +=nodalUnit[i].get_leaf()->get_potentialDWinc()/Rgleaf;
					if (nodalUnit[i].get_type()==2) dblPapicalleafgro +=nodalUnit[i].get_leaf()->get_potentialDWinc()/Rgleaf;
					if (nodalUnit[i].get_LAR() > 1) Bnewnode = true;
					if (nodalUnit[i].get_leaf()->isYoung())
					{
						dblPyoungleafgro += nodalUnit[i].get_leaf()->get_potentialDWinc()/Rgleaf;
						iyoungnodeNumber += 1;
					}
					if (nodalUnit[i].get_leaf()->isGrowing() && !(nodalUnit[i].get_leaf()->isYoung()))
					{
						dblPoldleafgro += nodalUnit[i].get_leaf()->get_potentialDWinc()/Rgleaf;
						ioldnodeNumber += 1;
					}
					this->leafageEffect += nodalUnit[i].get_leaf()->get_ageEffect()*(nodalUnit[i].get_leaf()->get_greenArea()/this->get_greenLeafArea());
				}			
				if (nodalUnit[i].get_leaf()->isTerminated())
				{
					C_deadpool += nodalUnit[i].get_leaf()->get_deadCHOpool();//note that this routine also sets leaf deadpool to 0
					N_dead += nodalUnit[i].get_leaf()->get_deadNmass(); //track N lost to plant
					double Ntemp = nodalUnit[i].get_leaf()->get_currentNitrogenAmount();
					NitrogenStatus.availableNitrogen += Ntemp; //assumes N from senesced leaf instantaneously available
					igreennodeNumber -= 1;
				}
				NitrogenStatus.leafNitrogenAmount += nodalUnit[i].get_leaf()->get_currentNitrogenAmount();
				NitrogenStatus.stemNitrogenAmount += nodalUnit[i].get_stem()->get_currentNitrogenAmount();
			}
		}
		// Additional senescence routine required to simulate shading effect of too large canopy when LAI greater than 4.0
		// Derived from SIMPOTATO
		dSenescedarea += this->get_greenLeafArea() * (1-(1.-0.008*(this->get_greenLAI()-4.)))/(24.*60./info.timeStep);
		if (dSenescedarea > 0) //force senescence of older leaves, 1 at a time, until senescence requirement is met and/or exceeded.
		{
			for (int i = 0; i < inodeNumber; i++)
			{
				if (!nodalUnit[i].get_leaf()->isTerminated()&&nodalUnit[i].isInitiated()&&dSenescedarea>0)
				{
					//force senescence of older leaves, 1 at a time, until senesced area requirement is met
					//if (dSenescedarea > nodalUnit[i].get_leaf()->get_greenArea())
					//{
						dSenescedarea = dSenescedarea - nodalUnit[i].get_leaf()->get_greenArea();
						nodalUnit[i].get_leaf()->senescence2(info);
						C_deadpool += nodalUnit[i].get_leaf()->get_deadCHOpool();//note that this routine also sets leaf deadpool to 0
						N_dead += nodalUnit[i].get_leaf()->get_deadNmass();
						double temp = nodalUnit[i].get_leaf()->get_deadNpool();
						NitrogenStatus.leafNitrogenAmount -= temp;
						NitrogenStatus.availableNitrogen += temp;
						igreennodeNumber -= 1;

					//}
				}
			}
		}

	}
	else if (initInfo.bigleaf == 1) //big leaf expansion starts here - majority of routine comes from SIMPOTATO//
	{
		if (!nodalUnit=='\0')
		{
			nodalUnit[0].update(iCur, develop, info, weather, NitrogenStatus, Tdaylag, Sradave, greenLAI, 0, 0., lwpd,swdf1,C_pool_room);
			
			dblPmainleafgro = nodalUnit[0].get_leaf()->get_potentialDWinc()/Rgleaf; //for big leaf model, use mainleaf demand only to simplify C allocation
			dblnostressPmainleafgro = nodalUnit[0].get_leaf()->get_nostresspotentialDWinc()/Rgleaf; // as above, except non0water stress ptoential
			
			N_dead += nodalUnit[0].get_leaf()->get_deadNmass(); //track N lost to senescence from plant
			N_dead += nodalUnit[0].get_stem()->get_deadNmass(); // as above
			NitrogenStatus.availableNitrogen += nodalUnit[0].get_leaf()->get_deadNpool();
			NitrogenStatus.leafNitrogenAmount = nodalUnit[0].get_leaf()->get_currentNitrogenAmount();
			NitrogenStatus.availableNitrogen += nodalUnit[0].get_stem()->get_deadNpool();
			NitrogenStatus.stemNitrogenAmount = nodalUnit[0].get_stem()->get_currentNitrogenAmount();
			//nodalUnit[0].get_leaf()->senescence2(info); //check for shaded canopy leaf senescence - already done in alternative senescence routine in leaf.cpp class
			C_deadpool += nodalUnit[0].get_leaf()->get_deadCHOpool();//note that this routine gets 50% senesced CHO from both routines and also sets leaf deadpool to 0
			//NitrogenStatus.availableNitrogen += nodalUnit[0].get_leaf()->get_deadNpool();
			//NitrogenStatus.leafNitrogenAmount = nodalUnit[0].get_leaf()->get_currentNitrogenAmount();
			C_deadpool += nodalUnit[0].get_stem()->get_deadCHOpool();
			//nodalUnit[0].get_leaf()->set_currentNitrogenAmount((nodalUnit[0].get_leaf()->get_currentNitrogenAmount()-(.25*templeafdead*NitrogenStatus.actualleafNitrogenRatio)));
			//this->totalDeadNitrogen += 0.25*templeafdead*NitrogenStatus.actualleafNitrogenRatio;
			this->totalDeadNitrogen += nodalUnit[0].get_leaf()->get_deadNmass();
			//nodalUnit[0].get_stem()->set_currentNitrogenAmount((nodalUnit[0].get_stem()->get_currentNitrogenAmount()-(.25*tempstemdead*NitrogenStatus.actualstemNitrogenRatio)));
			//this->totalDeadNitrogen += 0.25*tempstemdead*NitrogenStatus.actualstemNitrogenRatio;
			this->totalDeadNitrogen += nodalUnit[0].get_stem()->get_deadNmass();
			this->leafageEffect = 1.;
			//For big leaf model, N info is tracked in whole plant information in Nitrogen_Stress, Nitrogen_Balance routines
			NitrogenStatus.leafNitrogenAmount = nodalUnit[0].get_leaf()->get_currentNitrogenAmount();
			NitrogenStatus.stemNitrogenAmount = nodalUnit[0].get_stem()->get_currentNitrogenAmount();
		}
		else
		{
			dblnostressPmainleafgro = 0.; dblPmainleafgro = 0.; C_deadpool += 0.;
		}
	}
}

void CPlant::update_tubers(double Tubfraction)
{
	/******************/
	/* Method determines feedback effect of tuber sinkstrength on leaf gas exchange
	/* From modified Ng and Loomis*/
	/******************/
	if (this->get_tubers()->isInitiated())
	{
		Tubmod = this->get_tubers()->get_feedback(Tubfraction);
		//For Big leaf model, N status tracked only in whole plant routines, Nitrogen_Stress, Nitrogen_Balance
		NitrogenStatus.tuberNitrogenAmount = this->get_tubers()->get_currentNitrogenAmount();
	}
}

void CPlant::update_roots()
{
	/**************/
	/* Method updates N status of roots*/
	/******************************/
	//if ((this->get_roots()!=NULL)&&(this->get_rootMass()>0))
	//{
		//For Big leaf model, N status tracked only in whole plant routines, Nitrogen_Stress, Nitrogen_Balance
		NitrogenStatus.rootNitrogenAmount = this->get_roots()->get_currentNitrogenAmount();
	//}
}

void CPlant::C_allocation1(int iCur, const TWeather & w, const TInitInfo info, double lwpd)
{		

/********************************************************/
/* Routine primarily consists of modified SIMPOTATO code with some additions from MAIZESIM as noted in the code
/* Routine compares demand for C of different organ classes with C available for growth
/* Note that leaf demand was determined during update_nodes while demand for other organs determined during potential_growth() methods.
/* Modifies the demand if supply is insufficient based on fixed percentages and developmental stages
/* At end of routine, exact amounts of C to be distributed to each organ class is specified
/* An additional N balancing method is called to re-adjust CHO partitioning based on degree of N-stress
/*   -routine comes from SIMGUI and varies for pre or post tuber initiation
/* Needs:
	-For single organ model (bigleaf = 0), must incorporate N effects on partitioning - right now, assuming bigleaf SIMGUI approach, not sure how this will work with single leaf
/* *******************************************************/
	double tmprEffect = 1.;
	double temp=0.;
	double grofac = 0.2/(60./initInfo.timeStep); // translocation limitation and lag, assume it takes 1 hours to complete, 0.2=5 hrs
	grofac = 1;//try 30 minuter interval, too much C building up in plant which is typically C starved!?
	// this is where source/sink (supply/demand) valve can come in to play
   // 0.2 is value for hourly interval, Grant (1989)
	double excess = 0.; 
	double C_avail = 0., N_avail=0.;
	double shootPart_real = 0.; double rootPart_real = 0.; double shootPart_virtual = 0.;
	//double dt = info.timeStep/(24*60);

	//************************************************
	// Step 1: Determine available C supply
	//************************************************
	//a) 10% of seed reserve is available for growth per 24 hour period
	//b) 10% of C reserve is available for growth per 24 hour period
	//c) C_pool is mobile CHO used for growth and maint, residing in haulm; however, only portion can be accessed each time-step due to translocation lags
	//d) C_supply is actual fraction of C_pool available for growth and maint in curren time-step
	//e) C_supply_nostress is a virtual variable used to provide a minimum amount of CHO to root growth when leaf expansion ceases due to water-stress

	C_pool += assimilate*CH2O_MW/CO2_MW + C_deadpool;
	C_supply += __max(C_pool*tmprEffect*grofac,0.);
	C_supply_nostress = C_supply;
	C_pool -= __max(C_pool*tmprEffect*grofac,0.);
	C_pool_nostress = C_pool;

	//f) additional proportion of C_pool becomese available for below ground growth during water stress if CHO starts to build up in haulm
	//    don't have specific experimental data on labile C in plant except partitioning results at end-of-season confirm shift in tuber growth for water stressed plants
	//		-needs a bit of work still-
	if (this->PhotosyntheticFeedback < 0.95) {
		//to start, amount of existing CHO available for root growth is 10% of this C_room
		//double Croom = C_pool / (1 - this->PhotosyntheticFeedback);
		//double Cexcess = C_pool - 0.1*(C_pool/(1-this->PhotosyntheticFeedback)); //that is, anything over 10% of C_room is available for root or tuber
		//double Cadd = Cexcess;
		//C_supply_nostress = C_supply + Cadd * grofac;
		//C_pool_nostress -= Cadd*grofac;
		if (((develop->get_istage() == 2) || (develop->get_istage() == 3)) && (lwpd/10. <-0.05)) {//if tubers initiated and predawn lwp is < -0.1 MPa, allow excess CHO from pool to go to tubers
			C_demand = dblPleafgro+dblPstemgro+dblProotgro+dblPtubgro+maintRespiration;
			if (C_demand < C_supply)//don't increase allocation if supply is not in excess to support it
			{
				//double Ctub = (C_pool - Cadd*grofac) * grofac; //too much to add in one time-step!
				double Cadd = (C_supply - C_demand);
				double Ctub = 0.2*dblPtubgro;
				if (Ctub > Cadd) 
				{
					Ctub = Cadd;
				}

				dblPtubgro += Ctub; 
				//C_supply += Ctub;
				//C_supply_nostress += Ctub;
				//C_pool -= Ctub;
				//C_pool_nostress -= Ctub;
			}
			
		}
		if (C_pool_root > 0 && ((develop->get_istage() == 2) || (develop->get_istage() == 3))) { //add C_pool_root to tuber growth if it's not all used up
			//C_pool_root is not included in C_pool, so separate accounting like this will have no mass balance problems
			dblPtubgro += C_pool_root;
			C_supply += C_pool_root;
			C_supply_nostress += C_pool_root;
			C_pool_root = 0.;
		}
		dblnostressPtubgro = dblPtubgro;
	}
	
	/****************************************/
	/* Step 2: Determine available C_supply to meet C_demand*/
	/****************************************/
	// implicit in the organ growth variables is growth respiration costs which are accounted for once CHO is actually partitioned*/
	C_demand = dblPleafgro+dblPstemgro+dblProotgro+dblPtubgro+maintRespiration;
	C_demand_nostress = dblnostressPleafgro + dblnostressPstemgro + dblnostressProotgro+ dblnostressPtubgro + maintRespiration; //C demand if no water stress involved

	if (this->get_leafArea() < 400) //assume initial maintRespiration costs are met with tuber seedpiece
	{
		if (C_seed > maintRespiration)
		{
			C_supply += maintRespiration;
			C_seed -=maintRespiration;
			C_seed_used += maintRespiration;
		}
		else
		{
			C_supply += C_seed;
			C_seed_used += maintRespiration;
			C_seed = 0;
		}
	}

	if (C_supply < C_demand && this->get_leafArea() < 400) //not enough soluble C to meet demand, add C_reserve (C_seed) to C_supply when leaf area < 400cm2 (Ng and Loomis)
	{
		//as in SIMGUI< limit this amount to 1.5 g total each 24h
		if (C_seed > 1.5) {
			C_avail = 1.5/((24.*60.)/info.timeStep);
			C_avail = grofac*1.5/((24.*60.)/info.timeStep);//induce lag on storage reserves
			//C_avail = 1;
		}
		else {
			C_avail = C_seed / ((24.*60.)/info.timeStep);
			C_avail = grofac*C_seed / ((24.*60.)/info.timeStep);//induce lag on storage reserves
			//C_avail = C_seed;
		}
		if ((C_supply + C_avail) > C_demand) C_avail = C_demand - C_supply; //limit C and N from reserve to only what's needed
		C_supply += C_avail;
		C_seed_used += __max(C_avail,0.);
        C_seed -= __max(C_avail,0.);
	}
	if (C_supply > C_demand) //put excess C pack into pool
	{
		C_pool += (C_supply-C_demand);
		C_supply = C_demand;
	}
	if (C_supply_nostress > C_demand_nostress) //put excess C back into virtual pool
	{
		C_pool_nostress += (C_supply_nostress - C_demand_nostress);
		C_supply_nostress = C_demand_nostress;
	}
	if (C_supply < maintRespiration)
	{
		//add additional reserve to satisfy Mresp?
		temp = maintRespiration-C_supply;
		C_supply += temp;
		C_pool -= __max(temp,0.);
		if (C_pool < 0.) {
			C_pool = 0.;
			if (C_supply < 0) C_supply = 0.;
		}
	}
	C_pool_used = C_supply;
	C_supply -= maintRespiration;
	C_demand -= maintRespiration;
	C_supply_nostress -= maintRespiration;
	C_demand_nostress -= maintRespiration;
	C_new_organ = C_supply - C_demand;
	/*Step 3: Allocate C_supply*/
	/**************************************************/
	// * 1) Vegetative growth partitioning
	if (develop->get_istage()== 1){
        //* Does demand exceed or equal available assimilate
		// case a: excess assimilate
		if (C_demand <= C_supply) //supply meets demand
		{
			dblStemgro = dblPstemgro;
			dblRootgro = dblProotgro;
			//dblMainleafgro = dblPmainleafgro;
			//dblApicalleafgro = dblPapicalleafgro;
			//dblBasalleafgro = dblPbasalleafgro;
			//dblLeafgro = dblMainleafgro + dblApicalleafgro + dblBasalleafgro + maintRespleaf;
			dblLeafgro = dblPleafgro;
			//dblYoungleafgro = dblPyoungleafgro;
			//dblOldleafgro = dblPoldleafgro;
		}
		else
		{
			//* reduce demand equally for all organs ala SimPotato (or alter demand so that at least maintResp is satisifed)
			dblLeafgro = dblPleafgro * C_supply / (dblPleafgro+dblPstemgro+dblProotgro);
			dblStemgro = dblPstemgro * C_supply / (dblPleafgro+dblPstemgro+dblProotgro);
			dblRootgro = dblProotgro * C_supply / (dblPleafgro+dblPstemgro+dblProotgro);
		}
		// case a: excess assimilate for theoretical nonwater stress calculation
		if (C_demand_nostress <= C_supply_nostress)
		{
			dblnostressStemgro = dblnostressPstemgro;
			dblnostressRootgro = dblnostressProotgro;
			dblnostressLeafgro = dblnostressPleafgro;
		}
		else
		{
			//* reduce demand equally for all organs ala SimPotato (or alter demand so that at least maintResp is satisifed)
			dblnostressLeafgro = dblnostressPleafgro * C_supply_nostress / (dblnostressPleafgro+dblnostressPstemgro+dblnostressProotgro);
			dblnostressStemgro = dblnostressPstemgro * C_supply_nostress / (dblnostressPleafgro+dblnostressPstemgro+dblnostressProotgro);
			dblnostressRootgro = dblnostressProotgro * C_supply_nostress / (dblnostressPleafgro+dblnostressPstemgro+dblnostressProotgro);
	
		}
	}
	//* Tuber growth
	else if ((develop->get_istage() == 2) || (develop->get_istage() == 3)){
		if (info.G1 == 0){ // for crops that are 100% indeterminant, readjust tuber demand so that tubers have last priority on whatever assimilate is available
			if (dblProotgro+dblPleafgro+dblPstemgro <= C_supply)
				dblPtubgro = C_supply - dblProotgro - dblPleafgro - dblPstemgro;
			if (dblPtubgro > tubers->get_tubmax()) dblPtubgro = tubers->get_tubmax();
		}
		// Case a: excess assimilate
		if (C_demand <= C_supply){
			dblStemgro = dblPstemgro;
			dblRootgro = dblProotgro;
			dblLeafgro = dblPleafgro;
			dblTubgro = dblPtubgro;
		}
		if (C_demand_nostress <= C_supply_nostress){
			dblnostressStemgro = dblnostressPstemgro;
			dblnostressRootgro = dblnostressProotgro;
			dblnostressLeafgro = dblnostressPleafgro;
			dblnostressTubgro = dblnostressPtubgro;
		}
		// Case b: insufficient assimilate supply
		if (C_demand > C_supply){
			excess = 1.;
			if ((dblPtubgro >= C_supply)&&(dblPtubgro>0))
			{
				double diff = dblPtubgro - C_supply;
				if (C_pool >= diff) {
					C_supply += diff;
					C_pool -= diff;
				}
				else {
					C_supply += C_pool;
					C_pool = 0.;
				}
				dblPtubgro = C_supply;
				dblProotgro = dblPleafgro = dblPstemgro = 0.;
				dblTubgro = dblPtubgro;
				dblRootgro = dblProotgro;
				dblStemgro = dblPstemgro;
				dblLeafgro = dblPleafgro;
			}
			else 
			{
			//* reduce demand equally for all organs as in original SUBSTOR MODEL !!!! - coudln't get SIMGUI to work here
				dblTubgro = dblPtubgro;
				dblLeafgro = dblPleafgro * (C_supply-dblPtubgro) / (dblPleafgro+dblPstemgro+dblProotgro);
				dblStemgro = dblPstemgro * (C_supply-dblPtubgro) / (dblPleafgro+dblPstemgro+dblProotgro);
				dblRootgro = dblProotgro * (C_supply-dblPtubgro) / (dblPleafgro+dblPstemgro+dblProotgro);
			}	
		}
		//Case b: insufficient assimilate supply - nonwater-stressed example
		if (C_demand_nostress > C_supply_nostress){
			excess = 1.;
			if (dblnostressPtubgro >= C_supply_nostress)
			{
				double diff = dblnostressPtubgro - C_supply_nostress;
				if (C_pool_nostress >= diff) {
					C_supply_nostress += diff;
					C_pool_nostress -= diff;
				}
				else {
					C_supply_nostress += C_pool;
					C_pool = 0.;
				}
				dblnostressPtubgro = C_supply_nostress;
				dblnostressProotgro = dblnostressPleafgro = dblnostressPstemgro = 0.;
				dblnostressRootgro = dblnostressProotgro;
				dblnostressStemgro = dblnostressPstemgro;
				dblnostressLeafgro = dblnostressPleafgro;
			}
			else 
			{
			//* reduce demand equally for all organs as in original SUBSTOR MODEL !!!! - coudln't get SIMGUI to work here
				dblnostressTubgro = dblnostressPtubgro;
				dblnostressLeafgro = dblnostressPleafgro * (C_supply_nostress-dblnostressPtubgro) / (dblnostressPleafgro+dblnostressPstemgro+dblnostressProotgro);
				dblnostressStemgro = dblnostressPstemgro * (C_supply_nostress-dblnostressPtubgro) / (dblnostressPleafgro+dblnostressPstemgro+dblnostressProotgro);
				dblnostressRootgro = dblnostressProotgro * (C_supply_nostress-dblnostressPtubgro) / (dblnostressPleafgro+dblnostressPstemgro+dblnostressProotgro);
			}	
		}
	else if ((develop->get_istage() == 4) || (develop->get_istage() == 5) || (develop->get_istage() == 6)) return;
	else if (develop->get_istage() == 7){
		NitrogenStatus.rootNitrogendemand = NitrogenStatus.optimumrootNitrogenRatio*(this->get_rootMass()+dblRootgro)-NitrogenStatus.rootNitrogenAmount;
		if (NitrogenStatus.rootNitrogendemand < 0) NitrogenStatus.rootNitrogendemand = 0.;
		NitrogenStatus.totalNitrogendemand = NitrogenStatus.rootNitrogendemand;
		return;
		}
	}

	/* Verify N balance and adjust C partitioning to match if necessary*/
	/* Adjust C partitioning based on root C demand from previous time step in 2DSOIL as in MAIZSIM
	/* This is essentially the effect of waterstress on carbon partitioning in which roots get major priority*/
	/* Note: to simulate non-water stress, this needs to be commented out since the amount of C allocated to roots is assumed to be equal to rootgro at previous time-increment*/
	dblRootgro = max(dblRootgro,dblnostressRootgro); /*******************/
	shootPart_virtual = dblnostressLeafgro*Rgleaf + dblnostressStemgro*Rgorgan;

	if (w.pcrs > rootPart_old) 
	{
		//if true, this means 2DSOIL allocated more CHO than was originally allocated to the roots at prior time-step
		// therefore, need to subtract this extra amount to prevent mass balance issue
		// based on 2DSOIL's routine, since roots have 100% C priority during water stress, growth to all other organs is penalized
		// Not sure if tuber growth, hoewver, should be penalized, especially in determinate cultivars, so I don't include tuber C demand in any of these considerations
		// For Root CHO accounting purposes, dblRootgro is adjusted in crop.dll for root mass gained/lost at prior time-step in 2dsoil
		//    rootPart is what is sent to 2DSOIL and does not need to be adjusted for prior time-step!
		double excess = (w.pcrs - rootPart_old)/Rgorgan;
		if ((dblLeafgro + dblStemgro) < excess){ //for mass balance, need to pull additional amount from C_pool then
			C_pool -= (excess - (dblLeafgro + dblStemgro));
			dblLeafgro = 0.;
			dblStemgro = 0.;
			shootPart_real = 0.;
		}
		else { //enough current shoot growth to account for CHO moved to roots at last time-step
			shootPart_real = __max(0.,((dblLeafgro+dblStemgro)-(w.pcrs-rootPart_old)/Rgorgan)); //reduce leaf and stem growth for this iteration
			double temp = dblLeafgro + dblStemgro;
			dblLeafgro = shootPart_real * dblLeafgro / temp; //reduce actual growth for leaf and stem in proportion to their demand at current timestep - don't touch tuber / storage organ growth
			dblStemgro = shootPart_real * dblStemgro / temp;

		}
		rootPart = dblRootgro*Rgorgan;
		rootPart_real = dblRootgro + (w.pcrs-rootPart_old)/Rgorgan; //add additional C from 2DSOIL time-step to roots to maintain mass balance
		dblRootgro = rootPart_real;
	}
	else
	{
		//if false, a check needs to be made whether the CHO allocated to the root at prior time-step was used up.  If not, place in
		//root CHO reserve that will be used by 2DSOIL to grow additional root mass at night.  This may through off the hourly carbon balance in the model a bit, 
		// but seasonally should be okay
		C_pool_root += __max(0.,(rootPart_old - w.pcrs));
		shootPart_real = (dblLeafgro + dblStemgro);
		rootPart_real = dblRootgro - __max(0.,rootPart_old - w.pcrs)/Rgorgan; //subtract out CHO that was already sent to root pool

		rootPart = dblRootgro*Rgorgan;
		dblRootgro = rootPart_real;
	}

	//Now, C from prior time-step was accounted for, do N balance
	//nitrogen_balance_pretubers(info); //evaluate N demand in existing and new organs and adjust C partitioning where necessary
	//nitrogen_balance_posttubers(info); // as above, after TI
	if (info.Nitrogen_stress_off == 1)
	{ 
		nitrogen_stress_off(info); //set N demand N amounts equal
	}
	else
	{	
		nitrogen_balance_pretubers_substor(info); //evaluate N demand in existing and new organs and adjust C partitioning where necessary
		//nitrogen_balance_posttubers(info); // as above, after TI
		nitrogen_balance_posttubers2(info);
	}
	/* Now that C partitioning was adjusted based on N status, finalize overal partitioning information*/
	//Note: this doesn't account for CHO potentially available to roots from the nonstressed haulm growth variables
	//However, mass balance issue won't occure because this excess amount of CHO will be removed from C_pool at next time-step
	C_supply = C_supply - dblTubgro - dblStemgro - dblRootgro - dblLeafgro;
	C_demand = dblTubgro + dblStemgro + dblRootgro + dblLeafgro;
	if (C_supply < 0.00001) { //correct if nonwaterstressed root growth and/or tuber growth used instead of stressed version
		C_pool -= fabs(C_supply);
		C_pool_used += fabs(C_supply);
		C_supply = 0.;//avoid problems with imprecision of floats ?
	}
	if ((dblStemgro+dblRootgro+dblLeafgro)>0) Tubfraction = dblTubgro / (dblStemgro+dblRootgro+dblLeafgro);
	if (Tubfraction > 1) Tubfraction = 1.;
	if (C_supply > 0.) //not sure if this would actually occur?
	{
		C_pool += C_supply;
		C_pool_used -= C_supply;
		C_supply = 0.;

	}
} 

void CPlant::C_allocation2(int iCur, const TWeather &w, const TInitInfo info, double lwpd)
{
	//*********************************
	// This routine allocates overall leaf CHO among leaves within the canopy after growth adjusted for N deficiency
	// Notes:
	//	-assumes Rm of each leaf is satisfied first
	//		-However, in this iteration, Rm has no affect since it's already substracted out of the pool in C_allocation1
	//	-Previus versions assumed growth demand higher for sunlit than shaded leaves
	//		-Sunlit leaves are assumed to be composed of all young leaves, followed by highest level of medium leaves, etc.
	//		-within each leaf class, youngest and highest leaves in the canopy have priority
	//  -Current versions simply assumes all leaves have identical priority for assimilate
	//		-if assimilate limited, each leaf C demand is reduced by same percentage
	//  -N is allocated based on miminum N content for each leaf, then additional N allocated based on luxury N consumption
	//		-in current implementation, no practical affect on the simulation results, just storage for later knowledge gaps
	//**********************************
	double maintTemp = 0., Rl = 0.;
	double sun = 0., shade = 0.;
	double sundemand = 0., shadedemand = 0.;
	double Tempcho = 0., tmpLeafgro = 0., rFactor = 0.;
	double tmpYounggro = 0., tmpMedgro = 0., tmpOldgro = 0.;
	double dt = initInfo.timeStep/(24.*60.);
	double leafcho = 0.;
	double tempNamount = 0., minNreq = 0., deltaNreq = 0., luxNreq = 0.;

	tmpLeafgro = dblLeafgro;
	if (tmpLeafgro < 0) tmpLeafgro = 0.;
	if (initInfo.bigleaf == 0)
	{
		/* Step 0: Determine young, medium, and old demands first*/
		/* Notes:
			-should break this down according to branch type? 	*/
		for (int i = 0; i<inodeNumber; i++)
		{
			if (nodalUnit[i].get_leaf()->isYoung()) tmpYounggro += nodalUnit[i].get_leaf()->get_potentialDWinc()/Rgleaf;
			else if (nodalUnit[i].get_leaf()->isOld()) tmpOldgro += nodalUnit[i].get_leaf()->get_potentialDWinc()/Rgleaf;
			else tmpMedgro += nodalUnit[i].get_leaf()->get_potentialDWinc()/Rgleaf;
		}
		//* Step 1: Young leaf demand satisfied first (sunlit leaves approximated as composed first of all young leaves)
		/* Notes:
				-demand shared equally
				-should add component so that each leaf at least gets Rm
		*/
	/*	
		rFactor = tmpLeafgro / tmpYounggro;  //if < 1 then assimilate supply is limiting
		if (rFactor > 1) rFactor = 1;
		for (int i = 0; i< inodeNumber; i++)
		{
			if (nodalUnit[i].get_leaf()->isYoung()&& tmpLeafgro > 0)
			{		
				Tempcho = nodalUnit[i].get_leaf()->get_potentialDWinc()/Rgleaf*rFactor;
				if (Tempcho > tmpLeafgro)
				{
					Tempcho = tmpLeafgro;
				}
				nodalUnit[i].get_leaf()->import_CH2O(Rgleaf*Tempcho);
				nodalUnit[i].get_leaf()->expand(iCur, info, 1, NitrogenStatus, Rgleaf*Tempcho, develop);
				tmpLeafgro = tmpLeafgro - Tempcho;
				if (tmpLeafgro < 0.000000000000000000001) tmpLeafgro = 0;
				Rg+=(1-Rgleaf)*Tempcho;
			}
		}
		//* Step 2: Determine medium age leaf demand and satisfy that second*/
		/* Notes:
				-all medium leaves have equal priority
				-should add component so that each leaf at least gets Rm
				-should add component so medium aged leaves in different branches have different priorities
		*/
	/*
		if (tmpLeafgro >= tmpMedgro && tmpMedgro > 0) // each leaf can gro according to potential
		{
			for (int i = inodeNumber - 1;i>=0;i--)
			{
				if (!nodalUnit[i].get_leaf()->isYoung() && !nodalUnit[i].get_leaf()->isOld())
				{
					Tempcho = nodalUnit[i].get_leaf()->get_potentialDWinc()/Rgleaf;
					nodalUnit[i].get_leaf()->import_CH2O(Rgleaf*Tempcho);
					nodalUnit[i].get_leaf()->expand(iCur, info, 1, NitrogenStatus, Rgleaf*Tempcho, develop);
					tmpLeafgro = tmpLeafgro - Tempcho;
					if (tmpLeafgro < 0.0000000000000001) tmpLeafgro = 0;
					Rg+=(1-Rgleaf)*Tempcho;
				}
			}
		}
		else if (tmpMedgro > 0)//growth potential is limited, reduce growth of each medium leaf proportionately
		{
			rFactor = tmpLeafgro / tmpMedgro;
			if (rFactor > 1) rFactor = 1;
			for (int i = inodeNumber - 1;i>=0;i--)
			{
				if (!nodalUnit[i].get_leaf()->isYoung() && !nodalUnit[i].get_leaf()->isOld())
				{
					Tempcho = nodalUnit[i].get_leaf()->get_potentialDWinc()/Rgleaf*rFactor;
					nodalUnit[i].get_leaf()->import_CH2O(Rgleaf*Tempcho);
					nodalUnit[i].get_leaf()->expand(iCur, info, 1, NitrogenStatus, Rgleaf*Tempcho, develop);
					tmpLeafgro = tmpLeafgro - Tempcho;
					if (tmpLeafgro < 0.0000000000000001) tmpLeafgro = 0;
					Rg+=(1-Rgleaf)*Tempcho;
				}
			}
		}
		/*Step 4: Remaining assimilate partition evenly among old leaves
		/*Notes:
			-assume each leaf in canopy gets same priority for demand
			-should add component so that each leaf at least gets Rm
			-should add component so medium aged leaves in different branches have different priorities
			-should add component so that leaves that don't get Rm (if assimilate not available) start to senesce
		*/
	/*	
		if (tmpLeafgro > 0)
		{
			rFactor = tmpLeafgro / tmpOldgro;
			if (rFactor > 1) rFactor = 1;
			for (int i = inodeNumber-1;i>=0;i--)
			{
				if (nodalUnit[i].get_leaf()->isOld())
				{
					Tempcho = nodalUnit[i].get_leaf()->get_potentialDWinc()/Rgleaf*rFactor;
					nodalUnit[i].get_leaf()->import_CH2O(Rgleaf*Tempcho);
					nodalUnit[i].get_leaf()->expand(iCur, info, 1, NitrogenStatus, Rgleaf*Tempcho, develop);
					tmpLeafgro = tmpLeafgro - Tempcho;
					if (tmpLeafgro < 0.0000000000000001) tmpLeafgro = 0;
					Rg+=(1-Rgleaf)*Tempcho;
				}
			}
		}

		//instantResp += dblLeafgro * (1-Rgleaf) * initInfo.plantDensity / 30 *1000000 / (initInfo.timeStep*60);
		// Distribute CHO
		/* YOUNG LEAF ALTERNATIVE APPROACH - youngest leaves fet first priority*/
		for (int i=inodeNumber-1; i >= 0; i--)
		{
				Tempcho = nodalUnit[i].get_leaf()->get_potentialDWinc()/Rgleaf;
				if (Tempcho > tmpLeafgro)
				{
					Tempcho = tmpLeafgro;
				}
				if (tmpLeafgro > 0)
				{
					nodalUnit[i].get_leaf()->import_CH2O(Rgleaf*Tempcho);
					nodalUnit[i].get_leaf()->expand(iCur, info, 1, NitrogenStatus, Rgleaf*Tempcho, develop, lwpd, w);
					tmpLeafgro = tmpLeafgro - Tempcho;
					if (tmpLeafgro < 0) tmpLeafgro = 0.;
					Rg+=(1-Rgleaf)*Tempcho;
				}
		}
		//Distribute N
		tempNamount = NitrogenStatus.leafNitrogenAmount;
		//case 1
		//if (NitrogenStatus.actualleafNitrogenRatio > NitrogenStatus.optimumleafNitrogenRatio)
		for (int i = inodeNumber - 1; i >= 0; i--) //First, each leaf gets minimum N content
		{	
			if (!nodalUnit[i].get_leaf()->isTerminated())
			{
				minNreq = nodalUnit[i].get_leaf()->get_drymass() * NitrogenStatus.minimumleafNitrogenRatio;
				if (minNreq > tempNamount) minNreq = tempNamount;
				nodalUnit[i].get_leaf()->set_currentNitrogenAmount(minNreq);
				tempNamount -= minNreq;
			}
		}
		if (tempNamount > 0) //Second, starting with youngest leaf, N is added to optimum level
		{
			for (int i = inodeNumber - 1; i>=0; i--)
			{
				if (!nodalUnit[i].get_leaf()->isTerminated())
				{
					minNreq = nodalUnit[i].get_leaf()->get_drymass() * NitrogenStatus.minimumleafNitrogenRatio;
					deltaNreq = nodalUnit[i].get_leaf()->get_drymass() * (NitrogenStatus.optimumleafNitrogenRatio - NitrogenStatus.minimumleafNitrogenRatio);
					if (deltaNreq > tempNamount) deltaNreq = tempNamount;
					nodalUnit[i].get_leaf()->set_currentNitrogenAmount(deltaNreq+minNreq);
					tempNamount -= deltaNreq;
				}
			}
		}
		if (tempNamount > 0) //Third, luxury N proportioned according to leaf N content
		{
			for (int i = inodeNumber - 1; i>-0; i--)
			{
				if (!nodalUnit[i].get_leaf()->isTerminated())
				{
					minNreq = nodalUnit[i].get_leaf()->get_drymass() * NitrogenStatus.minimumleafNitrogenRatio;
					deltaNreq = nodalUnit[i].get_leaf()->get_drymass() * (NitrogenStatus.optimumleafNitrogenRatio - NitrogenStatus.minimumleafNitrogenRatio);
					luxNreq = (minNreq + deltaNreq) * 0.9;
					if (luxNreq > tempNamount) luxNreq = tempNamount;
					nodalUnit[i].get_leaf()->set_currentNitrogenAmount(minNreq+deltaNreq+luxNreq);
					tempNamount -= luxNreq;
				}
			}
			if (tempNamount > 0) //if still extra N sitting around, store it temporarily in newest leaf!
			{
				minNreq = nodalUnit[inodeNumber-1].get_leaf()->get_drymass() * NitrogenStatus.minimumleafNitrogenRatio;
				deltaNreq = nodalUnit[inodeNumber-1].get_leaf()->get_drymass() * (NitrogenStatus.optimumleafNitrogenRatio - NitrogenStatus.minimumleafNitrogenRatio);
				luxNreq = (minNreq + deltaNreq) * 0.9;
				nodalUnit[inodeNumber-1].get_leaf()->set_currentNitrogenAmount(minNreq + deltaNreq + luxNreq + tempNamount);
				tempNamount = 0.;
			}
		}

		/* OLD LEAF ALTERNATIVE APPROACH - older leaves get first priority*/
		/*
		{
			for (int i=0;i<inodeNumber;i++)
			{
				Tempcho = nodalUnit[i].get_leaf()->get_potentialDWinc()/Rgleaf;
				if (Tempcho > tmpLeafgro)
				{
					Tempcho = tmpLeafgro;
				}
				if (tmpLeafgro > 0)
				{
					nodalUnit[i].get_leaf()->import_CH2O(Rgleaf*Tempcho);
					nodalUnit[i].get_leaf()->expand(iCur, info, 1, NitrogenStatus, Rgleaf*Tempcho, develop);
					tmpLeafgro = tmpLeafgro - Tempcho;
					if (tmpLeafgro < 0) tmpLeafgro = 0;
					Rg+=(1-Rgleaf)*Tempcho;
				}
			}
		}
		*/

		/*ALTERNATIVE APPROACH*/
		/* All leaves have same priority regardless of branch, position, etc*/
		/*
		rFactor = tmpLeafgro / (tmpYounggro+tmpOldgro+tmpMedgro);
		if (rFactor > 1) rFactor = 1;
		if (rFactor >= 0)
		{
			for (int i=0;i<inodeNumber;i++)
			{
				Tempcho = nodalUnit[i].get_leaf()->get_potentialDWinc()/Rgleaf*rFactor;
				if (Tempcho > tmpLeafgro) Tempcho = tmpLeafgro;
				if ((Rgleaf*Tempcho) > nodalUnit[i].get_leaf()->get_potentialDWinc())
				{
					Tempcho = nodalUnit[i].get_leaf()->get_potentialDWinc() / Rgleaf;

				}
				//Tempcho = nodalUnit[i].get_leaf()->get_potentialDWinc()/Rgleaf; // only use this line if you want no C-limited leaf expansion
				nodalUnit[i].get_leaf()->import_CH2O(Rgleaf*Tempcho);
				nodalUnit[i].get_leaf()->expand(iCur, info, 1, NitrogenStatus, Rgleaf*Tempcho, develop);
				tmpLeafgro = tmpLeafgro - Tempcho;
				Rg+=(1-Rgleaf)*Tempcho;
			}
		}
		*/
		
		
	}
	else //big leaf model
	{
		if (!nodalUnit=='\0')
		{
			nodalUnit[0].get_leaf()->import_CH2O(Rgleaf*tmpLeafgro);
			nodalUnit[0].get_leaf()->expand(iCur, info, 1, NitrogenStatus, Rgleaf*tmpLeafgro, develop, lwpd, w);
			nodalUnit[0].get_leaf()->set_currentNitrogenAmount(NitrogenStatus.leafNitrogenAmount);
			Rg += (1-Rgleaf)*tmpLeafgro;
		}
	}
	this->totalLeafNitrogen = NitrogenStatus.leafNitrogenAmount;
	instantResp += dblLeafgro * (1-Rgleaf) * initInfo.plantDensity / 30. *1000000. / (initInfo.timeStep*60.);
	
}


void CPlant::C_allocation3(int iCur, const TWeather &w, const TInitInfo info)
{
	//************************************************
	// This third part alots C and N to all organs after growth adjusted for N deficiency
	// and removes growth respiration
	//************************************************
	if (initInfo.bigleaf == 0)
	{
		for (int i=0; i< inodeNumber; i++)
		{
			if (nodalUnit[i].isInitiated() && !nodalUnit[i].isTerminated())
			{
				nodalUnit[i].get_stem()->import_CH2O(Rgorgan*dblStemgro/inodeNumber);
				nodalUnit[i].get_stem()->set_currentNitrogenAmount(NitrogenStatus.stemNitrogenAmount/inodeNumber);
				
				//if (nodalUnit[i].get_leaf()->isYoung()) //CHO is partitioned evenly among young and old nodalunits
				//{
				//	nodalUnit[i].get_leaf()->import_CH2O(Rgleaf*dblYoungleafgro/iyoungnodeNumber);
				//	nodalUnit[i].get_leaf()->expand(1, dblYoungleafgro/iyoungnodeNumber);
				//}
				//else if (nodalUnit[i].get_leaf()->isGrowing())
				//{
				//	nodalUnit[i].get_leaf()->import_CH2O(Rgleaf*dblOldleafgro/(ioldnodeNumber));
				//	nodalUnit[i].get_leaf()->expand(1, dblOldleafgro/(ioldnodeNumber));
				//}
			}
		}
	}
	else
	{
		if (!nodalUnit=='\0') {
			nodalUnit[0].get_stem()->import_CH2O(Rgorgan*dblStemgro);
			nodalUnit[0].get_stem()->set_currentNitrogenAmount(NitrogenStatus.stemNitrogenAmount);
		}

	}
	this->totalStemNitrogen = NitrogenStatus.stemNitrogenAmount;
	get_tubers()->import_CH2O(Rgorgan*dblTubgro);
	get_tubers()->set_currentNitrogenAmount(NitrogenStatus.tuberNitrogenAmount);
	this->totalTuberNitrogen = NitrogenStatus.tuberNitrogenAmount;
	get_roots()->import_CH2O(Rgorgan*dblRootgro);
	get_roots()->set_currentNitrogenAmount(NitrogenStatus.rootNitrogenAmount);
	this->totalRootNitrogen = NitrogenStatus.rootNitrogenAmount;
	Rg += (1-Rgorgan)*(dblStemgro+dblTubgro+dblRootgro);
	instantResp += (dblTubgro+dblRootgro) * (1-Rgorgan) * initInfo.plantDensity / 30. *1000000. / (initInfo.timeStep*60.);
	shootPart = Rgleaf*dblLeafgro+Rgorgan*dblStemgro;//don't include tubers - assume CHO partitioninge coefficient to tubers unaffected by water status
	tubPart = Rgorgan*dblTubgro;
}


