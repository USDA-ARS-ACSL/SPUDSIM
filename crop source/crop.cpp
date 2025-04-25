// myplant.cpp : Defines the entry point for the DLL application.
//#define MYPLANT_EXPORTS

#include "stdafx.h"
#include "crop.h"
#include "controller.h"
#include "weather.h"

#include "time.h"
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <string>
#include <iomanip>
using namespace std;
#include <cmath>
#define endl "\n"
#define comma ","

// note that we have to dereference the variable in order to
// assign a value that can be passed back to 2DSOIL. This is 
// because the FORTRAN program expects a pointer rather than
// a value. I don't think this applies to structures as it does 
// variables that may be in the arguments list.
// note use of lower case names. Upper and lower case conversions between
// fortran and C++ don't matter here because these are arguments and not
// a function name. CROP must be upper case because it is a function
// name
#ifdef _WIN32
void _stdcall CROP(struct
#else
void crop(struct
#endif
	ShootCommon* SHOOTR,
	WeathCommon* Weather,
	GridCommon* grid_public,
	NodeCommon* node_public,
	ElementCommon* ele_public,
	BoundaryCommon* bound_public,
	TimeCommon* time_public,
	ModuleCommon* module_public,
	FileCommon* file_public
)

{

	// I think pLeaf can be local since the leaf array holds the list
	// and the pointers are not lost between invocations of the procedure
	//First read input data if start of simulation
	char* Buffer = (char*)calloc(256, sizeof(char));
	static int simulation_done = 0; //turn to 1 after maturity reached
	static CController* pSC; //SK, declare as static to ensure only one copy is instantiated during 2DSOIL execution
	if (time_public->lInput == 1) //SK, initialiing crop module
		//*****************************************

	{
		// Parse the file names from the FORTRAN strings passed from 2dsoil
		//KY looks like GNU Fortran handle linebreak differently, making filename detection unusable
		//KY this new macro based on std::string should work on both platforms with smaller code
#define SETSTR(s, n) std::string s(n, sizeof(n)); s.erase(s.find_last_not_of(" \n\r\t")+1);
		SETSTR(varFile, file_public->VarietyFile);
		SETSTR(GraphicFile, file_public->GraphicsFile);

		std::cout << "Initializing crop module..." << endl;
		initInfo.plantDensity = SHOOTR->PopArea;
		initInfo.latitude = Weather->LATUDE;
		initInfo.longitude = Weather->Longitude;
		initInfo.altitude = Weather->Altitude;
		initInfo.CO2 = Weather->CO2;
		initInfo.Seedreserve = SHOOTR->Plantingdrymass;
		initInfo.Plantingdepth = SHOOTR->Plantingdepth;
		initInfo.Sproutlength = SHOOTR->Plantinglength;
		if (time_public->DailyOutput == 1) {
			initInfo.outtimeStep = 1440;
		}
		else {
			initInfo.outtimeStep = 60;
		}
		initInfo.timeStep = time_public->TimeStep;
		initInfo.bigleaf = SHOOTR->Bigleaf;
		initInfo.beginDay = time_public->beginDay;
		initInfo.sowingDay = time_public->sowingDay;
		initInfo.emergeDay = time_public->emergeDay;
		initInfo.endDay = time_public->endDay;
		initInfo.year = time_public->Year;
		initInfo.Nitrogen_stress_off = SHOOTR->Nstressoff;
		initInfo.Water_stress_off = SHOOTR->Wstressoff;
		initInfo.Water_stress_simulation_type = SHOOTR->Wstresstype;

		time_public->iTime = 1;
		SHOOTR->LCAI = 0.0;
		SHOOTR->LAREAT = 0.0;
		SHOOTR->Height = 0.0;
		SHOOTR->Convr = 1.0; //would convert g CHO to g C in 2DSOIL, but all carbon currency should be on CHO basis
		SHOOTR->AWUPS = 0.0;  //initialize AWUPS, AWUPS_old and psil_ in 2DSOIL Yang 8/15/06
		SHOOTR->psil_ = -0.5;  //note, this is in bars!
		SHOOTR->PCRS = 0.0;
		SHOOTR->TRWU_SIM = 0.0; //dhf - total root water extractable by current roots system
		SHOOTR->PSILT_SIM = -15.0; //dhf - threshold LWP atwhich shoot growth can't occur
		SHOOTR->PSISM_SIM = -0.5; //dhf - initial average slab soil water potential
		SHOOTR->ET_demand = 0.0;
		//SHOOTR->HourlyCarboUsed = 0; //is now initialized in crop.h - AD-2011-12-09 hourly CHO used for root growth, also zeroed at initializationi in 2dsoil
		//Period=1.0/24.0;
		Period = (initInfo.timeStep / 60.0) / 24.0; // allow crop.dll to run at user specified time-step
		PopSlab = SHOOTR->PopRow / 100.0 * SHOOTR->EOMult; //note: this is really-> poprow * (100/row-spacing) * (row-spacing / 10000) * EOMult
		Popare = SHOOTR->PopRow * 100.0 / SHOOTR->RowSp;

		// A new plant model object is created and initialized (calls initialize function) here		
		//  ***************************************************************************
		pSC = new CController(varFile.c_str(), GraphicFile.c_str(), initInfo);
		//  ***************************************************************************

				//NitrogenUptake = pSC->getPlant()->get_N()*PopSlab; //initialize nitrogen uptake with what is already in the plant, mg N plant-1 * Popslab = mg N slab-1
		NitrogenUptake = 0; //should be zero at emergence since all plant N has come from seedpiece, DHF
		SHOOTR->NDemandError = 0;
		SHOOTR->CumulativeNDemandError = 0;
		time_public->RunFlag = 1;

		module_public->NumMod = module_public->NumMod + 1;
		ModNum = module_public->NumMod;
		time_public->tNext[ModNum - 1] = time_public->Time + Period;
		//time_public->tnext[Modnum-1]=pSC->getSowingDay();
	} //End Initialization

	if (simulation_done == 1) return; //DHF - check to make sure when crop is matured, dll is exited (avoids error message at end of simulation
	//DHF - Set weather data whether or not plant has emerged yet - this is needed for pre-emergence calculations in 2DSOIL
	TWeather wthr;
	{
		wthr.jday = Weather->JDAY;
		wthr.time = time_public->Time - Weather->JDAY;
		wthr.CO2 = Weather->CO2;
		wthr.airT = Weather->TAIR[time_public->iTime - 1];
		wthr.PFD = Weather->PAR[time_public->iTime - 1] * 4.6; // conversion from PAR in W m-2 to umol s-1 m-2
		wthr.solRad = Weather->WATTSM[time_public->iTime - 1]; //conversion from W m-2 total radiation to J m-2 in one hour 				
		double Es = (0.611 * exp(17.502 * wthr.airT / (240.97 + wthr.airT))); // saturated vapor pressure at airT
		wthr.RH = (1 - (Weather->VPD[time_public->iTime - 1] / Es)) * 100.; // relative humidity in percent
		wthr.rain = Weather->RINT[time_public->iTime - 1];
		wthr.wind = Weather->WIND * (1000.0 / 3600.0); // conversion from km hr-1 to m s-1
		wthr.psil_ = SHOOTR->psil_ / 10; // convert current leaf water potential from bars to MPa
		wthr.dayLength = Weather->DAYLNG;
		/******************************/
		if (initInfo.Water_stress_simulation_type == 1) {
			wthr.swdf1 = 1;	wthr.swdf2 = 1; wthr.swpartition = 0; //set to 0 for empirical water stress scenario*/
		}
		else
			wthr.swdf1 = 1; wthr.swdf2 = 1; wthr.swpartition = 1;
		/******************************/
	}
	//SK, Running the crop module step by step
	if (module_public->NShoot > 0) // cumulate nitrogen and water uptake for each interval prior to 60 minutes
	{
		TotalPotentialRootWaterUptake += SHOOTR->TRWU_SIM * time_public->Step;// g total potential extractable water slab-1 time-step-1
		WaterUptake += SHOOTR->AWUPS * time_public->Step; // g water slab-1 time-step-1 (will be per hour at end of each 60 minute)
		NitrogenUptake += SHOOTR->SIncrSink / 1.0e6; // MaizeSIM implentation, N taken up via roots in this time-step (mass, g N slab-1) in this time step
		// Note that SIncrSink has been multiplied by time step in the 2DSOIL solute uptake routine; the 1000000 scales ug to g
		HourlyCarboUsed = HourlyCarboUsed + SHOOTR->PCRS * time_public->Step; //AD 12-12-2011 put it into crop.cpp

	}
	if (fabs(time_public->Time - time_public->tNext[ModNum - 1]) < fabs(0.001 * time_public->Step)) //only loops through each hour, or 0.042 increment per day
	{
		//if((module_public->NShoot == 0) && (fabs(time_public->Time-pSC->getSowingDay()))<0.001)//If the sowing date has come and there is not plant, let the program know so other calculations are not done
		//switched from sowing to emergence - otherwise, initial root set up in 2DSOIL results in root aging which messes up N and W uptake - will need to fix once routine simulates sprout germination 5-21-2015
		if ((module_public->NShoot == 0) && (fabs(time_public->Time - pSC->getEmergenceDay())) < 0.001)//If the emergence date has come and there is not plant, let the program know so other calculations are not done

		{
			module_public->NShoot = 1;
		}
		if (!pSC->getPlant()->get_roots()->get_Initialized()) //This will only run once at first time-step after sowing (Not Emergence!) - DHF
		{
			// Transports initial carbon of root mass initialized in 2DSOIL (read in with element data) into root organ
			// These are the roots that grew from the seedpiece prior to the date of emergence
			pSC->getPlant()->get_roots()->set_initialized();
			pSC->getPlant()->get_roots()->set_EmergenceData(SHOOTR->TotalRootWeight / PopSlab);
			pSC->getPlant()->set_C_seed(SHOOTR->TotalRootWeight / PopSlab);
			pSC->getPlant()->set_N_reserve(pSC->getPlant()->get_roots()->get_currentNitrogenAmount());
			pSC->getPlant()->set_N_seedused(pSC->getPlant()->get_roots()->get_currentNitrogenAmount());
		}

		double Es = 0.;
		//CurrentNUptakeError = NitrogenUptake / PopSlab - pSC->getPlant()->get_CumulativeNitrogenDemand(); //MAIZSIM g N plant-1 h-1; calculate error for demand and actual uptake, if negative, demand is greater than uptake
		//CurrentNUptakeError = NitrogenUptake / PopSlab - pSC->getPlant()->get_HourlyNitrogenGrowthDemand(); // this should be comparison just for current hourly time-step - note NitrogenUptake reset to zero each hour, so above line was incorrect comparison

		CurrentNUptakeError = NitrogenUptake / PopSlab - pSC->getPlant()->get_HourlyNitrogenDemand(); // this should be comparison just for current hourly time-step - note NitrogenUptake reset to zero each hour, so above line was incorrect comparison
		CumulativeNUptakeError += CurrentNUptakeError;//MAIZSIM, g N plant-1 season-1; ie, this is season long N needs of plant that are not met. Positive value means no error (surplus)
		TWeather wthr;
		{
			wthr.HourlyOutput = time_public->HourlyOutput;
			wthr.DailyOutput = time_public->DailyOutput;
			wthr.jday = Weather->JDAY;
			wthr.time = time_public->Time - Weather->JDAY;
			wthr.CO2 = Weather->CO2;
			if (Weather->CO2 <= 0)
			{
				wthr.CO2 = initInfo.CO2; //can set CO2 in initials for specific simulations where CO2 is constant
			}
			wthr.airT = Weather->TAIR[time_public->iTime - 1];
			wthr.PFD = Weather->PAR[time_public->iTime - 1] * 4.6; // conversion from PAR in W m-2 to umol s-1 m-2
			wthr.solRad = Weather->WATTSM[time_public->iTime - 1]; //conversion from W m-2 total radiation to J m-2 in one hour 				
			Es = (0.611 * exp(17.502 * wthr.airT / (240.97 + wthr.airT))); // saturated vapor pressure at airT
			wthr.RH = (1 - (Weather->VPD[time_public->iTime - 1] / Es)) * 100.0; // relative humidity in percent
			wthr.rain = Weather->RINT[time_public->iTime - 1];
			wthr.wind = Weather->WIND * (1000.0 / 3600.0); // conversion from km hr-1 to m s-1
			wthr.dayLength = Weather->DAYLNG;
			wthr.psil_ = SHOOTR->psil_ / 10;  // convert current leaf water potential from bars to MPa
			/****************/
			if (pSC->getInitInfo().Water_stress_simulation_type == 1) { //SIMPOTATO version
				wthr.swdf1 = 1; //DHF - legacy SIMGUI soil water deficit factors 1 and 2.  A value of 1.0 indicates no stress
				wthr.swdf2 = 1; // as above, but not implemented here...
				wthr.swpartition = 0; // change to 1 if using root allocation from canopy during water stress
				wthr.trwu_simpotato = TotalPotentialRootWaterUptake * (1 / PopSlab); //need to check units, should be g or cm3 water plant-1 h-1, units of TRWU_sim were g H2O slab-1 d-1
				wthr.psilt_simpotato = SHOOTR->PSILT_SIM / 10;
				wthr.psism_simpotato = SHOOTR->PSISM_SIM / 10;

			}
			else //SPUDSIM AD version
				wthr.swdf1 = 1; wthr.swdf2 = 1; wthr.swpartition = 1.;
			/*****************/

			if (abs(wthr.time - 0.2083) < 0.0001) //If time is 5 am, then pass the leaf water potential (the predawn leaf water potential) from SHOOTR to the wthr object. YY
			{
				//lwpd=SHOOTR->psil_; //Here psil_ is in bar. Since the LWPeffect in leaf.cpp uses leaf water potentialin bar, so here lwpd is in bar, instead of being scaled to MPa. YY
				lwp_dawn = SHOOTR->psil_;
				//lwp_dawn = wthr.psil_; // predawn leaf water potential in MPa for consistency
			}
			// AD 12-12-2011 put HourlyCarboUsed back
			wthr.pcrs = HourlyCarboUsed / PopSlab; // pass CHO used for root growth in 2DSOIL back to plant, units g CH2O plant-1 h-1
			HourlyCarboUsed = 0.0;  //HourlyCarbonUsed is total amount of CHO used for root grwowth by 2dsoil in prior soil iteration loop, dividing by popslab coverts units to g plant-1 h-1
			wthr.TotalRootWeight = SHOOTR->TotalRootWeight / PopSlab; //2DSOIL update of current rootgrowth - g CHO plant-1 season-1
			wthr.MaxRootDepth = SHOOTR->MaxRootDepth;
			// Available water is cm per profile - should be divided by PopSlab
			wthr.ThetaAvailRZ = node_public->ThetaAvailRZ / PopSlab;

			if (NitrogenUptake > 0) //Pass through nitrogen uptake - this is  g N slab-1 hour-1
			{
				//Note, N mobilized from leaf senescence is done in senescence routine, DHF
				//pSC->getPlant()->set_AvailableNitrogen((NitrogenUptake - NitrogenUptakeOld)/PopSlab); //New g N plant-1 taken by roots in past time-increment, DHF
				pSC->getPlant()->set_AvailableNitrogen(NitrogenUptake / PopSlab); // DHF - different from MAIZSIM - logic with naming convention too difficult to follow - this should be current N uptake for the hour only
				pSC->getPlant()->set_CumulativeNitrogenSoilUptake(NitrogenUptake / PopSlab); // total g N plant-1 since emergence, including current time-step
			}
			//double HourlyActualNFromSoil = (NitrogenUptake-NitrogenUptakeOld)/PopSlab; //hourly rate per day, g N plant-1 hour-1
			double HourlyActualNFromSoil = NitrogenUptake / PopSlab; //hourly rate for current uptake, g N plant-1 h-1
			pSC->getPlant()->set_HourlyNitrogenSoilUptake(HourlyActualNFromSoil);
			//SLNmin base Specific leaf nitrogen content; for now assume it's 0.5 YY
			if (SHOOTR->LAI == 0)
			{
				wthr.ET_supply = 0;
				wthr.swdf1 = 1;
				wthr.swdf2 = 1;
			}
			else
			{
				//wthr.ET_supply = WaterUptake*((1/18.01)/((24*3600)*(SHOOTR->RowSp*SHOOTR->EOMult/10000)*SHOOTR->LAI));   //into MAIZESIM Yang 8/15/06
				wthr.ET_supply = WaterUptake / (SHOOTR->EOMult * SHOOTR->PopRow) * 100; // units are now area gram per plant per hour
				wthr.ET_supply = wthr.ET_supply * pSC->getInitInfo().plantDensity / pSC->getPlant()->get_greenLAI() / 18 / 3600; // units are now mol m-2 leaf s-1
				// ET_supply is used in energy balance calculation for leaf temperature and must be in units of mol m-2 leaf s-1
				// as far as I can tell, ET_supply is not being utilized anywhere else in MAIZSIM
				// In SPUDSIM, ratio of ET_supply and transpiration rate used to indicate reduce leaf expansion, therefore need to use same units for both
				//	To do this multiply by EOMULT to double slab width if plant is at the edge. Then multiply by 100/PopRow
				//	to get area inhabited by the plant. This provides a per plant estimate from area.
				// Finally, note that, from plant side, get_ET is in mmol h2o m-2 ground-1 s-1 and needs to be converted as well
			}
		}
		pSC->getPlant()->set_TotalSoilWaterUptake(WaterUptake / PopSlab); //g H2O plant-1
		pSC->getPlant()->set_HourlySoilWaterUptake(WaterUptake / PopSlab); //g H2O plant-1 h-1
		pot_transpiration_old = pSC->getPlant()->get_ET() / 1000 / pSC->getPlant()->get_greenLAI(); //g H2O plant-1 h-1
		pot_transpiration_old = pSC->getPlant()->get_ET() * 0.018 / pSC->getInitInfo().plantDensity * 3600;// g h2o plant-1 h-1
		pSC->getPlant()->set_HourlyTranspirationSupplyDemandRatio(wthr.ET_supply / pot_transpiration_old);
		/*****************************************/
		if (pot_transpiration_old > 0) {
			//wthr.swdf1 = 0.67* wthr.ET_supply / pot_transpiration_old;
			//wthr.swdf2 = wthr.ET_supply/pot_transpiration_old;
			wthr.swdf1 = min(wthr.trwu_simpotato / actual_transpiration_old, 1.0);
			if ((wthr.trwu_simpotato / actual_transpiration_old) < 1.5)
				wthr.swdf2 = min(0.67 * wthr.trwu_simpotato / actual_transpiration_old, 1.);
			else
				wthr.swdf2 = 1.0;
		}
		else {
			wthr.swdf1 = 1; wthr.swdf2 = 1;
		}
		/*****************************************/



		int end = grid_public->NumNP;
		double soilT = 0, maxY = 0;
		int count = 0;
		double LowerBoundary = 5;
		// first find top of grid
		for (int i = 0; i < end; i++)
		{
			if (grid_public->y[i] > maxY) maxY = grid_public->y[i];
		}
		LowerBoundary = maxY - LowerBoundary;
		// now find average temperature in layer between surface and lower boundary
		for (int i = 0; i < end; i++)
		{
			if (grid_public->y[i] >= LowerBoundary)
			{
				soilT = soilT + node_public->Tmpr[i];
				count++;
			}
		}
		wthr.soilT = soilT / count;
		int ier = pSC->getErrStatus(); //Iterate through crop model if no errors are found.  If plant isn't emerged yet, nothing should happen, DHF 12-10-2010
		if (ier == 0)
		{
			ier = pSC->run(wthr, lwp_dawn); //Pass both weather and leaf water potential into the "run" functionof the controller pSC YY								
		}
		// Assumes that germination takes place about halfway through the sowing date
		// Below simply spits out current jday just once each 24 hour period  DHF 12-10-2010
		if (wthr.time >= 0.49 && wthr.time <= 0.51)
		{
			if (!pSC->getPlant()->get_develop()->Germinated())
			{
				//cout << wthr.jday <<endl;
			}
		}
		if (pSC->getPlant()->get_develop()->Emerged()) // pass appropriate data to 2DSOIL file structures 
		{
			double pool = pSC->getPlant()->get_C_pool_root(); //Stores any carbon not used for root growth in previous time-step
			double haulm = pSC->getPlant()->get_shootPart();
			if (pSC->getPlant()->get_PhotosyntheticFeedback() < 0.9) { //check if virtual haulm pool is available for increased root growth
				double haulm1 = pSC->getPlant()->get_shootPart();
				double haulm2 = pSC->getPlant()->get_shootPart_virtual();
				haulm = max(haulm1, haulm2);
			}

			if ((pSC->getPlant()->get_C_pool_root() > 0) && (pSC->getPlant()->get_rootPart() < 0.00001)) //Assures root pool only used at night
				//DHF get_rootPART() criteria above might not work as intended - if C_pool is large enough, have root growth at night,
			//	plus, can have periods during sunlight where all CHO is allocated to other parts of plant
			//if ((pSC->getPlant()->get_C_pool_root()>0) && wthr.PFD < 10.0)
			{
				SHOOTR->PCRL = (pSC->getPlant()->get_rootPart() + pool) * 24 * PopSlab;

				SHOOTR->PCRQ = (pSC->getPlant()->get_rootPart() + pool + haulm * wthr.swpartition) * 24 * PopSlab;
				if (pSC->getInitInfo().Water_stress_off == 1) SHOOTR->PCRQ = (pSC->getPlant()->get_rootPart() + pool) * 24 * PopSlab;
				//SHOOTR->PCRQ = (pSC->getPlant()->get_rootPart() + pool)*24*PopSlab;
				pSC->getPlant()->set_C_pool_root(0.0);
			}
			else
			{
				SHOOTR->PCRL = (pSC->getPlant()->get_rootPart()) * 24 * PopSlab;
				SHOOTR->PCRQ = (pSC->getPlant()->get_rootPart() + haulm * wthr.swpartition) * 24 * PopSlab;
				if (pSC->getInitInfo().Water_stress_off == 1) SHOOTR->PCRQ = (pSC->getPlant()->get_rootPart()) * 24 * PopSlab;
				//SHOOTR->PCRQ = (pSC->getPlant()->get_rootPart() + pool)*24*PopSlab;
			}

			//ActualCarbonIncrement is calculated from "assimilate", which in turn is calculated from photosynthsis_net in
			//plant; the unit of assimilate then is in g/plant/hour, thus, at this point, pcrl has unit g/plant/hour
			// Multiply by 24 to get g plant-1 day-1; multiply by popslab to get g Carbo slab-1 day-1
			//SHOOTR->PCRQ=(pSC->getPlant()->get_rootPart())*24*PopSlab; //DHF = trial to prevent more C to go to roots no matter the waterstress
			SHOOTR->LCAI = pSC->getPlant()->get_greenLeafArea() * pSC->getInitInfo().plantDensity / (100 * 100);
			SHOOTR->Cover = 1.0 - exp(-0.79 * SHOOTR->LCAI);
			SHOOTR->Shade = (float)SHOOTR->Cover * SHOOTR->RowSp;
			SHOOTR->Height = min(SHOOTR->Shade, SHOOTR->RowSp);
			SHOOTR->ET_demand = (pSC->getPlant()->get_ET() * 0.018 * 3600 * 24 / pSC->getInitInfo().plantDensity);//pass ET demand from shoot to root, g plant-1 d-1
			actual_transpiration_old = pSC->getPlant()->get_HourlySoilWaterUptake(); //water uptake from soil g plant -1 d-1
			if (pSC->getInitInfo().Water_stress_off == 1)
			{
				//SHOOTR->ET_demand = 0; //keep water movement from soil turned on, just turn off physical stress effects on plant in above ground
			}
			/*In GasExchange, the unit of ET is mmol m-2(leaf) sec-1
			In Plant, ET is multiplied by LAI, so, the unit of ET changes to mmol m-2(ground) sec-1
			need to convert to grams plant-1
			Here, multiplying ET by 0.018 and 3600*24 converts it to g m-2(ground) day-1
			dividing it by plantdensity converts it to g plant-1 day-1 */
			SHOOTR->LAI = (float)pSC->getPlant()->get_greenLeafArea() * (float)pSC->getInitInfo().plantDensity / (100 * 100);
			//if (pSC->getPlant()->get_Water_stress_off()) SHOOTR->LAI=0.01;
			//shoot_weightPerM2 = (pSC->getPlant()->get_leafMass() + pSC->getPlant()->get_stemMass()) * pSC->getInitInfo().plantDensity; //total shoot mass
			//massIncrease = (shoot_weightPerM2 - old_shoot_weightPerM2); //calculated increase in above-ground biomass per m2
			//double shoot_weight = (pSC->getPlant()->get_stemMass()+pSC->getPlant()->get_leafMass())*pSC->getInitInfo().plantDensity/(100*100); //Calculate total shoot mass YY
			//The biomass  returned by getPlant()->get_shootMass() is the weight of each single plant (g/plant), 
			//to convert it into (g/m-2), it has to be 

			//set up accounting of N here.
			//first hourly/Actual and needed N uptake i nthe last hour per plant per day
			//double HourlyNitrogenDemand = max(U_P, 0)/pSC->getInitInfo().plantDensity/24.0; //determined N demand (eqn 1 Lindquist et al., 2007) in g plant-1
			double HourlyNitrogenDemand = pSC->getPlant()->get_HourlyNitrogenDemand();// This would be current N amount needed
			//pSC->getPlant()->set_CumulativeNitrogenSoilUptake(NitrogenUptake/PopSlab); // g N plant-1
			double OldNDemand = 0.;
			OldNDemand = SHOOTR->NitroDemand / PopSlab / 1e6 / 24;
			SHOOTR->NitroDemand = (float)HourlyNitrogenDemand * (float)PopSlab * 1e6 * 24.0; //pass the N demand into 2dsoil, units are ug N slab-1 d-1
			if (pSC->getInitInfo().Nitrogen_stress_off == 1) SHOOTR->NitroDemand = 0;
			//old_shoot_weightPerM2 = shoot_weightPerM2; //save the value of the above_ground biomass of this time-step
			NitrogenUptakeOld = NitrogenUptake; // save the cumulative N uptake from this time-step, g N slab-1
			SHOOTR->NDemandError = (float)CurrentNUptakeError; //but these units are g N plant-1 ?
			SHOOTR->CumulativeNDemandError = (float)CumulativeNUptakeError; //these units are g N plant-1 - why is nitroDemand in ug N slab-1?
		}

		if (pSC->getPlant()->get_develop()->Matured()) //End of hour increment - if plant is mature then end, if not continue
		{
			cout << "Completing crop simulation..." << endl;
			simulation_done = 1;
			module_public->NShoot = 0; //tell 2dsoil that crops harvested
			time_public->tNext[ModNum - 1] = 1e12; // set the next time to run the plant module far later so that it's not gonna run again this year, temporary fix
			if (pSC != NULL) // if matured points to nothing
			{
				delete pSC;
				pSC = NULL;
			}
			time_public->RunFlag = 0;
		}
		else
		{
			time_public->tNext[ModNum - 1] = time_public->Time + Period;
			WaterUptake = 0;
			NitrogenUptake = 0;
			TotalPotentialRootWaterUptake = 0;
		}

	} 	// end hourly calculations code
	// end if Nshoot > 0 section
	return;
}