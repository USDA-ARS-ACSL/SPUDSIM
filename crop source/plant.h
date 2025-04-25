#pragma once
#ifndef _PLANT_H_
#define _PLANT_H_
#include "organ.h"
#include "nodalunit.h"
#include "development.h"
#include "roots.h"
#include "tubers.h"
#include "gas_exchange_new.h"
#include "nitrogenstatus.h"

struct TStage //* Is this needed for potato? - have integer scale for development and double scale following tuber initiation
{
public:
	TStage() { V = 0.0; R = 0.0;}
	double V, R;
};

class CPlant
{
public:
	CPlant();
	CPlant(const TInitInfo);
	~CPlant();

	CNodalUnit* get_nodalUnit() {return nodalUnit;}
	CTubers* get_tubers() {return tubers;}
	CRoots* get_roots() {return roots;}
	CDevelopment * get_develop() {return develop;}
	CGas_exchange_new * get_sunlit() {return this->sunlit;}
	CGas_exchange_new * get_shaded() {return this->shaded;} // get access to the pointers that point to the sunlit/shade leaves Yang 8/29/06

	double lmass[32], page[32], cage[32], area[32]; //public variables for potato individual leaf output routine

	double test = 0.;
	int get_nodeNumber() {return inodeNumber;}
	int get_finalNodeNumber() {return finalNodeNumber;}
	int get_basalleafNo() {return basalleafNo;}
	int get_basalstemNo() {return basalstemNo;}
	int get_apicalleafNo() {return apicalleafNo;}
	int get_apicalstemNo() {return apicalstemNo;}
	int get_mainleafNo() {return mainleafNo;}
	double get_averagegs() {return average_gs;} //average 24hour canopy stomtatl conductance
	double calcPotentialLeafArea();
	double get_age() {return age;}
	double get_averageLWP() {return average_LWP;} // average 24hour bulk canopy water potential

	double get_CH2O() {return CH2O;}
	
	double get_CumulativeNitrogenDemand() {return NitrogenStatus.CumulativeNitrogenDemand;} //return total N demand over course of season, g N plant-1 (updated only in crop.cpp)
	double get_HourlyNitrogenDemand() {return NitrogenStatus.HourlyNitrogenDemand;} // return N demand for current time-step, g N plant-1 (updated in crop.cpp)
	double get_HourlyNitrogenSoilUptake() {return NitrogenStatus.HourlyNitrogenSoilUptake;}
	double get_HourlyNitrogenSeedUsed() {return N_hourlyseedused;}//return N used by plant N demand for current hour, g N plant- 1hour-1
	double get_HourlyNitrogenTranslocated() {return N_hourlytranslocated;} //return N translocated from organs, with buffers below optimum level, for current hour for growth, g N plant-1 hour-1
	double get_HourlyNitrogenGrowthDemand() {return N_hourlygrowthdemand;} // total N required for N growth and for existing organs to get to optimum level, g N plant-1 hour-1
	double get_HourlySoilWaterUptake() {return Hourly_soil_water_uptake;}
	double get_HourlyTranspirationSupplyDemandRatio() {return Hourly_transpiration_supply_demand_ratio;}

	double get_CumulativeNitrogenSoilUptake() {return NitrogenStatus.CumulativeNitrogenSoilUptake;}
	double get_Ndemand() {return NitrogenStatus.totalNitrogendemand;}//return N demand to 2DSOIL from current time step, g N plant-1
	double get_N() {return TotalNitrogen;} // Total N in plant g plant-1
	double get_reserveSeedNitrogen() {return N_reserve;} // remaining N content in seedpiece available for growth demand, g N plant-1
	double get_seedNitrogenUsed() {return N_seedused;} // how much Nseed was used in current time-stem, g N plant-1 h-1
	double get_Nstressfactor() {return NitrogenStressFactor;}
	double get_Nstressfactorone() {return NitrogenStatus.Nitrogendeficiencyone;} //photosynthetic reduction factor, 0 to 1
	double get_Nstressfactorone_24h() {return Nstressfactorone_24h;} //photosynthetic reduction factor, 0 to 1, 24h average
	double get_Nstressfactortwo() {return NitrogenStatus.Nitrogendeficiencytwo;} //root growth factor, 0 to 1
	double get_Nstressfactortwo_24h() {return Nstressfactortwo_24h;} //root growth factor, 0 to 1 24h average
	double get_psi_leafexpansion_stress_24h() {return psi_leafexpansion_stress_24h;}//effect of leaf water potential (0 to 1) on leaf expansion rate 24h average
	double get_psi_leafexpansion_stress() { return psi_leafexpansion_stress; }//effect of leaf water potential (0 to 1) on leaf expansion rate 24h average
	double get_NCRatioTotalPlant() {return NCRatioTotalPlant;}
	double get_NCRatioLeaves() {return NCRatioLeaves;}
	double get_NCRatioStems() {return NCRatioStems;}
	double get_NCRatioRoots() {return NCRatioRoots;}
	double get_NCRatioTubers() {return NCRatioTubers;}

	double get_C_pool() {return C_pool;}
	double get_C_pool_root() {return C_pool_root;}
	double get_C_deadpool_hourly() {return C_deadpool;}
	double get_C_deadpool_season() {return C_deadpool_season;}
	double get_C_reserve() {return C_reserve;}
	double get_C_seed() {return C_seed;}
	double get_C_supply() {return C_supply;}
	double get_C_seed_used() {return C_seed_used;}
	double get_C_pool_used() {return C_pool_used;}

	double get_conductance() {return conductance;} //stomatal conductance


	double get_Parave() {return Parave;}
	double get_Sradave() {return Sradave;}
	double get_Tdayave() {return Tdayave;}
	double get_Tdaylag() {return Tdaylag;}
	double get_Pg() {return photosynthesis_gross;}//instantaneous flux (umol C02 m-2 ground s-1)
	double get_Pn() {return Pnet;}
	double get_Pnleafshaded() {return photosynthesis_netshadedleaf;} //flux umol co2 m-2 leaf s-1
	double get_Pnleafsunlit() {return photosynthesis_netsunlitleaf;} //as above
	double get_Transpirationshaded() {return transpiration_shadedleaf;} // mmol H2O m-2 leaf s-1
	double get_Transpirationsunlit() {return transpiration_sunlitleaf;} // mmol H2O m-2 leaf s-1
	double get_Pg900() {return photosynthesis_gross900;}//instant flux between 850 and 950 during a day (umol co2 m-2 s-1)
	double get_Pn900() {return photosynthesis_net900;}//instant "
	double get_Pleaf900() {return photosynthesis_grossleaf900;}//as above put per m-2 leaf
	double get_Pgsunlitleaf() {return photosynthesis_grosssunlitleaf;} // umol m-2 leaf s-1 
	double get_Pgshadedleaf() {return photosynthesis_grossshadedleaf;} // umol m-2 leaf s-1
	double get_internalCO2() {return internal_CO2;}//ppm
	double get_PhotosyntheticFeedback() {return PhotosyntheticFeedback;}

	double get_Pn2() {return photosynthesis_net;}
	double get_Pg2() {return Pgross;}//cumulative g cho plant-1 d-1 (only resets each 24 h)
	double get_Rg() {return Rg;}//cumulative g cho plant-1 d-1 
	double get_Rm() {return Rm;}//"
	double get_Resp() {return instantResp;} // instantaneous Rm+Rg flux (umol C02 m-2 ground s-1)
	double get_assimilate() {return assimilate;}
	double get_ET() {return transpiration;} // instantaneous transpiration demand (mmol H2O m-2 ground s-1)
	double get_cumulativeET() {return cumulativeET;} // cumulative transpiration over 24 h period (g H2O m-2 ground d-1)
	double get_waterstressfactor() {return swdf1;} // ratio of soil water uptake to ET demand for 24 hour time-period
	
	double get_psistress_gs_factor() {return psistress_gs_factor;} //water stress factor based on leaf water potential effect on stomatal conductance, 0 to 1
	double get_psistress_gs_factor_24h() {return psistress_gs_factor_24h;} //water stress factor based on leaf water potential effect on stomatal conductance, 0 to 1, 24h average
	double get_tmpr() {return temperature;} // hourly canopy temp, C
	double get_tmpr2() {return Tcanopyave;}// average 24h canopy temp, C

	double get_heat_veg_stressfactor() {return dblheatstressVegetative; } // return heat stress indication during vegetative growth
	double get_heat_repro_stressfactor() {return dblheatstressReproductive; }
	
	double get_totalMass() {return totalMass;}
	double get_shootMass() {return stemMass + leafMass;}
	double get_totalLeafNitrogen() {return NitrogenStatus.leafNitrogenAmount;} //same as maizsim get_LeafN()
	double get_totalStemNitrogen() {return NitrogenStatus.stemNitrogenAmount;}
	double get_totalRootNitrogen() {return NitrogenStatus.rootNitrogenAmount;}
	double get_totalTuberNitrogen() {return NitrogenStatus.tuberNitrogenAmount;}
	double get_totalDeadNitrogen() {return totalDeadNitrogen;}
	double get_totalPlantNitrogen() {return TotalNitrogen;}
	double get_totalReserveNitrogen() {return totalReserveNitrogen;}
	double get_totalSoilWaterUptake() {return Total_soil_water_uptake;}
	double get_transpirationSupplyDemandRatio() {return Transpiration_supply_demand_ratio;}
	double get_temperature() {return temperature;}
	double get_stemMass() {return stemMass;}
	double get_leafMass() {return leafMass;}
	double get_deadMass() {return deadMass;}
	double get_deadLeafArea() {return deadLeafArea;}
	double get_leafArea() {return leafArea;}
	double get_greenLeafArea() {return greenLeafArea;}
	double get_greenLAI() {return greenLAI;}
	double get_LAI() {return LAI;}
	double get_sunLAI() {return sunLAI;}
	double get_shadeLAI() {return shadeLAI;}
	double get_sunPFD() {return sunPFD;}
	double get_shadePFD() {return shadePFD;}
	double get_PFD(){return PFD;}
	double get_shootPart() {return shootPart;}
	double get_shootPart_virtual() {return shootPart_virtual;}
	double get_tubPart() {return tubPart;}
	double get_tuberMass() {return tuberMass;}
	double get_rootMass() {return rootMass;}
	double get_rootPart() {return rootPart;}
	double get_rootGrowth() {return dblRootgro;}
	double get_mainstemMass() {return mainstemMass;}
	double get_basalstemMass() {return basalstemMass;}
	double get_apicalstemMass() {return apicalstemMass;}
	double get_mainleafMass() {return mainleafMass;}
	double get_basalleafMass() {return basalleafMass;}
	double get_apicalleafMass() {return apicalleafMass;}
	double get_mainleafArea() {return mainleafArea;}
	double get_basalleafArea() {return basalleafArea;}
	double get_apicalleafArea() {return apicalleafArea;}
	double get_potentialleafGro() {return dblPleafgro;}
	double get_potentialstemGro() {return dblPstemgro;}
	double get_potentialrootGro() {return dblProotgro;}
	double get_potentialtuberGro() {return dblPtubgro;}
	double get_SRAD() {return SRAD;}
	double get_VPD() {return VPD;}


	TStage get_stage() {return stage;}

	void reset();
	void set_C_seed(double x) {C_seed-=x;}
	void set_C_pool_root(double x) {C_pool_root = x;}
	void set_N_reserve(double x) {N_seed-=x;}
	void set_N_seedused(double x) {N_seedused += x;}
	void set_N_hourlyseedused(double x) {N_hourlyseedused = x;}
	void set_N_hourlytranslocated(double x) {N_hourlytranslocated = x;}
	void set_N_hourlygrowthdemand(double x) {N_hourlygrowthdemand = x;}
	void set_mass();
	void set_SRAD(double x) {SRAD=x;}
	void set_leafArea(double x) {leafArea=x;}
	void set_greenLeafArea(double x) {greenLeafArea=x;}
	void set_deadLeafArea(double x) {deadLeafArea = x;}
	void set_greenLAI(double x) {greenLAI=x;}
	void set_deadLAI(double x) {deadLAI=x;}
	void set_LAI(double x) {LAI = x;}
	void set_sunLAI(double x){sunLAI=x;}
	void set_shadeLAI(double x){shadeLAI=x;}
	void set_sunPFD(double x){sunPFD=x;}
	void set_shadePFD(double x){shadePFD=x;}
	void set_PFD(double x){PFD=x;}
	void set_sunPg(double x){sunPg = x;}
	void set_shadePg(double x){shadePg=x;}
	void set_age(double x) {age=x;}
	//void set_CH2O();
	void set_TotalPlantNitrogen(double x) {TotalNitrogen=x;}// g N plant-1 - total N in entire plant including available N
	void set_TotalSoilWaterUptake(double x) {Total_soil_water_uptake+=x;} // g H2O plant-1 from 2DSOIL
	void set_TranspirationSupplyDemandRatio(double x) {Transpiration_supply_demand_ratio=x;} // ratio of water from root to canopy transpiration 
	void set_HourlySoilWaterUptake(double x) {Hourly_soil_water_uptake=x;} // g h2o plant-1 h-1 from 2dsoil
	void set_HourlyTranspirationSupplyDemandRatio(double x) {Hourly_transpiration_supply_demand_ratio = x;} // hourly ratio of water supplied by roots to canopy transpiration demand
	void set_AvailableNitrogen(double x) {NitrogenStatus.availableNitrogen+=x;} //adds new g N plant-1 from 2DSOIL to pool available for growth
	void set_HourlyNitrogenDemand (double x) {NitrogenStatus.HourlyNitrogenDemand=x;}
	void set_CumulativeNitrogenDemand (double x) {NitrogenStatus.CumulativeNitrogenDemand+=x;}
	void set_HourlyNitrogenSoilUptake(double x) {NitrogenStatus.HourlyNitrogenSoilUptake=x;}
	void set_CumulativeNitrogenSoilUptake(double x) {NitrogenStatus.CumulativeNitrogenSoilUptake+=x;} //total g N plant-1 since emergence in plant
	void set_dayinc(double x) {dayinc=x;}
	void set_VPD(double x) {VPD = x;} //sets VPD from gas exchange routine
	
	void update(int iCur, const TWeather &, double lwpd, const TInitInfo, double dayinc);
	void update_nodes(int iCur, const TWeather&, const TInitInfo, double lwpd); //updates individual nodes
	void update_tubers(double Tubfraction);//update tuber sink strength on pgross parameters
	void update_roots();//N status
	void initiate(int iCur, const TWeather& weather);
	void calcGasExchange(const TWeather & weather, const TInitInfo info);
	void calcMaintRespiration(const TWeather&, const TInitInfo info);
	void calcMaxPoolSize(); //what is maximum soluble C that can be stored in leaves?  feedback inhibition effect as a result?
	void calcLeafArea(const TInitInfo);
	void nitrogen_balance_pretubers (const TInitInfo); //determine N demand of new and existing biomass - SIMPOTATO based
	void nitrogen_balance_pretubers_simplified (const TInitInfo); //N demand should be based on new biomass - SIMPOTATO heavily modified for SPUDSIM
	void nitrogen_balance_pretubers_substor (const TInitInfo); // determine N demand and allocation - SUBSTOR coefficients, modified routine
	void nitrogen_balance_posttubers (const TInitInfo); // as above, post-TI
	void nitrogen_balance_posttubers2 (const TInitInfo); // as above, post-TI, simplified version
	void nitrogen_balance_posttubers_substor (const TInitInfo); //



	//void nitrogen_balance_posttubers_simplified (const TInitInfo); // as above, post-TI, but heavily modified for SPUDSIM so N demand based on new biomass only
	void nitrogen_stress (const TInitInfo); //N stress factors, SIMPOTATO based
	void nitrogen_stress_off (const TInitInfo); //Turn off effects of N stress
	void nitrogen_stress_substor (const TInitInfo); //N stress factors, SUBSTOR based
	void dailyave(int iCur, const TWeather &, const TInitInfo, double dayinc); //daily weather values
	//void grow();
	void potential_growth(const TInitInfo, const TWeather &); //summarize potential dw gains for each organ class
	//void C_allocation(const TWeather&, double carbonDemand); //maizesim routine from YY
	void C_allocation1(int iCur, const TWeather&, const TInitInfo info, double lwpd); //compute amount partitioned to all organs
	void C_allocation2(int iCur, const TWeather&, const TInitInfo info, double lwpd); //divy up C among canopy leaves
	void C_allocation3(int iCur, const TWeather&, const TInitInfo info); //partition C to all organs
	void leaf_output(); //output information on individual leaves


private:
	TInitInfo initInfo;
	TNitrogenStatus NitrogenStatus;
	CNodalUnit * nodalUnit;
	CTubers * tubers;
	CRoots * roots;
	CDevelopment * develop;
	CGas_exchange_new * sunlit;
	CGas_exchange_new * shaded; //declare two pointers that point to a sunlit and a shaded leaf Yang 8/29/06 

	int counter = 0; //keeps track of current time incremeent during day
	TStage stage;

	// * Whole plant status indicators *
	// Carbon balance / tracking variables:
	double assimilate_old = 0.; // assimilate fixed from last time-step
	double C_pool = 0.; // shorterm C pool, g(CH2O), includes all assimlate + mobilized C from deadpool
	double C_pool_nostress = 0.; //as above, except this is virtual pool to adjust root growth when plant is under moderate to severe drought
	double C_new_organ = 0.; //CHO available for new organ growth based on difference between houlry C_supply and C_demand, reset to 0
	double C_pool_root = 0.; // storage for CHO allocated to roots but not used in the previous time-step
	double C_pool_room = 0.; // amount of room in leaves available for storage of additional unallocated CHO
	
	double C_deadpool = 0.; //temporary pool for dead leaf mass returned to assimilate after senescence, hourly
	double C_deadpool_season = 0.; //cumulative dead leaf mass returned to assimilate after senescence for entire season
	double C_reserve = 0.; // longterm C pool - not being used right now, 2009
	double C_maxPoolSize = 0.; //maximum size of C pool, g CH2O
	double C_seed = 0.; // seedpiece C pool (goes into C_reserve)
	double C_seed_used = 0.; //seedpiece CHO used for growth at current time-step only
	double C_pool_used = 0.; //soluble CHO used for growth at current time-step only
	double C_content = 0.;
	double C_demand = 0.; //current time-step demand for all organs
	double C_demand_old = 0.; //demand from last time step
	double C_demand_nostress = 0.; //time-step demand for CHO for all organs if no water stress on haulm simulated
	double C_supply_nostress = 0.; //as belwo, except growth demand based on no water stress effect
	double C_supply = 0.; //actual amount of CHO able to support growth demand for current time-step - this is estimated as portion of C_pool (reset to 0 each time)
	double CH2O = 0.; // glucose, g
	double Tubfraction = 0.; //fraction of assimilate going to tubers from previous time-step
	double Tubmod = 0.; //value that modifies tuber sink strength
	double shootPart = 0., rootPart = 0., tubPart = 0.; //g CHO partitioned to root or shoot or tuber per h-1 per plant-1 for current time-step
	double shootPart_virtual = 0.; //potential g CHO which would have been allocated to shoot for time-step if nostress on leaf expansion occurred - amount is virtual because it is used as upper limit for root growth only in 2DSOIL
	double shootPart_old = 0., rootPart_old = 0., tubPart_old = 0.; //as above but at previous timestep
	// ***********************************

	// *** whole canopy area related information ***
	double leafArea = 0.;
	double LAI = 0.; //total LAI of green canopy
	double maxLAI = 0.; //maximum canopy lai reached at any point in simulation
	double greenLAI = 0.;
	double sunLAI = 0., shadeLAI = 0., deadLAI = 0.;
	double greenLeafArea = 0., deadLeafArea = 0., potentialLeafArea = 0.;
	double temperature = 0.; // leaf T, hourly
	double leafageEffect = 0.; //from Ng and Loomis, influence of weight leaf age in canopy on Pgross
	double reference = 0.; //from Ng and Loomis, reference Pg value at an LAI of 3 for use in leaf Rm calculation
	// *********************************************

	// *** whole plant nitrogen variables ***
	double totalReserveNitrogen = 0.; // THis is an artificial variable contaiing N supply not used in plant but taken from soil...
	double totalLeafNitrogen = 0.; // total nitrogen content in leaf in g plant-1
	double totalStemNitrogen = 0.; // total nitrogen content in stem in g plant-1
	double totalRootNitrogen = 0.; // total nitrogen content in root in g plant-1
	double totalTuberNitrogen = 0.; // total nitrogen content in tuber in g plant-1
	double totalDeadNitrogen = 0.; // total nitrogen lost from plant due to senescence, in g plant -1
	double NitrogenStressFactor = 0.; // whole plant N stress from 0 to 1 (NitrogenStressfactor from SIMGUI)
	double N_dead = 0.; // N content lost from plant in each time-step, g N (set to 0 at beginning of each iteration)
	double N_seed = 0.; // total N content in tuber seedpiece, g plant-1, 1% of dry mass of seed piece
	double N_reserve = 0.; // total N content available from seed to meet N demand; based on portion of C_seed used to support C demand
	double N_seedused = 0.; //Total N used from seedpiece for growth, g N plant-1
	double N_hourlyseedused = 0.; //N used from seedpiece for current time step only, g N plant-1 hour-1
	double N_hourlytranslocated = 0.; //N moved from buffers, below the optimum range only, for current time-step only, g N plant-1 hour-1
	double N_hourlygrowthdemand = 0.; //N required for new growth and to maintain organs at optimum level, g N plant-1 hour-1
	double EffectiveCanopySLA = 0.; //canopy SLA when C_pool is included as part of leafmass
	double NCRatioTotalPlant = 0.; // N to C ratio
	double NCRatioLeaves = 0.;
	double NCRatioStems = 0.;
	double NCRatioRoots = 0.;
	double NCRatioTubers = 0.;
	double newleafNitrogen = 0.; //accounting purposes, to remove N from pool of new leaf
	double newstemNitrogen = 0.; //accounting purposes, to remove N from pool of new stem


	// NOTE: MOVED ALL other N Variables that were based on SIMPOTATO Approach INTO NITROGENSTATUS STRUCTURE - public structure
	// **************************************

	// *** whole plant nitrogen variables - maizsim approach, not part of NitrogenStatus Structure ***
	double TotalNitrogen = 0.; //This is the total nitrogen content of the plant in g plant-1, from MAIZESIM
	double leaf_NFraction = 0.; //records the fraction of nitrogen in leaves YY, 0 to 1
	double leaf_N = 0.; //total nitrogen in the leaves of a plant YY, g N leaves-1 plant-1
	double leaf_N_content = 0.; //leaf nitrogen content (per unit square) of a plant YY, g N m-2 leaf area
	

	// * End of time-increment values for output to data files *
	double totalMass = 0., stemMass = 0., leafMass = 0., deadMass = 0., rootMass = 0., tuberMass = 0.; // this is redundant, but for convenience of access
	double mainstemMass = 0., basalstemMass = 0., apicalstemMass = 0., mainleafMass = 0., basalleafMass = 0., apicalleafMass = 0.;
	double mainleafArea = 0., basalleafArea = 0., apicalleafArea = 0.;
	int basalstemNo = 0, apicalstemNo = 0, mainleafNo = 0, basalleafNo = 0, apicalleafNo = 0;
	double Pgross = 0., Pnet = 0., cumulativeET = 0., Rg = 0., Rm = 0.; //cumulative mol or mmol co2 or  g h2o m-2 ground time d-1
	// *********************************************************

	// * Gas Exchange Data *
	double photosynthesis_gross = 0.; // gross photosynthesis, umolCO2 m-2 ground s-1
	double photosynthesis_grosssunlitleaf = 0.;//pg umol co2 m-2 leaf s-1, sunlit leaves
	double photosynthesis_grossshadedleaf = 0.;//pg umol co2 m-2 leaf s-1, shaded leaves
	double photosynthesis_net = 0.; // gross photosynthesis, umolCO2 m-2 ground s-1
	double photosynthesis_gross900 = 0.; //Pg at 850-950 umol co2 m-2 s-1
	double photosynthesis_net900 = 0.; // Pn at "
	double photosynthesis_grossleaf900 = 0.; //leaf p.s. (m2 is per unit leaf area)
	double photosynthesis_netsunlitleaf = 0.; //net p.s., umol co2 m-2 leaf s-1
	double photosynthesis_netshadedleaf = 0.; //net p.s., umol co2 m-2 leaf s-1
	double psistress_gs_factor = 0.; //effect of leaf water potential on stomatal conductance, 0 to 1, hourly
	double psistress_gs_factor_24h = 0.; //effect of leaf water potential on stomatal conductance, 0 to 1 averaged over 24h
	double Nstressfactorone_24h = 0.; //photosynthetic reduction factor, 0 to 1, 24h average
	double Nstressfactortwo_24h = 0.; // 24h root growth
	double psi_leafexpansion_stress_24h = 0.; //effect of leaf water potential (0 to 1) on leaf expansion rate 24h average
	double psi_leafexpansion_stress = 0.; //effect of leaf water potential (0 to 1) on leaf expansion rate per hour
	double transpiration_sunlitleaf = 0.; //leaf et, mmol H2O m-2 leaf s-1
	double transpiration_shadedleaf = 0.; //leaf et, mmol h2o m-2 leaf s-1
	double PhotosyntheticFeedback = 0.; //feedback on photosynthetic rates when Cpool gets too high - try to implement this onto TP? 
	double assimilate = 0.; //assimilation flux, g CO2 per plant per timestep
	double transpiration = 0.; //instanteous transpiration, mmol H2O m-2 ground s-1
	double instantResp = 0.; // instantaneous Rg and Rm, umol CO2 m-2 s-1
	double maintRespiration = 0.;// g CH20 for whole plant
	double maintRespleaf = 0., maintRespstem = 0., maintResptuber = 0., maintResproot = 0.; // maintenance respiration costs (g CH2O time-1) for time increment
	double growthRespiration = 0.; // g CH2O
	double Rgleaf = 0., Rgorgan = 0.; // gCH2O formed per g dry matter (specific growth respiration?)
	double sunlit_LAI, shaded_LAI; // sunlit and shaded LAI values
	double sunlit_PFD, shaded_PFD; // sunlit and shaded PFD values
	double sunPg = 0., shadePg = 0.; //store Pg (umol co2 m-2 leaf s-1) estimated for all sunlit and shaded leaves
	double conductance = 0.; //weighted stomatal conductance from gas exchange class for current hour, mmol H2o m-2 s-1
	double average_gs = 0.; //daily weight average for the above
	double internal_CO2 = 0.; //hourly leaf sub-stomatal CO2 concetration (ppm)
	double VPD = 0.; //as above
	// ************************************

	// * Growth values for individual organs *
	double dblPtubgro = 0.; // potential growth of all tubers per time step, g ch2o plant-1
	double dblPleafgro = 0., dblPmainleafgro = 0., dblPbasalleafgro = 0., dblPapicalleafgro = 0.; // as above for main, basal, and apical leaves
	double dblPstemgro = 0., dblProotgro = 0.;
	double dblTubgro = 0., dblLeafgro = 0., dblMainleafgro = 0., dblBasalleafgro = 0., dblApicalleafgro = 0.;
	double dblStemgro = 0., dblRootgro = 0.;
	double dblPyoungleafgro = 0., dblYoungleafgro = 0.; //potential dry mass gain of all leaves 30% of final physiological age
	double dblPoldleafgro = 0., dblOldleafgro = 0.;
	double dblnostressPmainleafgro = 0., dblnostressPleafgro = 0., dblnostressPstemgro = 0., dblnostressProotgro = 0., dblnostressPtubgro = 0.; // what haulm growth would be if there was no water stress restriction
	double dblnostressLeafgro = 0., dblnostressStemgro = 0., dblnostressRootgro = 0., dblnostressTubgro = 0.;
	//double leaf_NFraction; //records the fraction of nitrogen in leaves YY
	//double leaf_N; //total nitrogen in the leaves of a plant YY
	//double leaf_N_content; //leaf nitrogen content (per unit square) of a plant YY
	double carbonRatio = 0.; //Calcualte the ratio between carbon supply and carbon demand for leaf expasion YY

	double dblheatstressVegetative = 0.; //index from (0 most) (1 none) 0 to 1 indicating heat stress affect on either photosynthesis or area expansion
	double dblheatstressReproductive = 0.; // as above, in SPUDSIM, I just set these equal to each other since vegetative growth can continue post TI
	// *****************************************

	// * Developmental status info *
    int finalNodeNumber = 0; //final number of nodes
	int inodeNumber = 0; // currently initiated number of nodes or nodalunits
	int iyoungnodeNumber = 0; // number of currently initiated nodes that are considered young
	int ioldnodeNumber = 0; // # as above that are old but still growing
	int igreennodeNumber = 0; //number of nodes that have not senesced
	bool Bnewnode = 0; // true if time for new node
	double sowingDay = 0.;
	double age = 0.;
	// *****************************

	// * 24-h average values for SIMGUI whole plant responses *
	double Tdayave = 0., Tdaymax = 0., Tdaymin = 0., Parave = 0., Sradave = 0., Photoave = 0.; //par is in mol m-2 d-1; srad is MJ m-2 d-1 total radiation, photo is h
	double Tdaylag = 0.; //yesterday's 24 h ave T for use in tuber growth routines that were developed for 24h temperatures
	double dayinc = 0.; //current time of day, from 0 to 24 hr
	double Tcanopyave = 0.; //average 24h canopy T, C
	double average_LWP = 0.; // average 24hour leaf water potential, MPa
	double swdf1 = 0.; //average 24 hour water stress response, from 0 to 1 as in SIMGUI
	// *********************************************************
	
	// Miscellaneous - model testing
	double PFD = 0., sunPFD = 0., shadePFD = 0., SRAD = 0.;//instantaneous incident, shaded and sunlit PAR values in umol PAR m-2 s-1, SRAD in W m-2
	double Total_soil_water_uptake = 0.; //from 2dsoil water uptake, cumulative g H2O plant-1
	double Hourly_soil_water_uptake = 0.; //from 2dsoil wateruptake, hourly g H2O plant-1
	double Transpiration_supply_demand_ratio = 0.; //ratio of actual root uptake to spudsim canopy prediction transpiration, daily
	double Hourly_transpiration_supply_demand_ratio = 0.; // as above, but hourly basis
	


};
#endif

