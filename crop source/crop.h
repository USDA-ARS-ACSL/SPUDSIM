// The following ifdef block is the standard way of creating macros which make exporting 
// from a DLL simpler. All files within this DLL are compiled with the PLANT_EXPORTS
// symbol defined on the command line. this symbol should not be defined on any project
// that uses this DLL. This way any other project whose source files include this file see 
// PLANT_API functions as being imported from a DLL, whereas this DLL sees symbols
// defined with this macro as being exported.
// See creating a pocket PC dll using C++ article for information
// on how to do this

#ifdef _WIN32
#ifdef MYPLANT_EXPORTS
#define PLANT_API __declspec(dllexport)
#else
#define PLANT_API __declspec(dllexport)
#endif
#else
#define PLANT_API
#endif

#include "initInfo.h"

// Other global definitions go here

    double p_VMAX = 0,
		  Popare = 0,			//plant population per m2, i.e. plant density
		  PopSlab = 0,			//area of the slab versus the ground area per plant, m2 m-2; or, in layman's terms, it's the amount of a plant that fits over the slab
		  Trel = 0,            // relative time since emergence
		  Period = 0,          // one hour of time
		  //Emergence;       // Day of year plant emerged - I removed this, already defined in Development.CPP - DHF 3-15-07, but changed to Emerge below
		  Emerge = 0;
	TInitInfo	initInfo;
	double lwp_dawn=-0.05;	// predawn leaf water potential, MPa
	double pot_transpiration_old = 0; //potential transpiration from plant at prior time-step    
	double actual_transpiration_old = 0; //actual transpiration, or soil water uptake, to plant at prior time-step, g plant- d-1
	double old_shoot_weightPerM2 = 0; //declare and initialize a variable to save the above-ground biomass in each time step YY
		  

	bool  bEmergence=false; //boolean for emergence occurring
	int   ModNum=0, //Module Number
		  errPlant=1    // error number to return to 2DSOIL
          ;
	const double CO2=370, Press=98;
	const double cm2perm2=10000;
	double  TotalTillerLeafArea =0, TotalLeafArea =0, TillerLeafArea =0, 
		    MainStemLeafArea=0;
	double netPhotosyn =0,respiration=0, GPhotosyn= 0;
	int nNodes = 0;
	int nLeaves = 0;
	int ntLeaves = 0;
	int nTiller = 0;
	const double maxAge=47;
	double WaterUptake = 0; //hourly uptake from 2dsoil
	double NitrogenUptake = 0; //nitrogen uptake value from 2dsoil accumulated between time steps, mg / plant

	double TotalPotentialRootWaterUptake=0; //hourly potential extractable water from 2DSOIL g slab-1 d-1
	double NitrogenUptakeOld = 0; //N uptake from prior time-step, mg / plant
    double HourlyCarboUsed=0;
	double CumulativeNitrogenDemand=0.0; //grams plant-1
	double CumulativeActualNFromSoil=0.0; //grams plant-1
	double U_N = 0, U_M = 0, U_P = 0, U_D = 0;	//U_N maximum observed N uptake rate (g N m-2 ground-1) (Lindquist et al., 2007)
								//U_M maximum uptake rate as limited by max N fraction per unit (eqn 2 in above ref)
								//U_P potential rate of N accumulation (g N m-2 ground d-1) (as above)
								//U_D update rate (g N m-2 d-1) as limited by difference between potential and actual amount of N
									//in existing biomass (eqn 3
	double d = 0.075; //d; shape coefficient in logisit function to simulate cum N uptake, equation 9 above ref
	double q_n = 0.032; //q_n max ratio of daily N uptake measured daily growth rate (g N g-1;)
	double a = 4.10; //max N concentration for C4 species, 4.1%
	double b = 0.5; //shape coefficient in calculation of N concentraiton in relation to above ground biomass (eqn 4, above ref)
    
//Common Structures defined here
	double massIncrease = 0; //increase in shoot biomass, leaf + stem, in each time step
	double shoot_weightPerM2 = 0; // current shoot biomass?
	double SLNmin = 0.5; //base specific leaf N content
	double CurrentNUptakeError = 0;
	double CumulativeNUptakeError = 0;

    float LAMDAS = 0, LAMDAC = 0;
	const int NumNPD=4000, NumElD=3500, NumBPD=600, NSeepD = 2,
              NumSPD= 30, NumSD =10, NDrainD=2, NumDR=30,
			  NumGD = 3, NumPlD=100, 
              NMatD=15, NumModD=20, MBandD=15,
	          NumSurfDatD=3+NumGD+NumSD;
#pragma pack(2)
 struct ShootCommon{
    double PCRL,PCRQ,PCRS,HourlyCarboUsed,ET_demand,
	LCAI,Cover, Convr;
	float MaxRootDepth,Shade,Height,LAI,AWUPS, NitroDemand;
	float xBStem,yBStem,SGT,PSIM,
    LAREAT,PopRow,RowSp,RowAng,PopArea,COREC,
    ESCS,AWUPSS,SolRad,Total_Eor,
    Total_Pcrs,SIncrSink,Psild,
    OsmFac, EOMult,psil_, NDemandError, CumulativeNDemandError, 
	TotalRootWeight, InitialRootCarbo,ConstI[2], constK[2], Cmin0[2];
	int isGerminated, isEmerged;
	float TRWU_SIM, PSILT_SIM, PSISM_SIM;
	float Plantingdrymass,Plantingdepth,Plantinglength;
	int Bigleaf,Nstressoff,Wstressoff,Wstresstype;
  };


 //Weather
 struct WeathCommon{
	 int   MSW1,MSW2,MSW3,MSW4,MSW5,MSW6,MSW7;
	 float BSOLAR,ETCORR,BTEMP,ATEMP,ERAIN,BWIND,BIR,WINDA,IRAV;
	 int   JDAY, NCD, JDLAST;
     float CLDFAC,DEL[24],RINT[24],RNS,RNC,RAIN,IR;
	 float WIND,CO2,TDUSK,TDUSKY,CPREC[NumSD],TAIR[24],VPD[24],
			ROUGH,RADINT[24],WATTSM[24],DIFINT[24],ROWINC[24],
		 CLOUD,SHADOW[24],DIFWAT[24],DIRINT[24],WATACT,WATRAT,
		 WATPOT,RNLU;
	 int   NumF[40],NumFP;
	 float hFur[40],QF;
	 int   IFUR;
	 float GAIR[NumGD],PG,LATUDE,Longitude,Altitude,RI,PAR[24],PARINT[24],DAYLNG;
	 float AutoIrrigAmt;
	 int   AutoIrrigateF;
 };

 //grid
 struct GridCommon{
	     int     NumNP, NumEl, IJ, KAT, MBand,Nmat, KX[4][NumElD];
		 float   x[NumNPD], y[NumNPD], Area[NumElD], nodArea[NumNPD];
 };

//nodal
 struct NodeCommon{
	        int NumSol,NumG,ListN[NumNPD],ListNE[NumNPD], MatNumN[NumNPD];
			float hNew[NumNPD], ThNew[NumNPD], Vx[NumNPD], Vz[NumNPD],
				Q[NumNPD], Conc[NumSD][NumNPD], g[NumGD][NumNPD],
				Tmpr[NumNPD], Con[NumNPD], TcsXX[NumNPD], RO[NumNPD],
				hNew_org[NumNPD], QAct[NumNPD], ThetaAvailRZ, ThetaFullRZ;
			float ThAvail[NumNPD], ThFull[NMatD];
			bool lOrt;
			float QGas[NumGD][NumNPD], ThetaAir[NumNPD];
 };

//elements
 struct ElementCommon {

	 float Sink[NumNPD], cSink[NumSD][NumNPD],
		 gSink[NumGD][NumNPD], tSink[NumNPD],
		 RTWT[NumNPD],
		 RMassM[NumNPD], RDenM[NumNPD],
		 RMassY[NumNPD], RDenY[NumNPD];
	 int    MatNumE[NumElD];
	 float gSink_OM[NumGD][NumNPD], cSink_OM[NumGD][NumNPD],
		 gSink_root[NumGD][NumNPD], gSink_rootY[NumGD][NumNPD],
		 gSink_rootM[NumGD][NumNPD], gSink_N2O[NumGD][NumNPD];
 };




 //boundary
 struct  BoundaryCommon {
	 int  NumBP, NSurf, NVarBW, NVarBS, NVarBT, NVarBG,
		 NumSurfDat, NSeep, NSP[NSeepD], NP[NumSPD][NSeepD],
		 NDrain, NDR[NDrainD], ND[NumDR][NDrainD],
		 KXB[NumBPD];
	 int CodeW[NumNPD], CodeS[NumNPD], CodeT[NumNPD],
		 CodeG[NumNPD], PCodeW[NumNPD];
	 float Width[NumBPD], VarBW[3][NumBPD],
		 VarBS[NumSD][NumBPD],
		 VarBT[4][NumBPD], VarBG[3][NumGD][NumBPD], EO, Tpot, gridWidth;
	 int PondingbyHead, PondingbyFlux;
 };

// Time
  struct TimeCommon{
	         double tNext[NumModD],dtMx[4],Time, Step, dtOpt,dtMin, 
				     dMul1, dMul2,tTDB[4],Tfin,tAtm;
			 float  Tinit;
			 int    lInput,Iter;
			 int   DailyOutput, HourlyOutput,RunFlag, DailyWeather, HourlyWeather;
			 int   beginDay, sowingDay, emergeDay, endDay, OutputSoilNo, OutPutSoilYes, Year;
		     int iTime,iDawn,iDusk;
			 double TimeStep;
 };

 //modules
 struct   ModuleCommon{
	              int NumMod,Movers[4], NShoot;
 };
 struct   ErrorCommon{
	              int errPlant;
 };

 struct FileCommon{
	 double starter;
	 char WeatherFile[256], TimeFile[256], BiologyFile[256],
		 ClimateFile[256], NitrogenFile[256], SoluteFile[256],
		 ParamGasFile[256], SoilFile[256],
		 ManagementFile[256], IrrigationFile[256], DripFile[256],
		 WaterFile[256], WaterBoundaryFile[256],
		 GraphicsFile[256], InitialsFile[256], VarietyFile[256],
		 NodeGraphics[256], ElemGraphics[256], NodeGeomFile[256],
		 GeometryFile[256], SurfaceGraphics[256],
		 FluxGraphics[256], MasssBalanceFile[256],
		 MassBalanceFileOut[256], LeafFileIn[256],
		 OrganicMatterGraphics[256],
		 RunFile[256], MassBalanceRunoffFileOut[256],
		 MulchFile[256], MassBalanceMulchFileOut[256];
 };

#pragma pack()


 
#ifdef __cplusplus
extern "C" {
#endif

// Your exported function headers go here
#ifdef _WIN32
PLANT_API void _stdcall CROP(struct ShootCommon *, WeathCommon *,
#else
PLANT_API void crop(struct ShootCommon    *, WeathCommon    *,
#endif
							         GridCommon     *, NodeCommon     *,
									 ElementCommon  *, BoundaryCommon *,
									 TimeCommon     *, ModuleCommon   *,
									 FileCommon     *);
#ifdef __cplusplus
}
#endif


 


