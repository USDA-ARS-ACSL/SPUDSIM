#include "StdAfx.h"
#include <cmath>
#include "development.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <stdlib.h>

//#using <mscorlib.dll>

using namespace std;


CDevelopment::CDevelopment(void)
{
	T_base = 0.0;  T_opt = 27.98; T_ceil = 40.46;
	dblCdtt17 = dblCdtt22 = dblCstt = 0;
	dblTind = dblPartub = 0;
	dblDAS = 0;
	dblLAR = 0;
	dblBAR = 0;
	iDAS = 0;
	setParms();

	//LvsInitiated = LvsAppeared = LvsExpanded = Anthesis=0;
    //GerminationRate = EmergenceRate = LvsInitiated = LvsAppeared = LvsExpanded = Anthesis =0;
	//GDDsum = GDDgrain = dGDD = P2 = phyllochronsFromTI = 0;
	//GTIsum = dGTI = 0;
	//GDD_rating = 1000;
	//Rmax_LIR = Rmax_LTAR = Rmax_Germination = Rmax_Emergence =0;
	//totLeafNo = juvLeafNo = 15;
	//initLeafNo =  finalLeafNo = 5;
	//curLeafNo =1; 
	//LvsAtTI = 1;
}

CDevelopment::CDevelopment(const TInitInfo& info)
{
	
	//LvsInitiated = LvsAppeared = LvsExpanded = Anthesis=0;
    //GerminationRate = EmergenceRate = LvsInitiated = LvsAppeared = LvsExpanded = Anthesis =0;
	//GDDsum = GDDgrain = dGDD= P2 = phyllochronsFromTI =0;
	//GTIsum = dGTI = 0;
	//GDD_rating = 1000;
	//Rmax_LIR = Rmax_LTAR = Rmax_Germination = Rmax_Emergence =0;
	//T_base = 0.0;  T_opt = 30.0; T_ceil = 40.0;
	//totLeafNo = juvLeafNo = 15;
	//initLeafNo =  finalLeafNo = 5;
	//curLeafNo =1; 
	//LvsAtTI = 1;
	dblCdtt17 = dblCdtt22 = dblCstt = 0;
	dblTind = dblPartub = 0;
	dblDAS = 0;
	dblLAR = 0;
	dblBAR = 0;
	iDAS = 0;
	initInfo = info;
	setParms();
}

CDevelopment::~CDevelopment(void)
{
}

void CDevelopment::response(int iCur, const TInitInfo info, const TWeather& weather, double Tdaymax, double Tdaymin, double Tdayave)
{
	/******************************************************/
	/* Legacy whole plant temperature response functions come from SIMGUI (modified for C++), 24 hour basis
	/* That is, they use prior 24 day T information only and have not been modified as yet */
	/*******************************************************/
	//if (fmod(iCur,((24*60)/initInfo.timeStep))==0)//runs at beginnning of day
	if ((fmod(iCur+1,(24*60)/initInfo.timeStep) < 0.001)||(iCur==0))//run this at end of each day, since estimates of par are set to zero at beginning of each day
	{
	//if (abs(weather.time - 1) <= info.timeStep/60/24){
 		dblDtt17 = fdtt17(Tdayave, Tdayave, Tdaymax, Tdaymin);
		dblDtt22 = fdtt22(Tdayave, Tdayave, Tdaymax, Tdaymin);
		dblStt = fstt(Tdayave, Tdayave, Tdaymax, Tdaymin);
		dblCdtt17 = dblCdtt17 + dblDtt17;
		dblCdtt22 = dblCdtt22 + dblDtt22;
		dblCstt = dblCstt + dblStt;
	}
}
//* The three f... functions were SIMGUI subroutine therm.  Original comments were removed see original SIMPOTATO for documentation
double CDevelopment::fdtt17(double dT, double tempm, double tempmx, double tempmn)
/*************************************************************/
/* Legacy thermal time based resposne with 17C optimum from SIMPOTATO
/* See SIMPOTATO code for original derviations */
/***********************************************/
{
	double TM1_17 = 0.0, TM2_17 = 0.0, TM3_17 = 0.0, TM4_17 = 0.0, TM5_17 = 0.0, topcnt = 0.0, trf = 0.0, ttmp = 0.0, TX1 = 0.0, TX2 = 0.0;
	double dtt; 
	double tmfac[9];
	int i = 0, Ihead = 0;
//C   INITIALIZE COEFFICIENTS FOR THERMAL TIME CALCULATIONS
	Ihead = 0;
	if (Ihead==0)
	{
//C   For Tuber Initiation Use Optimum Temperature of 17 C
		Ihead = 1;
		TM1_17=5;
		TM2_17=14;
		TM3_17=17;
		TM4_17=20;
		TM5_17=35;
		TX1=33;
		TX2=50;
	}
	dtt = 0;
	tempm=0;
	if (iStage <= 4) {
		//C   ACCUMULATE THERMAL TIME OVER 8 INTERPOLATED TEMPERATURE VALUES BETWEEN
		//C        TEMPMX & TEMPMN
        for (i=1; i<9;i++)
		{
			/* DHF - TMFAC INITIALIZED AND KEPT AS LOCAL VARIABLE IN C++ VERSION*/

			double x(i);
			tmfac[i] = 0.931+0.114*i-0.0703*pow(x,2) + 0.0053*pow(x,3);
			ttmp=tempmn+tmfac[i]*(tempmx-tempmn);
			//C **********************************************************************
			//C   NEW SIMPOTATO CODE -- TOM HODGES, FEB 1992                         *
			//C   THIS IS ON A 0. - 1.2 SCALE SO MULTIPLY BY 5/6 TO PUT ON 0-1 SCALE *
			//C **********************************************************************
			//C   DTT17 based on optimum temperature of tm3_17=17 C
			if ((ttmp > TM1_17) && (ttmp < TM2_17))
				topcnt = (ttmp - TM1_17) / (TM2_17 - TM1_17);
			else if (ttmp < TM3_17) 
				topcnt = 1. + (ttmp - TM2_17)*0.2/(TM3_17-TM2_17);
			else if (ttmp < TM4_17)
				topcnt = 1.2 - (ttmp-TM3_17)*0.2/(TM4_17-TM3_17);
			else if (ttmp < TM5_17)
				topcnt = 1. / (TM5_17 - TM4_17) * (TM5_17 - ttmp);
			else if ((ttmp <= TM1_17) || (ttmp >= TM5_17))
				topcnt = 0.0;
			//	C   REDUCE TO A 0. -> 1. SCALE
            topcnt=topcnt*5./6.;
			if (topcnt < 0.0)
               topcnt=0.0;
			else if (topcnt < 1.0)
               topcnt=1.0;

            dtt=dtt+topcnt/8.;
            tempm=ttmp/8.+tempm;
		}

		//C   EXTREMELY HIGH TEMPMX REDUCES THERMAL TIME ACCUMULATION
		//C     USE ABOUT 28 TO 33 C FOR TX1, AND 50 C FOR TX2
		if ((tempmx > TX1) && (tempmx < TX2)) {
            dtt=dtt*(1.-(tempmx-TX1)/(TX2-TX1));
		}
		else if (tempmx>=TX2) {
            dtt=0.;
		}
		//C    MULTIPLY dtt17 and dtt22 BY FACTOR BASED ON DAILY TEMPERATURE RANGE
		//C      AT <=10 C   " =1.
		//C        "  20 C   " =1.33
		//C        "  25 C   " =1.5
        trf=1.;
        if ((tempmx-tempmn) > 25) trf=1.5;
		if (((tempmx-tempmn) <= 25) & ((tempmx-tempmn) > 10.)) {
			trf=0.67 + 0.83*(tempmx-tempmn)/25.;
		}
		dtt = dtt*trf;
	}
	return dtt;
}

double CDevelopment::fdtt22(double dT, double tempm, double tempmx, double tempmn)
{
/*************************************************************/
/* Legacy thermal time based resposne with 22C optimum from SIMPOTATO
/* See SIMPOTATO code for original derviations */
/***********************************************/
//C  local variables
	double TM1_22 = 0.0, TM2_22 = 0.0, TM3_22 = 0.0, TM4_22 = 0.0, TM5_22 = 0.0, TOP_22 = 0.0, trf = 0.0, ttmp = 0.0, TX1 = 0.0, TX2 = 0.0;
	double tmfac[9], dtt22 = 0.0;
	int i = 0, Ihead = 0;
	Ihead = 0;
	if (Ihead==0)
	{
		TX1=33;
		TX2=50;
		TM1_22=6;
		TM2_22=17.5;
		TM3_22=22;
		TM4_22=26.5;
		TM5_22=36;
	}
	dtt22=0;
	tempm=0;
	if (iStage <= 4) {
        for (i=1; i<9;i++)
		{
			/* DHF - TMFAC INITIALIZED AND KEPT AS LOCAL VARIABLE IN C++ VERSION*/
			double x(i);
			tmfac[i] = 0.931+0.114*i-0.0703*pow(x,2) + 0.0053*pow(x,3);
			ttmp=tempmn+tmfac[i]*(tempmx-tempmn);
			//C  dtt22 based on optimum temperature of TM3_22=22 C
			if ((ttmp < TM1_22) && (ttmp < TM2_22)) 
               TOP_22=(ttmp-TM1_22)/(TM2_22-TM1_22);
			else if (ttmp < TM3_22) 
               TOP_22=1.+(ttmp-TM2_22)*.2/(TM3_22-TM2_22);
			else if (ttmp < TM4_22) 
               TOP_22=1.2-(ttmp-TM3_22)*.2/(TM4_22-TM3_22);
			else if (ttmp < TM5_22)
               TOP_22=1./(TM5_22-TM4_22)*(TM5_22-ttmp);
			else if ((ttmp <= TM1_22)||(ttmp >= TM5_22))
               TOP_22=0.0;
			//	C   REDUCE TO A 0. -> 1. SCALE
            TOP_22=TOP_22*5./6.;
            if (TOP_22 < 0.0)
               TOP_22=0.0;
			else if (TOP_22 > 1.0)
               TOP_22=1.0;
            dtt22=dtt22+TOP_22/8.;
			tempm=ttmp/8.+tempm;
		}
		if ((tempmx > TX1) && (tempmx < TX2)) {
            dtt22=dtt22*(1.-(tempmx-TX1)/(TX2-TX1));
		}
		else if (tempmx>=TX2) {
            dtt22=0.;
		}
		trf=1.;
		if ((tempmx-tempmn) > 25) trf=1.5;
		if (((tempmx-tempmn) <= 25) & ((tempmx-tempmn) > 10.)) {
			trf=0.67 + 0.83*(tempmx-tempmn)/25.;

		}
		dtt22 = dtt22*trf;
	}
	return dtt22;
}	

double CDevelopment::fstt(double dT, double tempm, double tempmx, double tempmn)
{
	/*************************************************************/
	/* Legacy soil thermal time based resposne with 22C optimum from SIMPOTATO
	/* See SIMPOTATO code for original derviations */
	/***********************************************/
	//C  local variables
	double ttmp = 0.0, tubcnt = 0.0;
	double tmfac[9];
	double stt = 0.0;
	stt = 0.0;
	int i = 0;
	tempm=0.0;
	if (iStage <= 4)	{
        for (i=1; i<9;i++)
		{
			double x(i);
			tmfac[i] = 0.931+0.114*i-0.0703*pow(x,2) + 0.0053*pow(x,3);
			ttmp=tempmn+tmfac[i]*(tempmx-tempmn);
			tubcnt=1.25-.0962*fabs(ttmp-15.);
			if (tubcnt < 0.0) 
				tubcnt = 0.;
			else if (tubcnt > 1)
				tubcnt =1.;
			stt = stt + tubcnt / 8.;
			tempm = ttmp/8. + tempm;
		}
	}
	else{
		stt = 1 - 0.058823*fabs(dT-22); //*should be fabs(ST(1)-22) once soil layer temperature is added DHF, Jan 05
		if (stt < 0) stt = 0.;
		if (stt > 1) stt = 1.;
	}
	return stt;
}


void CDevelopment::setParms() // dt in days
{
	/******************************/
	/* Determine initial status of plant (germinated, emerged, etc.)
	/* Ascribe phenological stage and initial canopy architecture
	/* As of 2008, user input files MUST start with the plant at emergence
	/* -eventually plan on using soil moisture and temperature to help with germination, but not implemented yet, DHF*/
	/*******************************/

	double dt = initInfo.timeStep/(24*60); //converting minute to day decimal, 1= a day
	dblSdepth = initInfo.Plantingdepth;
	dblSprlthp = initInfo.Sproutlength;
	iMatjd = int(initInfo.endDay);
	//Determine initial developmental stage based on planting and emergence date (if provided)
	//Note that with current SPUDSIM versions, plant is assumed to have emerged
	// In SIMPOTATO, there were 3 cases: at emergence (stage 7), at germination (stage 6), or pregermination (stage 5)
	// With SPUDSIM, no simulation of plant processes will occur until plant has emerged (istage 7)
	// Thus, cases depend on inputs for beginday, sowingday, and emergenceday
	// By default, nothing should happen until emergence day is reached
	if (initInfo.beginDay < initInfo.sowingDay)
	{
		preplanting.day = initInfo.beginDay;
		preplanting.done = true;
		iStage = 5;
		dblSprlth = 0.;
	}
	if ((initInfo.beginDay == initInfo.sowingDay) && (initInfo.beginDay != initInfo.emergeDay))
	{
		preplanting.day = initInfo.beginDay;
		preplanting.done = true;
		germination.day = initInfo.sowingDay;
		germination.done = true;
		iStage = 6;
		dblSprlth = 0.;
	}
	if (initInfo.beginDay == initInfo.emergeDay)
	{
		preplanting.day = initInfo.beginDay;
		preplanting.done = true;
		germination.day = initInfo.sowingDay;
		germination.done = true;
		emergence.day = initInfo.emergeDay;
		emergence.done = true;
		iStage = 7;
		dblSprlth = initInfo.Plantingdepth; //hardcoded for now
	}
	
	T_base = 0.0;  T_opt = 27.2; T_ceil = 39.5; Rmax_leaf = 0.96*dt;  //from Fleisher et al (2005) SPAR data
	Rmax_stem = 0.5 * dt; //from Ng and Loomis
	initLeafNo = 4; //nodal units to be instantiated upon sprout emergence
	LvsInitiated = LvsAppeared = 0.;
}

void CDevelopment::tubinit(int iCur, const TInitInfo info, const TWeather& weather, const TNitrogenStatus NitrogenStatus, double area, double Tdaymax, double Tdaymin, double Tdayave, double Sradave, double Photoave)
{
	/**************************************************/
	/* Routine is legacy code ported to C++ from SIMPOTATO
	/* Only runs at end of time-step since tuber initiation based on prior 24h min, max T, solar radiation, etc
	/* See SIMPOTATO for additional information on routine
	/* NOTE:
	/* Need to incorporate nitrogen factor*/
	/****************************************/


	//double ampfac, areafac, dlfac, nitfac, solfac, tday, tmn, tempfac;
	int count = 0;

	//if (!fmod(iCur,((24*60  )/initInfo.timeStep))==0) return; 
	if ((Sradave == 0)||(Photoave == 0)) return;
	if (!(fmod(iCur+1,((24*60)/initInfo.timeStep)) < 0.001)) return;//run this at end of each time-step, since estimates of par are set to zero at beginning of each day
	//if (!(abs(weather.time - 1) <= info.timeStep/60/24)) return;
	if ((iStage < 5) && (dblTind < 300)) {
		//if (dblTind < 225 ) dblTind = dblTind * 0.85; //3/3/2017 - DHF this reduction in Tind plus Areafac were making calibration nearly impossible - commented out works much better
		//if (dblTind >= 225) dblTind = dblTind * 0.97; // as above
		tmn=(Tdaymax + Tdaymin)*0.5;
		if (tmn == 15) 
			tempfac = info.A1;
		else if ((tmn >= 5) && (tmn <= 25))
			tempfac = info.A1 - fabs(15.-tmn)*info.A2;
		else if ((tmn < 5) && (tmn > 0))
			tempfac = 1. - (5. - tmn) / 5.;
		else if ((tmn > 25) && (tmn < 40.))
			tempfac = 1. + (25. - tmn) / 15.;
		else if ((tmn >= 40.) || (tmn <= 0.))
			tempfac = 0.;
		if ((area*info.plantDensity/10000.) > 0) areafac = pow(area/3.,0.5)/info.A3;
		if (areafac < 1) areafac = 1.;
		solfac = 1. + pow((Sradave/Photoave),0.5)/info.A4;
		//solfac = 1;
		if ((Tdaymax-Tdaymin) < 25.)
			ampfac = 1 + info.A5 * (Tdaymax - Tdaymin)/25.;
		else
			ampfac = info.A6;
		if (Photoave > 12) 
			dlfac = 1 + (18. - Photoave)*info.A7/6.;
		else
			dlfac = info.A8;
		//if (Photoave > 18) dlfac = 1;  //Original SPUDSIM since April 2019
		if (Photoave > 16) dlfac = info.A7;
		//if (Photoave > 16.5) dlfac = info.A7;
		double Noptimum = (NitrogenStatus.optimumleafNitrogenRatio + NitrogenStatus.optimumstemNitrogenRatio) / 2.;
		double Nminimum = (NitrogenStatus.minimumleafNitrogenRatio + NitrogenStatus.minimumstemNitrogenRatio) / 2.;
		double Nactual = (NitrogenStatus.leafNitrogenAmount + NitrogenStatus.stemNitrogenAmount) / (1./NitrogenStatus.actualleafNitrogenRatio*NitrogenStatus.leafNitrogenAmount + 1./NitrogenStatus.actualstemNitrogenRatio*NitrogenStatus.stemNitrogenAmount);
		if (Nactual >= (Noptimum + (Noptimum - Nminimum)))
			nitfac = info.A9;
		else if (Nactual > Noptimum)
			nitfac = 1. - (1.- info.A9) * (Nactual-Noptimum)/(Noptimum-Nminimum);
		else if (Nactual == Noptimum)
			nitfac = 1.;
		else if (Nactual <= Nminimum)
			nitfac = info.A10;
		else if (Nactual  < Noptimum)
			nitfac = 1. + (info.A10-1.)*(Noptimum-Nactual)/(Noptimum-Nminimum);
		if (info.Nitrogen_stress_off==1) nitfac =1.;
		//nitfac = 1; //*****2009PAPER
		//ampfac = info.A6; //* fix for growth chamber data - should remove Tamplitude from TI index


		//cout << setiosflags(ios::left) << setw(4) << tempfac << setw(4) << areafac << setw(4) << solfac << setw(4) << ampfac << setw(4) << dlfac setw(4) << nitfac << endl;
		tday = tempfac*areafac*solfac*ampfac*dlfac*nitfac;
		if (count == 0) {
			
		}
		if (dblTind < 225) {
			dblTind = dblTind + tday;
			/*
			ofstream tindout("c:\\SPUDSIM v2-0\\testsim\\2017 research\\calibration\\phenology\\tind1.dat",ios_base::app);
			if (tindout) {
				tindout << setw(6) << setprecision(3) << tmn
					<< setw(6) << setprecision(3) << Tdaymax
					<< setw(6) << setprecision(3) << Tdaymin
					<< setw(6) << setprecision(3) << Photoave
					<< setw(6) << setprecision(3) << Sradave
					<< setw(6) << setprecision(3) << tempfac
					<< setw(6) << setprecision(3) << solfac
					<< setw(6) << setprecision(3) << ampfac
					<< setw(6) << setprecision(3) << dlfac
					<< setw(6) << setprecision(3) << nitfac
					<< setw(6) << setprecision(3) << tday
					<< setw(6) << setprecision(3) << dblTind
					<< endl;
				//tindout.close();
			}
			*/
			if (dblTind > 235) dblTind = 235.;		


			//cout << tempfac << " " << solfac << " " << areafac << " " << ampfac << " " << dlfac << " " << nitfac << " " << tday << " " << dblTind << endl;

		}
		else

			//dblTind = dblTind + tday / 6; //DHF - original denominator here almost forces crop to never get out of early growth stage, plant doesn't age unless high coefficients are used, really difficult to calibrate
			dblTind = dblTind + tday / 2.;

		if ((iStage <= 4) && (dblTind >= 225.)) {
			dblPartub = (dblTind - 225.)/75.;
			if (dblPartub >= 1) dblPartub = 1.;
		}
		else
			dblPartub = 0.;
	}
}

void CDevelopment::update(int iCur, const TWeather& weather, const TInitInfo info, const TNitrogenStatus NitrogenStatus, double area, double greenLAI, double maxLAI, double dayinc, double Tdaymax, double Tdaymin, double Tdayave, double Tdaylag, double Sradave, double Photoave, double Cseed, double Cpool)
{
	/***************************************/
	/* Track phenological progress including tuber initiation
	/* Obtain temperature response functions related to whole plant development
	/* Estimate potential leaf and branch appearance rates for current time-step
	/* Estimate tuber sink-strength based on phenology, temperature
	/* Assign new ontological information when new stages are reached
	/* Notes:
	/*		-most routines are predominantly legacy recoded versions from SIMPOTATO
	/*		-model is currently restricted to start at emergence - therefore, ISTAGE will always initially 
				be equal to 7 as set in the SetParms method 
	/* Nitrogen needs:
			-Need N restriction on leaf and branch appearance rate computation
	/***************************************/

	double Jday = weather.jday;
	double daylen = weather.dayLength;
	double T_air = weather.airT;
	double dt = info.timeStep/(24.*60.); //converting minute to day decimal, 1= a day
	double cumdep;
	int flag = 0;
	double valve = 0.2 / (60./info.timeStep);
	
	if ((Jday == info.emergeDay) && (emergence.done != true))
	{
		iStage = 7;
		emergence.done = true;
		emergence.day = Jday;
	}

	if (iStage <= 7) //compute whole plant thermal time response using simgui 24h responses (routines will only run at end of day)
	{	
		response(iCur, info, weather, Tdaymax, Tdaymin, Tdayave);
		tubinit(iCur, info, weather, NitrogenStatus, area, Tdaymax, Tdaymin, Tdayave, Sradave, Photoave);
	}

	if (iStage == 5){ // Determine planting date
		if (Jday != info.sowingDay) return;
		iDAS = 0;
		dblDAS = 0.;
		dblTind = 0.;
		dblSprlth = 0.;
		iPhase = 1;
		dblXstage = 0.;
		flag = 0;
		cumdep = 0.;
		dblXstage = 0.;
		//for (int L = 1; intNlayr; L++)
		//{
		//	if (flag == 0) {
		//		cumdep = cumdep+dblDlayr[ L ];
		//		intL0 = intL0 + 1;
		//	}
		//	if (dblSdepth < cumdep) flag = 1;
		//}
		if (dblSdepth < info.Plantingdepth) flag =1; // jump out of routine if error
		iJinit = 0;
	}
	else if (iStage == 6){ // Determine germination date
		//TO DO: Germination based on soilmoisture.  Will need to link with 2dSOIL...
		// original simgui routine is commented out below
		//if (dblSw[ intL0 ] <= dblLl[ intL0 ]){
		//	swsd = (dblSw[intL0]-dblLl[intL0])*0.65+(dblSw[intL0+1]-dblLl[intL0+1])*0.35;
		//	iDas= iDas + 1;
		if (Jday == info.emergeDay) {
			iPhase = 1;
		}
			if ((dblDAS >= 30) && (Jday != info.emergeDay)) {
				iPhase = 1; //DHF - force emergence
				//iStage = 7; //
				//dblDplants = 0;
				//write(STOUT,1000); need to fix this
				//write(41,1000);
				return;
			}
		//	if (swsd < 0.02) return;
		//}
		if ((dblCstt < 7.35) && (dblSprlthp == 0) && (Jday != info.emergeDay)) return;
		//iPhase = 1;
	}
	else if (iStage == 7){ // seedling emergence - will move parts of this routine to growth routine at later point, for now commented out
		//actually, growth prior to this stage should already be converd in stage 6, germination
		dblDAS = dblDAS + dt;
		iDAS = int(dblDAS);
		//if (dblSplthp == 0){
		//	dblGrospr = 17 * dblStt;
		//}
		//else {
		//	dblSpgrof = (0.167*(1-0.82* exp(-0.65*dblSplthp)))/0.0344;
		//	dblSGrospr = 2.84 * dblStt * dblSpgrof;
		//}
		//dblSprlth = dblSprlth + dblGrospr;
		//sprwt = dblGrospr*0.0272*dblStt;
		//dblTsprwt=sprwt+dblTsprwt;
		//dblGrort = sprwt;
		//dblSeedrv = dblSeedrv - (dblGrort+sprwt);
		//dblRtwt = dblRtwt + dblGrort;
		//if (intEmerge == 0) {
	//		if (dblSprlth/1000 < dblSdepth) return;
	//	}
	//	else
	//		if (intDay != intEmerge) return;
		iPhase = 1;
		dblXstage = 1.;
//		dblTempm = 6; NEED TO FIX THIS! and add to structure Weather variable list?
	}
	else if (iStage == 1) { //Vegetative growth, post emergence
		dblDAS = dblDAS + dt;
		iDAS = int(dblDAS);
		if (fmod(iCur,((24.*60.)/initInfo.timeStep))==0) dblXstage = 1 + dblTind/225.;

		//dblLAR = beta_fn(T_air, Rmax_leaf, T_opt, T_ceil);
		dblLAR = beta_fn(Tdaylag, Rmax_leaf, T_opt, T_ceil);
		//if ((Cpool+Cseed)*valve > 0.0525) dblBAR = beta_fn(Tdaylag, Rmax_stem, T_opt, T_ceil); //basal branch initiation
		if ((Cpool*valve)>0.0525) dblBAR = beta_fn(Tdaylag, Rmax_stem, T_opt, T_ceil);
// NEED TO ADD N RESTRICTION HERE AS WELL
		if ((Jday >= iMatjd)||(area <= 0.0000001) || ((Jday < info.sowingDay)&&((Jday+365)>=iMatjd)) || (greenLAI < maxLAI*0.1)){
			//skip to next stage if canopy severely senesced or other developmental problem - DHF and SIMGUI
			iStage = 4;
			iRet = 1;
			iPhase = 1;
		}
		if (dblTind < 225) return; //tubers initiate when TIND exceeds 225
		iJinit = int(dblDAS);
		iPhase = 1;

	}
	else if (iStage == 2){ //determine beginning of linear tuber growth
		if (fmod(iCur,((24.*60.)/initInfo.timeStep))==0) dblXstage = 2. + (dblTind-225.)/75.;
		dblDAS = dblDAS + dt;
		iDAS = int(dblDAS);
		if ((Cpool*valve) > 0.003) dblLAR = beta_fn(Tdaylag, Rmax_leaf, T_opt, T_ceil);
		if ((Cpool*valve) > 0.0525) dblBAR = beta_fn(Tdaylag, Rmax_stem, T_opt, T_ceil);
		if (greenLAI < 0.4 * maxLAI) { //Assume no new vegetative growth points once canopy has senesced 60% of maximum
			dblLAR = 0. ; dblBAR = 0.;
		}

// NEED TO ADD N RESTRICTION HERE AS WELL
		if ((Jday >= iMatjd) || (area <= 0.0000001) || ((Jday < info.sowingDay) && ((Jday+365) >= iMatjd)) || (greenLAI < maxLAI*0.1)){
			//skip to next stage if canopy severely senesced or other developmental problem - DHF and SIMGUI
			iStage = 3;
			iRet = 1;
			iPhase = 1;
		}
		if (dblTind < 300) return;
		iPhase = 1;
	}
	else if (iStage == 3) { //dominatn tuber growth
		if (fmod(iCur,((24.*60.)/initInfo.timeStep))==0)
		{
			//dblXstage = dblXstage + dblDtt17/9;
			//dblXstage = dblXstage + dblDtt17/8;//life-span too long under cool temperatures with original SIMGUI values
			dblXstage = dblXstage + dblDtt17 / 6.;//spudsimv2-0 change, still too-long under most conditions, couples with senescence rate
		}
		dblDAS = dblDAS +dt;
		iDAS = int(dblDAS);
		if ((Cpool*valve) > 0.003) dblLAR = beta_fn(Tdaylag, Rmax_leaf, T_opt, T_ceil);
		if ((Cpool*valve) > 0.0525) dblBAR = beta_fn(Tdaylag, Rmax_stem, T_opt, T_ceil);
		if (greenLAI < 0.4 * maxLAI) { //Assume no new vegetative growth points once canopy has senesced 60% of maximum
			dblLAR = 0. ; dblBAR = 0.;
		}

// NEED TO ADD N RESTRICTION HERE AS WELL
		if ((Jday >= iMatjd) || (area <= 0.0000001) || ((Jday < info.sowingDay) && ((Jday+365) >= iMatjd))){
			iRet = 1;
			iPhase = 1;
		}
	}
	else if (iStage == 4) { // end of tuber growth
		iRet = 1;
		iPhase = 1;
	}
	else {
		//write (STOUT,*)' ISTAGE= ', ISTAGE;
		//STOP 888;
	}

// This is  PHASEI
	if (iPhase == 1){
		if (iStage == 1){
			iStage = 2;
			iPhase = 0;
			tuberlinear.day = dblDAS;
			tuberlinear.done = true;
			return;
		}
		else if (iStage == 2) {
			iStage = 3;
			iPhase =0;
			tuberdominant.day = dblDAS;
			tuberdominant.done = true;
			return;
		}
		else if (iStage == 3){
			iStage = 4;
			iPhase =0;
			maturity.day = dblDAS;
			maturity.done = true;
			return;
		}
		else if (iStage == 4){ //reset all variables for next simulation?
			iStage = 5;
			iPhase =0;
			//for (i=1; intNlayr; i++){
				//RLV(I)=0.0
				//RWU(I)=0.0
			//}
			return;
		}
		else if  (iStage == 5){ //move to pre-germination phase
			iStage =6;
			iPhase =0;
			germination.day = dblDAS;
			germination.done = true;
			//numnode=0.0
			//numleaf=0.0
			//STEMLEN=0.0
			dblCdtt17=0.;
			dblCdtt22=0.;
			dblCstt=0.;
			dblTsprwt=0.;
			//for (i=1; intNlayr; i++){
				//RLV(I)=0.0
				//RWU(I)=0.0
			//}
			//if (dblSprlthp > 0)
			//	iStage = 7;
			//else
				return;
		}
		else if (iStage == 6){
			iStage = 7;
			emergence.day = dblDAS;
			emergence.done = true;
			iPhase =0;
			dblCdtt17=0.;
			dblCdtt22=0.;
			return;
		}
		else if (iStage == 7){ // set initial values after emergence
			iStage = 1;
			vegetative.day = dblDAS;
			vegetative.done = true;
			iPhase =0;
			if (info.bigleaf == 0) 
			{
				dblCdtt17 = 0.; dblCdtt22 = 0.;dblCstt = 0.;
			}
			if (info.bigleaf == 1)
			{
				dblDtt17=0.; dblDtt22=1.; dblStt = 1.;dblCdtt17=1.; dblCdtt22 = 1.; dblCstt=1.;
			}
			LvsAppeared = 4.; //set 4 initial leaves at emergence along mainstem
			dblSprlth = initInfo.Plantingdepth; //hardcoded for now
			return;
		}
	}
	return;
}

double CDevelopment::beta_fn(double t, double R_max, double t_m, double t_e)
{
/*******************************************/
//Generalized Temperature Response Model based on beta function - Fleisher et al., 2006
//Used for potential leaf appearance and basal branch appearance rates
/********************************/
	double deltaT1, deltaT2, alpha;
	const double t_b = 0;
	deltaT1 = t_m - t_b;
	deltaT2 = t_e - t_m;
	alpha = (deltaT1/deltaT2);
	if (deltaT1 < 0. || deltaT2 < 0. ) return -1;
	if (t < 0) return 0; //avoid situation where negative temperature results in undefined result - physiologically, means no cell division for new leaf
	else return __max(0.,R_max*pow((t-t_b)/deltaT1,alpha)*(t_e-t)/deltaT2);
}

