#include "stdafx.h"
#include "controller.h"
#include "initInfo.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <time.h>
#include <stdlib.h>
#ifndef FLOAT_EQ
#define EPSILON 0.001   // floating point comparison tolerance
#define FLOAT_EQ(x,v) (((v - EPSILON) < x) && (x <( v + EPSILON)))
#endif
#define comma ","
#ifdef _WIN32
const char* pathSymbol = "\\";
#else
const char* pathSymbol = "/";
#endif

inline double E_sat(double T) { return 0.6105 * exp(17.27 * T / (237.7 + T)); }

//#using <mscorlib.dll>
using namespace std;


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

CController::CController(const char* filename, const char* outfile, TInitInfo iniInfo)
{
	int size = 0;
	//string a = "architect.crp", l = "leaf.crp", n = "nitrogen.crp", p = "plantstress.crp"; // WS add new output file for interface
	string n = "nitrogen.crp", p = "plantstress.crp";
	time			 = NULL;
	weather 		 = NULL;
	plant			 = NULL;

	strcpy_s(varietyFile, filename);
	strcpy_s(cropFile, outfile);

	char* pathSymbol = (char*)calloc(256, sizeof(char));
	pathSymbol = (char*)"/\\"; //for both linux and windwos
	std::string plantstressFile = "plantstress.crp";
	std::string nitrogenFile = "nitrogen.crp";
	std::string basePath;
	std::size_t found;
	std::string cropFileAsString = cropFile;
	// find last path separator and break path from file name	
	found = cropFileAsString.find_last_of(pathSymbol);
	basePath = cropFileAsString.substr(0, found);
	// now get filename to use as a root
	std::string fileName = cropFileAsString.substr(found + 1);
	std::string root = fileName.substr(0, fileName.length() - 3);

	//size = strlen(outfile);
	//for (int i=size ;i>0;i--) //extract file location for graphics folder
	//{
	//	if (outfile[i] == '\\')
	//	{
	//		strncpy_s(directory, outfile, i+1);
			//strcpy_s(archFile, directory);
			//strcpy_s(leafFile, directory);
	//		strcpy_s(nitrogenFile, directory);
	//		strcpy_s(plantstressFile, directory);
			//string temp = archFile;
			//temp.insert(strlen(archFile),a);
			//strcpy_s(archFile,temp.c_str());
			//temp = leafFile;
			//temp.insert(strlen(leafFile),l);
			//strcpy_s(leafFile,temp.c_str());
	//		string temp = nitrogenFile;
	//		temp.insert(strlen(nitrogenFile),n);
	//		strcpy_s(nitrogenFile,temp.c_str());
	//		temp = plantstressFile;
	//		temp.insert(strlen(plantstressFile), p);
	//		strcpy_s(plantstressFile, temp.c_str());
	//		break;
	//	}
	// }

	iCur = 0;
	weatherFormat = ICASA;
	initInfo = iniInfo;
    firstDayOfSim = 0;
	lastDayOfSim = 365;
    initialize();
	errorFlag = 0;
}
CController::~CController()
{
	if ( time != NULL )
		delete time;
	if ( weather != NULL )
		delete [] weather ;
	if ( plant != NULL )
		delete plant;
}
//**********************************************************************
void CController::initialize()
{
	Timer dConvert; //object to convert dates
	int mm = 0, dd = 0, yy = 0; //for calendar dates
	char* Buffer=(char*)calloc(256,sizeof(char)); //to hold duumy strings from variety file
	cout <<setiosflags(ios::left) << endl
		<< " ***********************************************************" << endl
		<< " *          SpudSim: A Simulation Model for Potato         *" << endl
		<< " *                     VERSION  2.2                        *" << endl
		<< " *   USDA-ARS, Adaptive Cropping Systems Laboratory        *" << endl
		<< " ***********************************************************" << endl
		<< endl << endl;
// output headings of files

// Plant file for output
		ofstream cropOut(cropFile, ios::out); 
		cropOut << setiosflags(ios::right) << setiosflags(ios::fixed)
 			<< setw(4) << "date" << comma
 			<< setw(4) << "jday" << comma
			<< setw(4) << "time" << comma
			<< setw(5) << "LA/pl" << comma
			<< setw(3) << "LAI" << comma
			<< setw(3) << "PFD" << comma
			<< setw(6) << "SolRad" << comma
			<< setw(4) << "Tair" << comma
			<< setw(4) << "Tcan" << comma
			<< setw(6) << "Pgross" << comma
			<< setw(5) << "Rg+Rm" << comma
			<< setw(6) << "Tr-Pot" << comma
			<< setw(6) << "Tr-Act" << comma
			<< setw(5) << "Stage" << comma
			<< setw(7) << "totalDM" << comma
			<< setw(6) << "leafDM" << comma
			<< setw(6) << "stemDM" << comma
			<< setw(6) << "rootDM" << comma
			<< setw(7) << "tuberDM" << comma
			<< setw(6) << "deadDM" << comma
			<< setw(5) << "Cdead" << comma
			<< setw(5) << "Cpool" << comma
			<< setw(5) << "LWPpd" << comma
			<< setw(6) << "LWPave" << comma
			<< setw(6) << "gs_ave" << comma
			<< setw(8) << "Nstress1" << comma
			<< setw(8) << "Nstress2" << comma
			<< setw(8) << "Wstress1" << comma
			<< setw(8) << "Wstress2" << comma
			<< setw(8) << "Wstress3" 
		    << endl;
		//ofstream canOut(archFile, ios::out);
		//canOut << setiosflags(ios::right) << setiosflags(ios::fixed)
		//	<< setw(4) << "date" << comma
		//	<< setw(4) << "jday" << comma
		//	<< setw(4) << "time" << comma
		//	<< setw(9) << "mbrstemDM" << comma
		//	<< setw(9) << "mbrleafDM" << comma
		//	<< setw(8) << "mbrleafA" << comma
		//	<< setw(8) << "mbrnodes" << comma
		//	<< setw(9) << "bbrstemDM" << comma
		//	<< setw(9) << "bbrleafDM" << comma
		//	<< setw(8) << "bbrleafA" << comma
		//	<< setw(8) << "bbrstems" << comma
		//	<< setw(8) << "bbrnodes" << comma
		//	<< setw(9) << "abrstemDM" << comma
		//	<< setw(9) << "abrleafDM" << comma
		//	<< setw(8) << "abrleafA" << comma
		//	<< setw(8) << "abrstems" << comma
		//	<< setw(8) << "abrnodes" 
		//	<< endl;
		//ofstream leafout(leafFile, ios::out);
		//leafout << setiosflags(ios::right) << setiosflags(ios::fixed)
		//	<< setw(4) << "date" << comma
		//	<< setw(4) << "jday" << comma
		//	<< setw(4) << "time" << comma
		//	<< setw(5) << "L1:dm" << comma
		//	<< setw(7) << "L1:area" << comma
		//	<< setw(7) << "L1:page" << comma
		//	<< setw(7) << "L1:cage" << comma
		//	<< setw(5) << "L2:dm" << comma
		//	<< setw(7) << "L2:area" << comma
		//	<< setw(7) << "L2:page" << comma
		//	<< setw(7) << "L2:cage" << comma
		//	<< setw(5) << "L3:dm" << comma
		//	<< setw(7) << "L3:area" << comma
		//	<< setw(7) << "L3:page" << comma
		//	<< setw(7) << "L3:cage" << comma
		//	<< setw(5) << "L4:dm" << comma
		//	<< setw(7) << "L4:area" << comma
		//	<< setw(7) << "L4:page" << comma
		//	<< setw(7) << "L4:cage" << comma
		//	<< setw(5) << "L5:dm" << comma
		//	<< setw(7) << "L5:area" << comma
		//	<< setw(7) << "L5:page" << comma
		//	<< setw(7) << "L5:cage" << comma
		//	<< setw(5) << "L6:dm" << comma
		//	<< setw(7) << "L6:area" << comma
		//	<< setw(7) << "L6:page" << comma
		//	<< setw(7) << "L6:cage" << comma
		//	<< setw(5) << "L7:dm" << comma
		//	<< setw(7) << "L7:area" << comma
		//	<< setw(7) << "L7:page" << comma
		//	<< setw(7) << "L7:cage" << comma
		//	<< setw(5) << "L8:dm" << comma
		//	<< setw(7) << "L8:area" << comma
		//	<< setw(7) << "L8:page" << comma
		//	<< setw(7) << "L8:cage" << comma
		//	<< setw(5) << "L9:dm" << comma
		//	<< setw(7) << "L9:area" << comma
		//	<< setw(7) << "L9:page" << comma
		//	<< setw(7) << "L9:cage" << comma
		//	<< setw(6) << "L10:dm" << comma
		//	<< setw(8) << "L10:area" << comma
		//	<< setw(8) << "L10:page" << comma
		//	<< setw(8) << "L10:cage"
		//	<< endl;
		ofstream nitrogenout(nitrogenFile, ios::out);
		nitrogenout << setiosflags(ios::right) << setiosflags(ios::fixed)
			<< setw(4) << "date" << comma
			<< setw(4) << "jday" << comma
			<< setw(4) << "time" << comma
			<< setw(5) << "tot_N" << comma
			<< setw(6) << "leaf_N" << comma
			<< setw(6) << "stem_N" << comma
			<< setw(6) << "root_N" << comma
			<< setw(7) << "tuber_N" << comma
			<< setw(6) << "dead_N" << comma
			<< setw(7) << "tot_N_C" << comma
			<< setw(8) << "leaf_N_C" << comma
			<< setw(8) << "stem_N_C" << comma
			<< setw(8) << "root_N_C" << comma
			<< setw(8) << "tubr_N_C" << comma
			<< setw(8) << "N_uptake" << comma
			<< setw(8) << "N_demand" << comma
			<< setw(6) << "Seed N" << comma
			<< setw(7) << "Nstress"
			<< endl;
		ofstream plantstressout(plantstressFile, ios::out);
		plantstressout << setiosflags(ios::right) << setiosflags(ios::fixed)
			<< setw(4) << "date" << comma
			<< setw(4) << "jday" << comma
			<< setw(4) << "time" << comma
			<< setw(18) << "waterstressfactor" << comma
			<< setw(15) << "PSIEffect_leaf" << comma
			<< setw(13) << "NEffect_leaf" << comma
			<< setw(13) << "PSIEffect_Pn" << comma
			<< setw(11) << "NEffect_Pn" << comma
			<< setw(10) << "Dev_stage" << comma
			<< setw(9) << "Heat_veg" << comma
			<< setw(11) << "Heat_repre"
			<< endl;
	try
	{
		ifstream cfs(varietyFile, ios::in);
		if (!cfs)
		{
			throw "Variety File not found.";
		}
		cfs.getline(initInfo.description, sizeof(initInfo.description),'\n');
		cfs.getline(initInfo.cultivar, sizeof(initInfo.cultivar),'\n');
		cfs.getline(Buffer, 255,'\n');
		cfs.getline(Buffer, 255,'\n');
		cfs >> initInfo.A1 >> initInfo.A2 >> initInfo.A3 >> initInfo.A4 >> initInfo.A5 >> initInfo.A6 >> initInfo.A7 >> initInfo.A8 >> initInfo.A9 >> initInfo.A10;
		cfs >> initInfo.G1 >> initInfo.G2 >> initInfo.G3 >> initInfo.G4;
		
		dConvert.caldat(initInfo.sowingDay, mm, dd, yy);
		if (cfs.eof()) cfs.close();
		cout << "Read variety file : " << varietyFile << endl <<endl;
		cout << "Simulation information for:" << endl;
		cout << setiosflags(ios::left)
			<< setw(20) << "Description: " << initInfo.description << endl << endl
			<< setw(10)	<< "Cultivar: " << initInfo.cultivar << endl
			<< setw(6) << "Year: " << initInfo.year << endl
			<< setw(6) << "Sowing day: " << initInfo.sowingDay << endl 
			<< setw(6) << "Emergence day:" << initInfo.emergeDay << endl
			<< setw(6) << "TimeStep (min): " << initInfo.timeStep << endl
			<< setw(6) << "Outfile Timestep (min): " <<initInfo.outtimeStep <<endl
			<< setw(6) << "Leaf model (0-individual), (1-big): "<<initInfo.bigleaf <<endl
			<< setw(6) << "Nitrogen stress turned off? " <<initInfo.Nitrogen_stress_off <<endl
			<< setw(6) << "Nater stress turned off? " <<initInfo.Water_stress_off << endl
			<< setw(6) << "Water stress simulation method: " <<initInfo.Water_stress_simulation_type << endl
			;
	}
	catch(const char* message)
	{
		cerr << message << "\n";
		exit(1);
	}

	firstDayOfSim = int(initInfo.beginDay);
	lastDayOfSim = int(initInfo.endDay);
	SowingDay = int(initInfo.sowingDay);
	EmergenceDay = int(initInfo.emergeDay);
	cropEmerged = false;
	cropHarvested = false;
	dConvert.caldat(firstDayOfSim,mm,dd,yy);
    time = new Timer(dd, mm, yy, initInfo.timeStep/60.0); // Timer class gets stepsize in hours
	int dim = (int)(((lastDayOfSim+1)-firstDayOfSim)*(24*60/initInfo.timeStep)); // counting total records of weather data 
	weather = new TWeather[dim];
	plant	= new CPlant(initInfo);
}

void CController::readWeatherFrom2DSOIL(const TWeather & wthr)
{
	/***********/
	/* Assigns weather structure values (assigned in crop.cpp) 
	/* to corresponding weather array element*/

	weatherFormat = DDSOIL;
	weather[iCur] = wthr;
	ET_supply = wthr.ET_supply;
	weather[iCur].daytime = wthr.jday + wthr.time;
}


void CController::createOutputFiles()
{
    char dname[120];
	// find output file extension (e.g. ***.out)  
	int i=0;
	for (int i=0;i <= sizeof(outputFile);i++)
    {
    	if (outputFile[i] == '\0')
    		break;
    }
	int q=0;
    for (int q=i;q >= 0;q--)
    {
    	if (outputFile[q] == '.')
    	{
    		break;
    	}
    }
	
	// find output file name (e.g. name.***)
    for (int y=0;y<=q;y++)
    {
    	dname[y] = outputFile[y];
    }
    dname[q+1] = '\0';


	// use same name for crop and log files and add appropriate extensions
    strcpy_s(cropFile, dname);
	strcpy_s(logFile, dname);
	strcat_s(cropFile, "crp");
    strcat_s(logFile, "log");
	ofstream ostr(logFile, ios::out);
	ostr << "*** Notes and Warnings ***\n";
	ostr.close();
}

int CController::run(const TWeather & wthr, double lwp_dawn) // this version used with 2dSpudSim
{
/* Current run method used with 2DSPUDSIM
/* Transfer current weather data from 2DSOIL and places in local Weather array
/* Access plant growth and development methods through plant.update method
/* Outputs crop module specific information to CropFile
/* Adjust simulation time-counter*/
	readWeatherFrom2DSOIL(wthr);
	this->lwp_dawn = lwp_dawn;
    if (weather[iCur].jday >= initInfo.sowingDay && weather[iCur].jday <= lastDayOfSim)
	{
		plant->update(iCur,weather[iCur], lwp_dawn, initInfo, time->get_hour());
		RootWeightFrom2DSOIL=wthr.TotalRootWeight;
		MaxRootDepth=        wthr.MaxRootDepth;
		AvailableWater=      wthr.ThetaAvailRZ;
		//plant->leaf_output();
		//outputToCropFile();
//		if (FLOAT_EQ(weather[iCur].time,0.5))
		if (((iCur==0) || (fmod(iCur+1,(initInfo.outtimeStep/initInfo.timeStep)) < 0.001))) /*output at end of desired timestep*/
		{
			//if (initInfo.bigleaf == 0) plant->leaf_output(); //output data for individual leaf information
//			cout << weather[i].jday <<endl;
			outputToCropFile();
//			cout << weather[iCur].jday <<endl;
		}
		iCur++;
		time->step();
	} 
    return 0;
}


void CController::outputToCropFile()
{

	int mm = 0, id = 0, iyyy = 0;
	string DateForOutput;
	double av_gs = plant->get_conductance();
	if (!plant->get_develop()->Emerged())
	{
		av_gs = 0;
		return;
	}
	double vpd = plant->get_VPD();
	if (vpd<0)
	{
		vpd = 0;
	}
	time->caldat(weather[iCur].jday, mm,id,iyyy);
#if 0
	DateForOutput.Format("%.2d/%.2d/%4i",mm,id,iyyy);
#else
	char DateForOutputBuff[16];
	sprintf(DateForOutputBuff,"%.2d/%.2d/%4i",mm,id,iyyy);
	DateForOutput = DateForOutputBuff;
#endif
	if (initInfo.outtimeStep>1400) //24 hour output
	{
		ofstream ostr(cropFile, ios::app);
		ostr << setiosflags(ios::right) 
			<< setiosflags(ios::fixed)
			<< setw(9) << DateForOutput << comma
 			<< setw(6) << weather[iCur].jday << comma
			<< setw(8) << setprecision(3) << int(weather[iCur].time*24+0.1) << comma
			<< setw(8) << setprecision(2) << plant->get_develop()->get_LvsAppeared() << comma
			<< setw(8) << setprecision(2) << plant->get_greenLAI() << comma
			<< setw(8) << setprecision(2) << plant->get_Parave() << comma//for 24h output, mol PAR m-2 d-1
			//<< setw(8) << setprecision(2) << weather[iCur].PFD << comma//for hourlyoutput, umol m-2 s-1
			<< setw(8) << setprecision(2) << plant->get_Sradave() << comma//for 24h output, MW m-2 d-1
			//<< setw(8) << setprecision(2) << weather[iCur].solRad << comma//for hourly output, W m-2
			<< setw(8) << setprecision(2) << plant->get_Tdayave() << comma//for 24h output, average 24 hour air temperature, C
			//<< setw(8) << setprecision(2) << weather[iCur].airT << comma//for hourly output, current air T, C
			<< setw(8) << setprecision(2) << plant->get_tmpr2() << comma//for 24h output, average 24 hour canopy temperature, C
			//<< setw(8) << setprecision(2) << plant->get_tmpr() << comma//for hourly output, canopy temperature at last time-step, C
			//<< setw(8) << setprecision(2) << weather[iCur].psil_ << comma//for hourly output, leaf water potential at current time-step, MPa
			//<< setw(8) << setprecision(2) << lwp_dawn << comma//predawn leaf water potential, MPa
			<< setw(8) << setprecision(2) << plant->get_Pg2() << comma //for 24h output, gross photosynthetic rate in g CHO plant-1 day-1 (this is also whats used in CLASSIM)
			//<< setw(8) << setprecision(2) << plant->get_Pg() << comma//for hourly output, gross photosynthetic rate in umol co2 m-2 ground s-1
			<< setw(8) << setprecision(3) << plant->get_Rg()+plant->get_Rm() << comma// for 24 hour output, total respiration in g CHO plant-1 day-1
			//<<setw(8) << setprecision(3) << plant->get_Resp() << comma//for hourly output, umol co2 m-2 ground day-s
			<< setw(8) << setprecision(2) << plant->get_cumulativeET() << comma//for 24 hour output, plant transpiration (not ET!), in g H2O plant-1 day-1
			//<< setw(8) << setprecision(2) << plant->get_ET() << comma//for hourly output, transpiration rate, in mmol h2o m-2 ground s-1
			<< setw(8) << setprecision(2) << plant->get_totalSoilWaterUptake() << comma// for 24 hour output, plant transpiration per 2DSOIL, g H2o plant-1 day-1
			//<< setw(8) << setprecision(2) << plant->get_HourlySoilWaterUptake() << comma // for hourly output, plant transpiration per 2dsoil, g H2O plant-1 hour-1
			<< setw(8) << setprecision(2) << plant->get_develop()->get_xstage() << comma
			<< setw(8) << setprecision(3) << plant->get_totalMass() << comma
			<< setw(8) << setprecision(3) << plant->get_leafMass() << comma
			<< setw(8) << setprecision(3) << plant->get_stemMass() << comma
			<< setw(8) << setprecision(3) << plant->get_rootMass() << comma
			<< setw(8) << setprecision(3) << plant->get_tuberMass() << comma
			<< setw(8) << setprecision(3) << plant->get_deadMass() << comma//total amount of seneseced leaf
			<< setw(8) << setprecision(3) << plant->get_C_deadpool_season() << comma//for 24 hour, total C mobilized from dead leaf
			<< setw(8) << setprecision(3) << plant->get_C_pool() << comma
			
			//<< setw(8) << setprecision(3) << plant->get_C_pool_root()	 << comma		
			<< setw(8) << setprecision(3) << lwp_dawn/ 10 << comma//pre dawn leaf water potential in MPa
			<< setw(8) << setprecision(3) << plant->get_averageLWP() << comma//get daily averaged lwp in MPa
			//<< setw(8) << setprecision(3) << weather->psil_ << comma//get hourly psileaf
			<< setw(8) << setprecision(3) << plant->get_averagegs() << comma// get daily averaged stomatal conductance
			<< setw(8) << setprecision(3) << plant->get_Nstressfactorone() << comma// photosynthetic reduction
			<< setw(8) << setprecision(3) << plant->get_Nstressfactortwo() << comma // root increase factor
			<< setw(8) << setprecision(3) << plant->get_waterstressfactor() << comma// canopy senescence factor
			//<< setw(8) << setprecision(3) << plant->get_develop()->get_partub() << comma
			<< setw(8) << setprecision(3) << plant->get_psistress_gs_factor() << comma//bulk leaf water potential stress on Pn
			<< setw(8) << setprecision(3) << plant->get_nodalUnit()->get_leaf()->get_psi_leafexpansion_stress() // bulk leaf water potential stress on leaf expansion
			<< endl;
	}
	else //hourly output
	{
		ofstream ostr(cropFile, ios::app);
		ostr << setiosflags(ios::right) 
			<< setiosflags(ios::fixed)
			<< setw(9) << DateForOutput << comma
 			<< setw(6) << weather[iCur].jday << comma
			<< setw(8) << setprecision(3) << int(weather[iCur].time*24+0.1) << comma
			<< setw(8) << setprecision(2) << plant->get_develop()->get_LvsAppeared() << comma
			<< setw(8) << setprecision(2) << plant->get_greenLAI() << comma
			<< setw(8) << setprecision(2) << weather[iCur].PFD << comma//for hourlyoutput, umol m-2 s-1
			<< setw(8) << setprecision(2) << weather[iCur].solRad << comma//for hourly output, W m-2
			<< setw(8) << setprecision(2) << weather[iCur].airT << comma//for hourly output, current air T, C
			<< setw(8) << setprecision(2) << plant->get_tmpr() << comma//for hourly output, canopy temperature at last time-step, C
			//<< setw(8) << setprecision(2) << plant->get_Pg() << comma//for hourly output, gross photosynthetic rate in umol co2 m-2 ground s-1
			<< setw(8) << setprecision(2) << plant->get_Pg()*3600./initInfo.plantDensity*30./1000000. << comma//for hourly output, gross photosynthetic rate in g CHO plant / h for CLASSIM
			//<< setw(8) << setprecision(3) << plant->get_Resp() << comma//for hourly output, umol co2 m-2 ground day-s
			<< setw(8) << setprecision(3) << plant->get_Resp()*3600./initInfo.plantDensity*30./1000000. << comma//for hourly output, in g CHO plant / h for CLASSIM
			<< setw(8) << setprecision(2) << plant->get_ET()*0.018*3600./initInfo.plantDensity << comma//for hourly output, transpiration rate, in g h2o plant-1 h-1
			<< setw(8) << setprecision(2) << plant->get_HourlySoilWaterUptake() << comma // for hourly output, plant transpiration per 2dsoil, g H2O plant-1 hour-1
			<< setw(8) << setprecision(2) << plant->get_develop()->get_xstage() << comma
			<< setw(8) << setprecision(3) << plant->get_totalMass() << comma
			<< setw(8) << setprecision(3) << plant->get_leafMass() << comma
			<< setw(8) << setprecision(3) << plant->get_stemMass() << comma
			<< setw(8) << setprecision(3) << plant->get_rootMass() << comma
			<< setw(8) << setprecision(3) << plant->get_tuberMass() << comma
			<< setw(8) << setprecision(3) << plant->get_deadMass() << comma
			<< setw(8) << setprecision(3) << plant->get_C_deadpool_season() << comma
			<< setw(8) << setprecision(3) << plant->get_C_pool() << comma
			
			//<< setw(8) << setprecision(3) << plant->get_sunPFD()<< comma
			//<< setw(8) << setprecision(3) << plant->get_shadePFD()<< comma
			//<< setw(8) << setprecision(3) << plant->get_sunLAI()<< comma
			//<< setw(8) << setprecision(3) << plant->get_shadeLAI()<< comma
			//<< setw(8) << setprecision(3) << plant->get_Pnleafsunlit()<< comma
			//<< setw(8) << setprecision(3) << plant->get_Pnleafshaded()<< comma
			//<< setw(8) << setprecision(3) << plant->get_Transpirationsunlit()<< comma
			//<< setw(8) << setprecision(3) << plant->get_Transpirationshaded()<< comma
			
			//<< setw(8) << setprecision(3) << plant->get_C_pool_root()	<< comma					
			<< setw(8) << setprecision(3) << lwp_dawn/ 10 << comma//pre dawn leaf water potential in MPa
			//<< setw(8) << setprecision(3) << plant->get_rootGrowth()<< comma
			//<< setw(8) << setprecision(3) << plant->get_shootPart()<< comma
			//<< setw(8) << setprecision(3) << weather[iCur].pcrs<< comma
			//<< setw(8) << setprecision(3) << weather[iCur].TotalRootWeight<< comma
			//<< setw(8) << setprecision(3) << plant->get_C_seed_used()<< comma
			//<< setw(8) << setprecision(3) << plant->get_transpirationSupplyDemandRatio()<< comma
			//<< setw(8) << setprecision(3) << plant->get_develop()->get_xstage()<< comma
			//<< setw(8) << setprecision(3) << plant->get_C_pool_used()<< comma
			//<< setw(8) << setprecision(3) << plant->get_averageLWP() << comma//get daily averaged lwp in MPa
			<< setw(8) << setprecision(3) << weather[iCur].psil_ << comma//get hourly psileaf
			<< setw(8) << setprecision(3) << plant->get_conductance() << comma//hourly stomatal conductance
			//<< setw(8) << setprecision(3) << plant->get_develop()->get_xstage() << comma
			<< setw(8) << setprecision(3) << plant->get_Nstressfactorone() << comma// photosynthetic reduction
			<< setw(8) << setprecision(3) << plant->get_Nstressfactortwo() << comma // root increase factor
			<< setw(8) << setprecision(3) << plant->get_waterstressfactor() << comma// canopy senescence factor
			//<< setw(8) << setprecision(3) << plant->get_develop()->get_partub()<< comma
			<< setw(8) << setprecision(3) << plant->get_psistress_gs_factor() << comma //bulk leaf water potential stress on Pn
			<< setw(8) << setprecision(3) << plant->get_nodalUnit()->get_leaf()->get_psi_leafexpansion_stress() // bulk leaf water potential stress on leaf expansion



			<< endl;
	}
		//ofstream ostr2(archFile, ios::app);
		//ostr2 << setiosflags(ios::right)
		//	<< setiosflags(ios::fixed)
		//	<< setw(9) << DateForOutput << comma
		//	<< setw(6) << weather[iCur].jday << comma
		//	<< setw(8) << setprecision(3) << weather[iCur].time << comma
		//	<< setw(10) << setprecision(2) << plant->get_mainstemMass() << comma
		//	<< setw(10) << setprecision(2) << plant->get_mainleafMass() << comma
		//	<< setw(9) << setprecision(2) << plant->get_mainleafArea() << comma
		//	<< setw(9) << setprecision(0) << plant->get_mainleafNo() << comma
		//	<< setw(10) << setprecision(2) << plant->get_basalstemMass() << comma
		//	<< setw(10) << setprecision(2) << plant->get_basalleafMass() << comma
		//	<< setw(9) << setprecision(2) << plant->get_basalleafArea() << comma
		//	<< setw(9) << setprecision(0) << plant->get_basalstemNo() << comma
		//	<< setw(9) << setprecision(0) << plant->get_basalleafNo() << comma
		//	<< setw(10) << setprecision(2) << plant->get_apicalstemMass() << comma
		//	<< setw(10) << setprecision(2) << plant->get_apicalleafMass() << comma
		//	<< setw(9) << setprecision(2) << plant->get_apicalleafArea() << comma
		//	<< setw(9) << setprecision(0) << plant->get_apicalstemNo() << comma
		//	<< setw(9) << setprecision(0) << plant->get_apicalleafNo()
		//	<< endl;
	
		// output  leaf growth from specified nodes
		//if (plant->lmass2 < 0){ //leaf 8 not initiated yet
		//ofstream ostr3(leafFile, ios::app);
		//ostr3 << setiosflags(ios::right) << setiosflags(ios::fixed)
		//	<< setw(9) << DateForOutput << comma
		//	<< setw(6) << weather[iCur].jday << comma
		//	<< setw(8) << setprecision(3) << weather[iCur].time << comma
		//	<< setw(8) << setprecision(2) << plant->area[0] << comma
		//	<< setw(8) << setprecision(2) << plant->area[1] << comma
		//	<< setw(8) << setprecision(2) << plant->area[2] << comma
		//	<< setw(8) << setprecision(2) << plant->area[3] << comma
		//	<< setw(8) << setprecision(2) << plant->area[4] << comma
		//	<< setw(8) << setprecision(2) << plant->area[5] << comma
		//	<< setw(8) << setprecision(2) << plant->area[6] << comma
		//	<< setw(8) << setprecision(2) << plant->area[7] << comma
		//	<< setw(8) << setprecision(2) << plant->area[8] << comma
		//	<< setw(8) << setprecision(2) << plant->area[9] << comma
		//	<< setw(8) << setprecision(2) << plant->area[10] << comma
		//	<< setw(8) << setprecision(2) << plant->area[11] << comma
		//	<< setw(8) << setprecision(2) << plant->area[12] << comma
		//	<< setw(8) << setprecision(2) << plant->area[13] << comma
		//	<< setw(8) << setprecision(2) << plant->area[14] << comma
		//	<< setw(8) << setprecision(2) << plant->area[15] << comma
		//	<< setw(8) << setprecision(2) << plant->area[16] << comma
		//	<< setw(8) << setprecision(2) << plant->area[17] << comma
		//	<< setw(8) << setprecision(2) << plant->area[18] << comma
		//	<< setw(8) << setprecision(2) << plant->area[19] << comma
		//	<< setw(8) << setprecision(2) << plant->area[20] << comma
		//	<< setw(8) << setprecision(2) << plant->area[21] << comma
		//	<< setw(8) << setprecision(2) << plant->area[22] << comma
		//	<< setw(8) << setprecision(2) << plant->area[23] << comma
		//	<< setw(8) << setprecision(2) << plant->area[24] << comma
		//	<< setw(8) << setprecision(2) << plant->area[25] << comma
		//	<< setw(8) << setprecision(2) << plant->area[26] << comma
		//	<< setw(8) << setprecision(2) << plant->area[27] << comma
		//	<< setw(8) << setprecision(2) << plant->area[28] << comma
		//	<< setw(8) << setprecision(2) << plant->area[29] << comma
		//	<< setw(8) << setprecision(2) << plant->area[30] << comma
		//	<< setw(8) << setprecision(2) << plant->area[31]

			/*
			<< setw(8) << setprecision(2) << plant->lmass[0]<< comma
			<< setw(8) << setprecision(2) << plant->area[0]<< comma
			<< setw(8) << setprecision(2) << plant->page[0]<< comma
			<< setw(8) << setprecision(2) << plant->cage[0]<< comma
			<< setw(8) << setprecision(2) << plant->lmass[1]<< comma
			<< setw(8) << setprecision(2) << plant->area[1]<< comma
			<< setw(8) << setprecision(2) << plant->page[1]<< comma
			<< setw(8) << setprecision(2) << plant->cage[1]<< comma
			<< setw(8) << setprecision(2) << plant->lmass[2]<< comma
			<< setw(8) << setprecision(2) << plant->area[2]<< comma
			<< setw(8) << setprecision(2) << plant->page[2]<< comma
			<< setw(8) << setprecision(2) << plant->cage[2]<< comma
			<< setw(8) << setprecision(2) << plant->lmass[3]<< comma
			<< setw(8) << setprecision(2) << plant->area[3]<< comma
			<< setw(8) << setprecision(2) << plant->page[3]<< comma
			<< setw(8) << setprecision(2) << plant->cage[3]<< comma
			<< setw(8) << setprecision(2) << plant->lmass[4]<< comma
			<< setw(8) << setprecision(2) << plant->area[4]<< comma
			<< setw(8) << setprecision(2) << plant->page[4]<< comma
			<< setw(8) << setprecision(2) << plant->cage[4]<< comma
			<< setw(8) << setprecision(2) << plant->lmass[5]<< comma
			<< setw(8) << setprecision(2) << plant->area[5]<< comma
			<< setw(8) << setprecision(2) << plant->page[5]<< comma
			<< setw(8) << setprecision(2) << plant->cage[5]<< comma
			<< setw(8) << setprecision(2) << plant->lmass[6]<< comma
			<< setw(8) << setprecision(2) << plant->area[6]<< comma
			<< setw(8) << setprecision(2) << plant->page[6]<< comma
			<< setw(8) << setprecision(2) << plant->cage[6]<< comma
			<< setw(8) << setprecision(2) << plant->lmass[7]<< comma
			<< setw(8) << setprecision(2) << plant->area[7]<< comma
			<< setw(8) << setprecision(2) << plant->page[7]<< comma
			<< setw(8) << setprecision(2) << plant->cage[7]<< comma
			<< setw(8) << setprecision(2) << plant->lmass[8]<< comma
			<< setw(8) << setprecision(2) << plant->area[8]<< comma
			<< setw(8) << setprecision(2) << plant->page[8]<< comma
			<< setw(8) << setprecision(2) << plant->cage[8]<< comma
			<< setw(8) << setprecision(2) << plant->lmass[9]<< comma
			<< setw(8) << setprecision(2) << plant->area[9]<< comma
			<< setw(8) << setprecision(2) << plant->page[9]<< comma
			<< setw(8) << setprecision(2) << plant->cage[9]<< comma
			<< setw(8) << setprecision(2) << plant->lmass[10]<< comma
			<< setw(8) << setprecision(2) << plant->area[10]<< comma
			<< setw(8) << setprecision(2) << plant->page[10]<< comma
			<< setw(8) << setprecision(2) << plant->cage[10]<< comma
			*/
			//<< endl;
	if (initInfo.outtimeStep>1400) //24 hour output
	{
		ofstream ostr4(nitrogenFile, ios::app);
		ostr4 << setiosflags(ios::right) << setiosflags(ios::fixed)
			<< setw(10) << DateForOutput << comma
			<< setw(6) << weather[iCur].jday << comma
			<< setw(8) << setprecision(3) << int(weather[iCur].time * 24 + 0.1) << comma
			<< setw(8) << setprecision(2) << this->getPlant()->get_totalPlantNitrogen() << comma //total N content in plant, g N plant-1
			<< setw(8) << setprecision(2) << this->getPlant()->get_totalLeafNitrogen() << comma// total leaf N,  g N plant-1
			<< setw(8) << setprecision(2) << this->getPlant()->get_totalStemNitrogen() << comma// total stem N,  g N plant-1
			<< setw(8) << setprecision(2) << this->getPlant()->get_totalRootNitrogen() << comma// total root N,  g N plant-1
			<< setw(8) << setprecision(2) << this->getPlant()->get_totalTuberNitrogen() << comma// total tuber N, g N plant-1
			<< setw(8) << setprecision(2) << this->getPlant()->get_totalDeadNitrogen() << comma // total dead N lost to senescecnce, g N plant-1
			<< setw(9) << setprecision(3) << this->getPlant()->get_NCRatioTotalPlant() << comma // N:C ratio of whole plant g N g-1
			<< setw(9) << setprecision(3) << this->getPlant()->get_NCRatioLeaves() << comma // N:C ratio of leaves, g N g-1
			<< setw(9) << setprecision(3) << this->getPlant()->get_NCRatioStems() << comma// N:C ratio of stems, g N g-1
			<< setw(9) << setprecision(3) << this->getPlant()->get_NCRatioRoots() << comma // N:C ratio of roots, g N g-1
			<< setw(9) << setprecision(3) << this->getPlant()->get_NCRatioTubers() << comma // N:C ratio of tubers, g N g-1
			<< setw(9) << setprecision(3) << this->getPlant()->get_CumulativeNitrogenSoilUptake() << comma//total N uptake, cumulative, g N plant-1
			<< setw(9) << setprecision(3) << this->getPlant()->get_CumulativeNitrogenDemand() << comma//total N demand, cumulative, g  plant-1
			<< setw(9) << setprecision(3) << this->getPlant()->get_seedNitrogenUsed() << comma// cumulative N translocated from seedpiece, g g N plant-1
			<< setw(8) << setprecision(1) << this->getPlant()->get_Nstressfactor() // 0 to 1 stress factor, less than 1 indicates some level of N stress
			<< endl;
	}
	else //hourly output
	{
		ofstream ostr4(nitrogenFile, ios::app);
		ostr4 << setiosflags(ios::right) << setiosflags(ios::fixed)
			<< setw(10) << DateForOutput << comma
			<< setw(6) << weather[iCur].jday << comma
			<< setw(8) << setprecision(3) << int(weather[iCur].time * 24 + 0.1) << comma
			<< setw(8) << setprecision(2) << this->getPlant()->get_totalPlantNitrogen() << comma//total N content in plant, g N plant-1
			<< setw(8) << setprecision(2) << this->getPlant()->get_totalLeafNitrogen() << comma// total leaf N,  g N plant-1
			<< setw(8) << setprecision(2) << this->getPlant()->get_totalStemNitrogen() << comma// total stem N,  g N plant-1
			<< setw(8) << setprecision(2) << this->getPlant()->get_totalRootNitrogen() << comma // total root N,  g N plant-1
			<< setw(8) << setprecision(2) << this->getPlant()->get_totalTuberNitrogen() << comma // total tuber N,g N plant-1
			<< setw(8) << setprecision(2) << this->getPlant()->get_totalDeadNitrogen() << comma// total dead N lost to senescecnce, g N plant-1
			<< setw(9) << setprecision(3) << this->getPlant()->get_NCRatioTotalPlant() << comma// N:C ratio of whole plant g N g-1
			<< setw(9) << setprecision(3) << this->getPlant()->get_NCRatioLeaves() << comma // N:C ratio of leaves, g N g-1
			<< setw(9) << setprecision(3) << this->getPlant()->get_NCRatioStems() << comma// N:C ratio of stems, g N g-1
			<< setw(9) << setprecision(3) << this->getPlant()->get_NCRatioRoots() << comma// N:C ratio of roots, g N g-1
			<< setw(9) << setprecision(3) << this->getPlant()->get_NCRatioTubers() << comma// N:C ratio of tubers, g N g-1
			<< setw(9) << setprecision(3) << this->getPlant()->get_CumulativeNitrogenSoilUptake() << comma //total N uptake, cumulative, g N plant-1
			<< setw(9) << setprecision(3) << this->getPlant()->get_CumulativeNitrogenDemand() << comma//total N demand, cumulative, g  plant-1
			<< setw(9) << setprecision(3) << this->getPlant()->get_seedNitrogenUsed() << comma// cumulative N translocated from seedpiece, g g N plant-1
			<< setw(8) << setprecision(1) << this->getPlant()->get_Nstressfactor() // 0 to 1 stress factor, less than 1 indicates some level of N stress
			<< endl;
		ofstream ostr5(plantstressFile, ios::app);
		ostr5 << setiosflags(ios::right) << setiosflags(ios::fixed)
			<< setw(10) << DateForOutput << comma
			<< setw(6) << weather[iCur].jday << comma
			<< setw(8) << setprecision(3) << int(weather[iCur].time * 24 + 0.1) << comma
			<< setw(20) << setprecision(3) << plant->get_waterstressfactor() << comma
			<< setw(20) << setprecision(3) << plant->get_psi_leafexpansion_stress() << comma
			<< setw(20) << setprecision(3) << plant->get_Nstressfactortwo() << comma// From SIMPOTATO (NDEF2), 0 to 1; reduces whole plant leaf expansion, increases senescence
			<< setw(20) << setprecision(3) << plant->get_psistress_gs_factor() << comma //
			<< setw(20) << setprecision(3) << plant->get_Nstressfactorone() << comma // photosynthetic reduction
			<< setw(20) << setprecision(3) << plant->get_develop()->get_xstage() << comma
			<< setw(20) << setprecision(3) << plant->get_heat_veg_stressfactor() << comma // heat stress at vegetative growth set 1 for now
			<< setw(20) << setprecision(3) << plant->get_heat_repro_stressfactor()
			//need to add heat

			<< endl;
	}
	if (initInfo.outtimeStep > 1400) //24 hour output
	{
		ofstream ostr5(plantstressFile, ios::app);
		ostr5 << setiosflags(ios::right) << setiosflags(ios::fixed)
			<< setw(2) << DateForOutput << comma
			<< setw(8) << weather[iCur].jday << comma
			<< setw(8) << setprecision(3) << int(weather[iCur].time * 24 + 0.1) << comma
			<< setw(20) << setprecision(3) << plant->get_waterstressfactor() << comma
			<< setw(20) << setprecision(3) << plant->get_psi_leafexpansion_stress_24h() << comma
			<< setw(20) << setprecision(3) << plant->get_Nstressfactortwo_24h() << comma// From SIMPOTATO (NDEF2), 0 to 1; reduces whole plant leaf expansion, increases senescence
			<< setw(20) << setprecision(3) << plant->get_psistress_gs_factor_24h() << comma //
			<< setw(20) << setprecision(3) << plant->get_Nstressfactorone_24h() << comma // photosynthetic reduction
			<< setw(20) << setprecision(3) << plant->get_develop()->get_xstage() << comma // developmental stage
			<< setw(20) << setprecision(3) << plant->get_heat_veg_stressfactor() << comma // heat stress at vegetative growth set 1 for now
			<< setw(20) << setprecision(3) << plant->get_heat_repro_stressfactor()

			//total N content in plant, mg N plant-1
			<< endl;
	}
	return;
}
