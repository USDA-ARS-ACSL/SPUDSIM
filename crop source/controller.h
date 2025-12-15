//Class Controller
//
// Controller.h
//
//Based on CPM simulation_controller and Maizesim 0.0.01
#pragma once
#ifndef _CONTROLLER_H_
#define _CONTROLLER_H_
#include "timer.h"
#include "development.h"
#include "plant.h"
#include "weather.h"
#include "initInfo.h"
#ifndef FLOAT_EQ
#define EPSILON 0.001 //floating point comparison tolerance
#define FLOAT_EQ(x,v) (((v - EPSILON) < x) && (x <( v + EPSILON)))
#endif


class CController
{
private:
	enum InputDataFormat {DDSOIL, ICASA, SPAR, FACE}; // input data format
	bool cropEmerged, cropHarvested;
	int	firstDayOfSim = 0,	lastDayOfSim = 0, SowingDay = 0, EmergenceDay = 0;
	char varietyFile[218], outputFile[11], cropFile[218], logFile[216];
	double lwp_dawn;
	char archFile[133], leafFile[133], directory[133];
	//char weatherFile[120], varietyFile[120], outputFile[120], cropFile[120], logFile[120];
	std::string plantstressFile;
	std::string nitrogenFile;
	int iCur = 0 , // current record number
		errorFlag = 0; // 1 if problem, 0 if not

	TInitInfo			initInfo;  
	CPlant*             plant;
	Timer*				time;
	TWeather*			weather;
	CDevelopment*       develop;
	InputDataFormat     weatherFormat;

public:
	CController(const char*, const char*, TInitInfo);
    ~CController();
	void setErrStatus(int ier) {errorFlag = ier;}
	int getErrStatus() {return errorFlag;}
	//char* getWeatherFile() {return weatherFile;}
	char* getInitFile() {return varietyFile;}
	//char* getOutputFile() {return outputFile;}
	char* getLogFile() {return logFile;}
	CPlant * getPlant() {return plant;}

	void addOutputMessage(char*);
	void readWeatherFrom2DSOIL(const TWeather &);
	void createOutputFiles();
	void outputToCropFile();

//	CWeather* getWeather() {return weather;}
//	COutput* getOutput() {return output;}
	Timer* getTime() {return time;}

	double getFirstDayOfSim() {return firstDayOfSim;}
	double getLastDayOfSim() {return lastDayOfSim;}

	double getSowingDay() {return SowingDay;}
	double getEmergenceDay() {return EmergenceDay;}
	TInitInfo getInitInfo() {return initInfo;}
//	TWeather * getCurrentWeather() {return weather;}
	void initialize();
	int run(const TWeather &, double lwpd);
	double ET_supply; //actual supply of water to plant mol m-2 (leaf) s-1
	double RootWeightFrom2DSOIL;
	double MaxRootDepth, AvailableWater;
};
#endif
