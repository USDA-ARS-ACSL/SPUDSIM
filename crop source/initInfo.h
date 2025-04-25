#pragma once
#ifndef _INITINFO_H_
#define _INITINFO_H_
#define MINUTESPERDAY (24.0*60.0)
#ifndef FLOAT_EQ
#define EPSILON 0.001   // floating point comparison tolerance
#define FLOAT_EQ(x,v) (((v - EPSILON) < x) && (x <( v + EPSILON)))
#endif

// Contains planting, some field, and cultivar information

struct TInitInfo
{
public:
	TInitInfo()
	{
	}
		char description[255];
		char cultivar[255];
		double plantDensity = 0; // plants m-2 (tubers m-2 with single mainstem per tuber assumed)
		double A1 = 0, A2 = 0, A3 = 0, A4 = 0, A5 = 0, A6 = 0, A7 = 0, A8 = 0, A9 = 0, A10 = 0; // genetic info for tuber initiation
		double G1 = 0, G2 = 0, G3 = 0, G4 = 0; // genetic info for growth and leaf expansion
		double latitude = 0, longitude = 0, altitude = 0;
		int year = 0;
		double beginDay = 0, sowingDay = 0, emergeDay = 0, endDay = 0;
		double Seedreserve = 0; //CHO in tuber seedpiece at planting, g 
		double Plantingdepth = 0, Sproutlength = 0; //(depth in cm, length in cm)
		double timeStep = 0;  //min
		double outtimeStep = 0; //min
		double CO2 = 400; //daily or hourly CO2 value in ppm
		int bigleaf = 0; // = 0 for single leaf model, 1 for big leaf, SIMPOTATO canopy growth model
		int Nitrogen_stress_off = 0; // 1 for no stress
		int Water_stress_off = 0; //1 for no stress
		int Water_stress_simulation_type = 0; //0 for AD version; 1 for SIMPOTATO (micrometeorolgoical) apporach
};
#endif