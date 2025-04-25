#pragma once
#ifndef _WEATHER_H_
#define _WEATHER_H_

struct TWeather
{
public:
	TWeather()
	{
		year=2004,jday=1,time=0.0, CO2=370.0, airT=20.0, PFD=0.0, solRad=0.0,
		RH=50.0, wind=1.5, rain=0.0, dayLength=12.0, soilT = airT;
	}
	int year = 0;
	int jday = 0;
	double time = 0., daytime = 0.;
	double CO2 = 0., airT = 0., PFD = 0., solRad = 0., RH = 0., wind = 0., rain = 0., dayLength = 0., soilT = 0.;
	double ET_supply = 0.;
	double psil_ = 0.; //leafwater potential, MPa
	double pcrs = 0.; // rate of C supplied to roots, g h-1 plant-1
	double pcrl = 0., pcrq = 0., TotalRootWeight = 0.; //3 variables added to do mass balance checks with root allocation in 2DSOIL as per Dennis
	double swdf1 = 0., swdf2 = 0.; //soil water deficiency factors 1 and 2 from SIMGUI; swdf1 midiifes photosynthesis; swdf2 modifies leaf expansion
	double swpartition = 0.; //ranges from 0 to 1, with 1 meaning roots can get 100% priority of canopy biomass during h2o deficiency
	double trwu_simpotato = 0.; //total potential root water extraction rate from soil for hour, as in SIMPOTATO
	double psism_simpotato = 0.; //average soil water potential from slab for hourly potential extraction from 2DSOIL carbon partitioning
	double psilt_simpotato = 0.; //threshold bulk leaf water potential below which shoot growth gets restricted from 2DSOIL carbon partitioning
	double MaxRootDepth = 0.,ThetaAvailRZ = 0.; //for autoirrigate functions from MAIZSIM
	int DailyOutput = 0, HourlyOutput = 0;
};
#endif
