#pragma once
#include "initinfo.h"
class CThermalTime
{
public:
	CThermalTime(void);
	virtual ~CThermalTime(void);
	double get_Tcur() {return Tcur;}
	double get_Tbase() {return Tbase;}
	double get_Topt() {return Topt;}
	double get_Tmax() {return Tmax;}
	double get_sum() {return sum;}
	double get_dTmpr() {return dTmpr;}
	double get_timeStep() {return timeStep;}
	double get_actualAge() {return actualAge;}
	double Agefraction(const TInitInfo info, double age, double T, double duration, double Tlagleaf);
	double beta_fn(double t, double R_max, double t_m, double t_e);
	double Tfactor(const TInitInfo info, double T, double Tlagleaf);
	double Cage(const TInitInfo info); // method calcs. leaf chronological age
	double Page(const TInitInfo info, double T, double Tlagleaf); //method calcs. leaf physiological age
	void set_Tcur(double x) {Tcur = x;}
	void set_Tbase(double x) {Tbase = x;}
	void set_Topt(double x) {Topt=x;}
	void set_Tmax(double x) {Tmax=x;}
	void set_timeStep(double x) {timeStep=x;}
	void set_temperatures(double Tb, double To, double Tm) {Tbase = Tb; Topt = To; Tmax = Tm;}
	void add(double x);
    void update(double Tmpr, double step);
	void initialize(double step);


private:
	TInitInfo initInfo;
	double Tcur = 0.; //previous 24 H air temperature surrounding organ - leaf expansion model (Fleisher and Timlin 2006)
	double Tbase = 0.; //minimum air T below which no leaf expansion occurs (Fleisher and Timlin, 2006)
	double Topt = 0.; //optimum T at which max expansion occurs "
	double Tmax = 0.; //max T at whihc no xpansion occurs "
	double sum = 0.; //old MAIZEIM variable - not used?
	double dTmpr = 0.; //old MAIZEIM variable - not used?
	double timeStep = 0.; //old MAIZEIM variable - not used?
	double actualAge = 0.; //old MAIZEIM variable - not used?

};
