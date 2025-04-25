#pragma once
#include "organ.h"

class CTubers :
	public COrgan
{
public:
	CTubers(void);
	bool isInitiated() {return initiated;}
	void setInitiated() {initiated = true; this->set_currentNitrogenAmount(0);}
	double get_tubmax() {return tubmax;}
	double get_feedback(double Tubfraction);
	double get_maintResp(double age);
	double Growth(CDevelopment* dv, double Tlag, const TInitInfo info);
	virtual ~CTubers(void);
private:
	double etgt = 0., tubmax = 0.; //max potential tuber growth factor in simgui
	bool initiated = 0; //true if tubers initiation has started
	double dblSolids = 0.; //percent solids in tubers as per simgui
	double maintCoeff = 0.; // maintenance respiration coefficient, g ch2o g-1 dm-1 d-1
	double maintResp = 0.; //g CH2O d-1 lost to maintenance respiration
	//double dblPartub; //portion of tuber capable of growth as per simgui
};