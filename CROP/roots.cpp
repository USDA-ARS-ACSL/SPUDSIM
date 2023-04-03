#include "stdafx.h"
#include "roots.h"

CRoots::CRoots(void):COrgan()
{
}

CRoots::CRoots(const TInitInfo& info)
{
	this->set_CH2O(0);
	this->set_currentNitrogenAmount(0);
	this->set_drymass(0);
	this->set_drymass_at_emergence(0);
	this->set_maintCoeff(0.0206); //gCH2O g-1DM day-1 at 25C as in Ng and Loomis
	set_actualNitrogenRatio(0.0215);
	set_optimumNitrogenRatio(0);
	set_minimumNitrogenRatio(0);
	set_currentNitrogenAmount(0);
	set_optimumNitrogenRatio(0.0215);//SUBSTOR
	set_actualNitrogenRatio(0.0215);//SUBSTOR
	set_minimumNitrogenRatio(0.014);//SUBSTOR
	//this->set_initialized(); //just a flag for 2DSOIL to use to pass up initial carbon content
}

void CRoots::set_EmergenceData(double mass)
/***********************************************/
/* Method needed if plant starts at emergence*/
/* Need to estimate initial root mass and N content*/
/* In SIMPOTATO, root growth at emergence is assumed equal to
/*		the weigth of the sprout underground.  
/* In SPUDSIM, sprout weight currently not simulated
/* For simplicity, set root mass equal to initial stem mass at emergence, 0.2 g CH2O plant-1
/* However, this initial root mass is set in the ElemGeomFile from 2DSOIL
/**********************************************/
{
	this->set_drymass(mass);
	this->set_CH2O(mass);
	this->set_drymass_at_emergence(mass);
	this->set_currentNitrogenAmount(mass * this->get_actualNitrogenRatio());
}


CRoots::~CRoots(void)
{
}