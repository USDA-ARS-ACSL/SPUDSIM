#include "stdafx.h"
#include "math.h"
#include "nodalunit.h"
#include "weather.h"

CNodalUnit::CNodalUnit()
{
	leaf = NULL;
	stem = NULL;
	dblLAR = 0.;
	dblBAR = 0.;
	initiated = growing = aging = terminated = false;
}
void CNodalUnit::initialize(const TInitInfo info, int n, int x, int y, int z, CDevelopment * dv)
{
	id.val = n; type.val = x; location.val = y; node.val = z;
	leaf = new CLeaf(info, n, z);
	stem = new CStem(info, n);
	leaf->initialize();
	stem->initialize();
	initiated = true;
	sunlit = true;
	// Stagger age of 1st four leaves that initiate by 0.5 days to smooth out growth at emergence//
	/*if (n < 4)
	{
		this->get_leaf()->set_age(1.-n*0.25);
		this->get_stem()->set_age(1.-n*0.25);
		this->get_leaf()->set_physAge(1.-n*0.25);
		this->get_stem()->set_physAge(1.-n*0.25);
	}*/
}

CNodalUnit::~CNodalUnit()
{
	if (leaf !=NULL) delete leaf;
	if (stem !=NULL) delete stem;
//	delete sheath;
//	delete internode;
}

void CNodalUnit::set_LAR()
{
	//if (type.val==0) dblLAR -= 1; //reset appearance rate counter to 0 for all branches after new leaf except mainstem
	//else
		dblLAR = 0.;
}

void CNodalUnit::update(int iCur, CDevelopment * dv, const TInitInfo info, const TWeather& weather, const TNitrogenStatus nitrogen, double Tdaylag, double Sradave, double gLAI, int potential, double CH2O, double lwpd, double swdf1, double C_pool_room)
{
	/*********************************************/
	/* Determines aging, potential, and actual growth for nodal units*/
	/* Also determines senescence.
	/* Method is by the PLANT class Update_Nodes method.  the first time determines potential expansion
	/* Note if bigleaf model used, there is also an estimate for stem senescence rate
	/********************************************/

	double Tlagleaf = Tdaylag;
	double area = leaf->get_greenArea();
	if (Tlagleaf < -20.) Tlagleaf = 20.; //set default T at first day of simulation if no previous T exists
	//Tlagleaf = 20; // use this line only to test leaf expansion at constant T
	//each stem keeps its own leaf appearance rate counter stored at node 0
	if ((type.val == 0) && (isInitiated())) //main stem node
	{
		if (node.val == 0) dblLAR += dv->get_LAR(); //tabulate only on node 0 of each stem
		if (node.val == 0) dblBAR += dv->get_BAR(); //tabulate branch initiation rate on node 0 of ms stem only; assume same T response as leaf
	}
	if ((type.val == 1) && (isInitiated())) //basal stem node
	{
		if (node.val == 0) dblLAR += dv->get_LAR();
	}
	if ((type.val == 2) && (isInitiated())) //apical stem nodes
	{
		if (node.val == 0) dblLAR += dv->get_LAR();
	}
	if ((fmod(iCur+1,((24.0*60.0)/info.timeStep)) < 0.001))//set leaf SLA at beginning of day
	{
		leaf->set_SLA(info, Sradave, gLAI);
	}
	leaf->update(iCur,dv,info,weather,nitrogen, Tlagleaf,potential,CH2O,lwpd,swdf1,C_pool_room);
	if (info.bigleaf == 1) stem->update(area,iCur,dv,info,weather,Tlagleaf);
	if (potential==1)
	{
		mass = leaf->get_drymass() + stem->get_drymass();
	}
}


