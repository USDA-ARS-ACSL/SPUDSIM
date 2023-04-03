#pragma once
#include "organ.h"

class CRoots :
	public COrgan
{
public:
	CRoots(void);
	CRoots(const TInitInfo& info);
	void set_EmergenceData(double);
	//void set_initialized() {initialized = true;}
	virtual ~CRoots(void);

private:
	//bool initialized; //set to true once emerged in 2DSOIL - should allow 2DSOIL to do simulations prior to emergence of crop
};
