#ifndef __NODALUNIT_H__
#define __NODALUNIT_H__ 

#include "organ.h"
#include "leaf.h"
#include "stem.h"
#include "development.h"


struct TRank 
{
public:
	TRank() {val = 0;}
	int val;  
};

class CNodalUnit
{
public:
	CNodalUnit();
	~CNodalUnit();
	//int get_rank() {return rank;}
	bool isInitiated() {return initiated;}
	bool isGrowing() {return growing;}
	bool isProlific() {return prolific;}
	bool isAging() {return aging;}
	bool isTerminated() {return terminated;}
	bool isSunlit() {return sunlit;}
	
    CLeaf * get_leaf() {return leaf;}
	CStem * get_stem() {return stem;}
//	CSheath * get_sheath() {return sheath;}
//	CInternode * get_internode() {return internode;}
	void set_leaf(CLeaf * x) {leaf=x;}
	void set_stem(CStem * x) {stem=x;}
	void set_Sunlit(bool x) {sunlit = x;}
//	void set_sheath(CSheath * x) {sheath=x;}
//	void set_internode(CInternode * x) {internode=x;}
	void update(int iCur, CDevelopment *, const TInitInfo, const TWeather &, const TNitrogenStatus, double Tdaylag, double Sradave, double gLAI, int potential, double CH2O, double lwpd, double swdf1, double C_pool_room);
	void initialize(const TInitInfo, int, int, int, int, CDevelopment * dv);

	double get_leafLength(int rank);
	double get_LAR() {return dblLAR;}
	double get_BAR() {return dblBAR;}
	void set_LAR();
	void set_BAR() {dblBAR = 0;}
	int get_node(){return node.val;}
	int get_location(){return location.val;}
	int get_type(){return type.val;}
	bool get_sunlit(){return sunlit;}
	

private:
	TInitInfo initInfo;
	//int rank;
	TRank id;  // overall id, 0 through ... based on order of instantiation
	TRank type;  //(0,1,2 for main, basal, apical stem)
	TRank location;  //(leaf position from mainstem that branch originates from)-describes branch, not node
	TRank node;  //(location of unit along current stem as counted from zero at origin) - describes node

	bool initiated, growing, prolific, aging, terminated;
	CLeaf * leaf;
	CStem * stem;
	double mass = 0;
	double dblLAR = 0;
	double dblBAR = 0;
	bool sunlit = 0; // 0 for sunny, 1 for shaded
};
#endif

