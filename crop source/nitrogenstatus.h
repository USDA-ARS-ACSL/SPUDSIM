#pragma once
#ifndef _NITROGENSTATUS_H_
#define _NITROGENSTATUS_H_

struct TNitrogenStatus //Provide structure to pass all whole plant N related information to different classes as needed
{
public:
	TNitrogenStatus()
	{
		availableNitrogenRatio = 0; //5-24-2021
		availableNitrogen=0;optimumleafNitrogenRatio=0;minimumleafNitrogenRatio=0;actualleafNitrogenRatio=0;
		optimumstemNitrogenRatio=0;minimumstemNitrogenRatio=0;actualstemNitrogenRatio=0;
		optimumtuberNitrogenRatio=0;minimumtuberNitrogenRatio=0;actualtuberNitrogenRatio=0;
		optimumrootNitrogenRatio=0;minimumrootNitrogenRatio=0;actualrootNitrogenRatio=0;
		leafNitrogenAmount=0;stemNitrogenAmount=0;rootNitrogenAmount=0;tuberNitrogenAmount=0;
		leafNitrogenbuffer=0;stemNitrogenbuffer=0;tuberNitrogenbuffer=0;
		minimumvegetativeNitrogenRatio=0;
		totalNitrogendemand=0;leafNitrogendemand=0;stemNitrogendemand=0;rootNitrogendemand=0;tuberNitrogendemand=0;
		Nitrogenstressfactor=0;Nitrogendeficiencyone=0;Nitrogendeficiencytwo=0;Nitrogendeficiencythree=0;
		HourlyNitrogenDemand = 0;CumulativeNitrogenDemand = 0; HourlyNitrogenSoilUptake = 0; CumulativeNitrogenSoilUptake = 0;
	}
	double availableNitrogenRatio = 0;//5-24-2021
	double availableNitrogen = 0; //From SIMPOTATO (AVAILN), current amoutn of N available to the plant, g N plant-1, is this same as TotalNitrogen?
	// below are state N variables //
	double optimumleafNitrogenRatio = 0; //From SIMPOTATO, optimum or critical leaf N, in g N g leaf-1 required to achieve maximum growth rate
	double minimumleafNitrogenRatio = 0; //From SIMPOTATO, minimum leaf N
	double actualleafNitrogenRatio = 0; //From SIMPOTATO, actual leaf N
	double optimumstemNitrogenRatio = 0, minimumstemNitrogenRatio = 0, actualstemNitrogenRatio = 0; // as above, for stem
	double optimumtuberNitrogenRatio = 0, minimumtuberNitrogenRatio = 0, actualtuberNitrogenRatio = 0; //", for tubers
	double optimumrootNitrogenRatio = 0, minimumrootNitrogenRatio = 0, actualrootNitrogenRatio = 0;// as above, for roots
	double leafNitrogenAmount = 0, stemNitrogenAmount = 0, rootNitrogenAmount = 0, tuberNitrogenAmount = 0; //total N content in all leaves, stems, or roots (g N plant-1)
	double leafNitrogenbuffer = 0, stemNitrogenbuffer = 0, tuberNitrogenbuffer = 0; //SIMPOTATO (leafNitrogenbuffer, BUFST) - N above min or optimum level in leaves or stems (g N plant-1)
	double minimumvegetativeNitrogenRatio = 0; // g N g-1, SIMPOTATO
	// demand N variables
	double totalNitrogendemand = 0; //From SIMPOTATO, total N demand (ANDEM), initially determined for existing biomass and potential new biomass, g N plant-1
	double leafNitrogendemand = 0; // From SIMPOTATO (TNDEM), g N plant-1
	double stemNitrogendemand = 0; // From SIMPOTATO (stemNitrogendemand)
	double rootNitrogendemand = 0; // From SIMPOTATO (RNDEM)
	double tuberNitrogendemand = 0; // From SIMPOTATO (TUBDEM)
	// whole plant N stress variables
	double Nitrogenstressfactor = 0; // From SIMPOTATO (NFAC), 0 to 1
	double Nitrogendeficiencyone = 0; // From SIMPOTATO (NDEF1), 0 to 1; reduces whole plant p.s. rate
	double Nitrogendeficiencytwo = 0; // From SIMPOTATO (NDEF2), 0 to 1; reduces whole plant leaf expansion, increases senescence
	double Nitrogendeficiencythree = 0; // From SIMPOTATO (NDEF3), 0 to 1; ?
	// crop.dll / 2dspudsim.for linkages (legacy from MAIZSIM - these are needed to keep same structure as corn model)
	double HourlyNitrogenDemand = 0; //g N plant-1 h-1 needed for non-stress conditions at current time-step
	double CumulativeNitrogenDemand = 0; //g N plant-1 season-1, amount of N needed since emergence
	double HourlyNitrogenSoilUptake = 0; //g N plant-1 h-1, pulled from soil at current time-step
	double CumulativeNitrogenSoilUptake = 0; //g N plant-1 season-1, pulled from soil since emergence

};
#endif