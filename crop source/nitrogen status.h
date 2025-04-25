#pragma once
#ifndef _INITINFO_H_
#define _INITINFO_H_
#define MINUTESPERDAY (24*60);

// Contains whole plant nitrogen information and stress factors

struct TNitrogenStatus
{
public:
	TNitrogenStatus()
	{
	}
	double availableNitrogenRatio; //From SIMPOTATO (AVAILN), current amoutn of N available in plant, g N g plant-1, is this same as TotalNitrogen?
	// below are state N variables //
	double optimumleafNitrogenRatio; //From SIMPOTATO, optimum or critical leaf N, in g N g leaf-1 required to achieve maximum growth rate
	double minimumleafNitrogenRatio; //From SIMPOTATO, minimum leaf N
	double actualleafNitrogenRatio; //From SIMPOTATO, actual leaf N
	double optimumstemNitrogenRatio, minimumstemNitrogenRatio, actualstemNitrogenRatio; // as above, for stem
	double optimumtuberNitrogenRatio, minimumtuberNitrogenRatio, actualtuberNitrogenRatio; //", for tubers
	double optimumrootNitrogenRatio, minimumrootNitrogenRatio, actualrootNitrogenRatio;// as above, for roots
	double leafNitrogenAmount, stemNitrogenAmount, rootNitrogenAmount; //total N content in all leaves, stems, or roots (g N plant-1)
	double leafNitrogenbuffer, stemNitrogenbuffer; //SIMPOTATO (leafNitrogenbuffer, BUFST) - N above min or optimum level in leaves or stems (g N plant-1)
	double minimumvegetativeNitrogen; // g N g-1, SIMPOTATO
	// demand N variables
	double totalNitrogendemand; //From SIMPOTATO, total N demand (ANDEM), initially determined for existing biomass and potential new biomass, g N plant-1
	double leafNitrogendemand; // From SIMPOTATO (TNDEM), g N plant-1
	double stemNitrogendemand; // From SIMPOTATO (stemNitrogendemand)
	double rootNitrogendemand; // From SIMPOTATO (RNDEM)
	// whole plant N stress variables
	double Nitrogenstressfactor; // From SIMPOTATO (NFAC), 0 to 1
	double Nitrogendeficiencyone; // From SIMPOTATO (NDEF1), 0 to 1; reduces whole plant p.s. rate
	double Nitrogendeficiencytwo; // From SIMPOTATO (NDEF2), 0 to 1; reduces whole plant leaf expansion, increases senescence
	double Nitrogendeficiencythree; // From SIMPOTATO (NDEF3), 0 to 1; ?
	// crop.dll / 2dspudsim.for linkages (legacy from MAIZSIM - these are needed to keep same structure as corn model)
	double HourlyNitrogenDemand; //g N plant-1 h-1 needed for non-stress conditions at current time-step
	double CumulativeNitrogenDemand; //g N plant-1 season-1, amount of N needed since emergence
	double HourlyNitrogenSoilUptake; //g N plant-1 h-1, pulled from soil at current time-step
	double CumulativeNitrogenSoilUptake; //g N plant-1 season-1, pulled from soil since emergence
};
#endif