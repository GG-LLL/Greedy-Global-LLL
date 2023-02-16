/*******************************************************************************
 * Program:     collate
 ********************************************************************************/
#include "common.h"

// The main function
int main (int argc, char * argv[]) {

	// argv[0]: executable name
	// argv[1]: m
	// argv[2]: end basis number

	if (argc != 3) {
		cout << "Usage: <executable> <m> <end_basis_#>" << endl;
		exit(-1);
	}
	
	int m = atoi(argv[1]);
	int basis_number_end = atoi(argv[2]);

    int nIterations = basis_number_end;
	RR RRIterations; conv(RRIterations, nIterations);

	char filename[20000];
	char dir[15000];
	char algo_name[10000];
	struct stat stats;

   	long long int Swaps = 0;
 	long long int Reductions = 0;
 	long double Time = 0.0;
  	RR V1Norm; V1Norm = 0;
	RR SVPFactor; SVPFactor = 0;
	RR RHF; RHF = 0;
 	RR OD; OD = 0;
 	RR Pot; Pot = 0;
   	RR SS; SS = 0;
   
	long long int totalSwaps = 0;
 	long long int totalReductions = 0;
 	long double totalTime = 0.0;
  	RR totalV1Norm; totalV1Norm = 0;
	RR totalSVPFactor; totalSVPFactor = 0;
	RR totalRHF; totalRHF = 0;
 	RR totalOD; totalOD = 0;
 	RR totalPot; totalPot = 0;
   	RR totalSS; totalSS = 0;
	
	map<string, long double> mapAverageSwaps;
	map<string, long double> mapAverageReductions;
	map<string, long double> mapAverageTime;
   	map<string, double> mapAverageV1Norm;
   	map<string, double> mapAverageSVPFactor;
   	map<string, double> mapAverageRHF;
   	map<string, double> mapAverageOD;
   	map<string, double> mapAveragePot;
   	map<string, double> mapAverageSS;
        
	int algoNumbers[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9};

	for (int i : algoNumbers) {

		switch (i) {
			case 1:
				sprintf (algo_name, "%s", _LLL_);
				break;
			case 2:
				sprintf (algo_name, "%s", _SS_LLL_);
				break;
			case 3:
				sprintf (algo_name, "%s", _SS_GGLLL_);
				break;
			case 4:
				sprintf (algo_name, "%s", _Pot_LLL_);
				break;
			case 5:
				sprintf (algo_name, "%s", _Pot_GGLLL_);
				break;
			case 6:
				sprintf (algo_name, "%s", _BKZ_08_);
				break;
			case 7:
				sprintf (algo_name, "%s", _BKZ_10_);
				break;
			case 8:
				sprintf (algo_name, "%s", _BKZ_12_);
				break;
			case 9:
				sprintf (algo_name, "%s", _BKZ_20_);
				break;
			default:
				cout << "No such basis reduction algorithm." << endl;
				exit(i);
		}

		sprintf (dir, "../SVP_Outputs/%s", algo_name);
		cout << dir << endl;
		stat(dir, &stats);
		if (!S_ISDIR(stats.st_mode)){
			cout << "Input directory doesn't exist." << endl;
			exit(4);
		}
   	
		totalSwaps = 0;
 		totalReductions = 0;
 		totalTime = 0.0;
  		totalV1Norm = 0;
		totalSVPFactor = 0;
		totalRHF = 0;
 		totalOD = 0;
 		totalPot = 0;
   		totalSS = 0;

		for (int i=1; i<=basis_number_end; i++) {
			sprintf (filename, "%s/out_%d_%06d.txt", dir, m, i);
			cout << filename << endl;
	
			ifstream inFile (filename);
	
   			inFile >> Swaps;
 			inFile >> Reductions;
 			inFile >> Time;
  			inFile >> V1Norm;
			inFile >> SVPFactor;
			inFile >> RHF;
 			inFile >> OD;
 			inFile >> Pot;
   			inFile >> SS;
   	
			totalSwaps += Swaps;
 			totalReductions += Reductions;
 			totalTime += Time;
  			totalV1Norm += V1Norm;
			totalSVPFactor += SVPFactor;
			totalRHF += RHF;
 			totalOD += OD;
 			totalPot += Pot;
   			totalSS += SS;
			
			inFile.close();
        
		}

		if (totalSwaps != 0) {		
			mapAverageSwaps[algo_name] = (long double) (totalSwaps) / (long double) (nIterations);
		}
		
		if (totalReductions != 0) {
			mapAverageReductions[algo_name] = (long double) (totalReductions) / (long double) (nIterations);
		}
		mapAverageTime[algo_name] = (long double) (totalTime) / (long double) (nIterations);
   		RR RRAverageV1Norm; RRAverageV1Norm =  (totalV1Norm) / (RRIterations);
   		conv(mapAverageV1Norm[algo_name], RRAverageV1Norm);
		RR RRAverageSVPFactor; RRAverageSVPFactor = (totalSVPFactor) / (RRIterations);
   		conv(mapAverageSVPFactor[algo_name], RRAverageSVPFactor); 
        RR RRAverageRHF; RRAverageRHF = (totalRHF) / (RRIterations);
        conv(mapAverageRHF[algo_name], RRAverageRHF);
   		RR RRAverageOD; RRAverageOD = (totalOD) / (RRIterations);
   		conv(mapAverageOD[algo_name], RRAverageOD);
   		RR RRAveragePot; RRAveragePot = (totalPot) / (RRIterations);
   		conv(mapAveragePot[algo_name], RRAveragePot);
   		RR RRAverageSS; RRAverageSS =  (totalSS) / (RRIterations);
   		conv(mapAverageSS[algo_name], RRAverageSS);
		
	}

   	printf("\n--------------------------------------------\n");
   	printf("--------------------------------------------\n");
   	printf ("Averages for %d basis vectors:", nIterations);
   	printf("\n--------------------------------------------\n");
   	printf("--------------------------------------------\n");
   	
	printf("\n--------------\nSWAPS\n--------------\n");
   	map_sort_long_double_fix(mapAverageSwaps);
   	printf("\n-------------------\nREDUCTIONS\n-------------------\n");
   	map_sort_long_double_fix(mapAverageReductions);
    printf ("\n-------------\nTIME\n-------------\n");
   	map_sort_long_double_fix(mapAverageTime);
   	printf("----------------------\nLENGTH OF FIRST VECTOR\n----------------------\n");
   	map_sort_double_fix(mapAverageV1Norm);
   	printf("\n----------------------\nSVP CHALLENGE CONSTANT\n----------------------\n");
   	map_sort_double_fix(mapAverageSVPFactor);
   	printf("\n-------------------\nROOT HERMITE FACTOR\n-------------------\n");
   	map_sort_double_fix(mapAverageRHF);
   	printf("\n--------------------\nORTHOGONALITY DEFECT\n--------------------\n");
   	map_sort_double_fix(mapAverageOD);
   	printf("\n---------\nPOTENTIAL\n---------\n");
   	map_sort_double_fix(mapAveragePot);
   	printf("\n-----------\nSQUARED SUM\n-----------\n");
   	map_sort_double_fix(mapAverageSS);
	
	exit(0);

}
