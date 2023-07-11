/*******************************************************************************
 * Program:     standalone_reduction
*************************************************/
#include "common.h"

using namespace std;

#define EXECUTABLE_NAME "standalone_reduction"

#define OUTPUT_FILE_NAME "_out_"

// The main function
int main (int argc, char * argv[]) {

	// argv[0]: executable name
	// argv[1]: m
	// argv[2]: start basis number
	// argv[3]: end basis number
	// argv[4]: Algorithm

	if (argc != 5) {
		cout << "Usage: <executable> <m> <start_basis_#> <end_basis_#> <algorithm>" << endl;
		cout << "1: " << _LLL_ << endl;
		cout << "2: " << _SS_LLL_ << endl;
		cout << "3: " << _SS_GGLLL_ << endl;
		cout << "4: " << _Pot_LLL_ << endl;
		cout << "5: " << _Pot_GGLLL_ << endl;

		exit(-1);
	}
	
	int m = atoi(argv[1]);
	int n = m;
	int basis_number_start = atoi(argv[2]);
	int basis_number_end = atoi(argv[3]);
	int algo = atoi(argv[4]);

        int nIterations = basis_number_end - basis_number_start + 1;
	RR RRIterations; 
	conv(RRIterations, nIterations);

        mat_ZZ ppOriginalBasis;
	ppOriginalBasis.SetDims(m,m);
	ZZ volume;

	char filename[20000];
	char dir[10000];
	struct stat stats;

	long double delta_int = 0.999;
	long double delta_threshold_int = 0.999;
   	RR delta_threshold, delta;
   	delta_threshold = 0.999;
   	delta = 0.999;
        int a = 999;
	int b = 1000;

   	// The value of eta used in S^2: 1 - 10^(-6)
	long double eta_int = 0.999999;
   	RR eta;
   	eta = 0.999999;

   	long long int totalSwaps = 0;
 	long long int totalReductions = 0;
 	long double totalTime = 0.0;
  	RR totalV1Norm; totalV1Norm = 0;
	RR totalSVPFactor; totalSVPFactor = 0;
	RR totalRHF; totalRHF = 0;
 	RR totalOD; totalOD = 0;
 	RR totalPot; totalPot = 0;
   	RR totalSS; totalSS = 0;
   	RR V1NormInput; V1NormInput = 0; 
   	RR V1NormOutput; V1NormOutput = 0;
	RR RRSVPConstant; RRSVPConstant = 0;
	RR SVPInputFactor; SVPInputFactor = 0;
	RR SVPOutputFactor; SVPOutputFactor = 0;
   	RR RHFInput; RHFInput = 0;
   	RR RHFOutput; RHFOutput = 0;
   	RR ODInput; ODInput = 0;
   	RR ODOutput; ODOutput = 0;
   	RR PotInput; PotInput = 0;
   	RR PotOutput; PotOutput = 0;
   	RR SSInput; SSInput = 0;
   	RR SSOutput; SSOutput = 0;

      	long long int nNonReducedPairs = 0;
      	long long int nMinDeltaReductions = 0;
        
	setPrecision(m);

	sprintf (dir,"%s/%s",_PARENT_INPUT_DIR_, _INPUT_DIR_);
	cout << dir << endl;
	stat(dir, &stats);
	if (!S_ISDIR(stats.st_mode)){
		cout << "Input directory doesn't exist." << endl;
		exit(4);
	}

	switch (algo) {
		case 1:
			sprintf (dir, "%s/%s",_PARENT_OUTPUT_DIR_, _LLL_);
			break;
		case 2:
			sprintf (dir, "%s/%s",_PARENT_OUTPUT_DIR_, _SS_LLL_);
			break;
		case 3:
			sprintf (dir, "%s/%s",_PARENT_OUTPUT_DIR_, _SS_GGLLL_);
			break;
		case 4:
			sprintf (dir, "%s/%s",_PARENT_OUTPUT_DIR_, _Pot_LLL_);
			break;
		case 5:
			sprintf (dir, "%s/%s",_PARENT_OUTPUT_DIR_, _Pot_GGLLL_);
			break;
		default:
			cout << "No such basis reduction algorithm." << endl;
			exit(algo);
	}

	stat(dir, &stats);
	if (!S_ISDIR(stats.st_mode)){
		if (!mkdir (dir, S_IRWXO|S_IRWXU|S_IRGRP)) {
			cout << "Directory created successfully." << endl;
		} else {
			cout << "Directory could not be created." << endl;
			exit(5);
		}
	} else {
		cout << "Directory already exists." << endl;
	}

	for (int i=basis_number_start; i<=basis_number_end; i++) {
		
		sprintf (filename, "%s/%s/%s_%d_%06d.txt", _PARENT_INPUT_DIR_ , _INPUT_DIR_, _INPUT_FILE_PREFIX_, m, i);
		cout << filename << endl;

		ifstream inFile (filename);
		inFile >> m; 

		for (int i=0; i<m; i++) {
			for (int j=0; j<m; j++) {
	      			inFile >> ppOriginalBasis[i][j]; 
			}
		}
		
		inFile.close();

		volume = ppOriginalBasis[0][0];
		RRSVPConstant = SVP_Challenge_Factor(n, volume);

      		long long int nSwapCount = 0;
      		long long int nReductionCount = 0;
      		long long int nMinDeltaReductions = 0;
      		long double time_taken = 0.0;
      		V1NormInput = NTL_vector_norm (ppOriginalBasis[0], m);
   		SVPInputFactor = V1NormInput / RRSVPConstant;
      		RHFInput = NTLrootHermiteFactor (ppOriginalBasis, n, m);
      		ODInput = NTLorthogonalityDefect (ppOriginalBasis, n, m);
      		PotInput = NTLbasisPotential (ppOriginalBasis, n, m);
      		SSInput = NTLsquaredSum (ppOriginalBasis, n, m);

      		clock_t begin_time;
      		clock_t end_time;
		switch (algo) {
			case 1: // _LLL_
      				begin_time = clock();
      				ppOriginalBasis = LLL_mat_ZZ (ppOriginalBasis, n, m, delta, &nSwapCount, &nReductionCount);
      				end_time = clock();
				break;
			case 2: // _SS_LLL_
      				begin_time = clock();
      				ppOriginalBasis = NTL_SS_LLL (ppOriginalBasis, m, n, eta, &nSwapCount, &nReductionCount);
      				end_time = clock();
				break;
			case 3: // _SS_GGLLL_
      				begin_time = clock();
      				ppOriginalBasis = NTL_SSGGLLL (ppOriginalBasis, m, n, eta, &nSwapCount, &nReductionCount);
      				end_time = clock();
				break;
			case 4: // _Pot_LLL_
      				begin_time = clock();
      				ppOriginalBasis = NTL_PotLLL (ppOriginalBasis, m, n, delta, &nSwapCount, &nReductionCount);
      				end_time = clock();
				break;
			case 5: // _Pot_GGLLL_
      				begin_time = clock();
      				ppOriginalBasis = NTL_PotGGLLL (ppOriginalBasis, m, n, delta, &nSwapCount, &nReductionCount);
      				end_time = clock();
				break;
			default:
				cout << "No such basis reduction algorithm." << endl;
				exit(algo);
		}
      		
		time_taken = (long double) (end_time - begin_time) / CLOCKS_PER_SEC;

		totalSwaps += nSwapCount;
      		totalReductions += nReductionCount;
      		totalTime += time_taken;
      		V1NormOutput = NTL_vector_norm (ppOriginalBasis[0], m);
      		totalV1Norm += V1NormOutput;
   		SVPOutputFactor = V1NormOutput / RRSVPConstant;
		totalSVPFactor += SVPOutputFactor;
      		RHFOutput = NTLrootHermiteFactor (ppOriginalBasis, m, n);
      		totalRHF += RHFOutput;
      		ODOutput = NTLorthogonalityDefect (ppOriginalBasis, m, n);
      		totalOD += ODOutput;
      		PotOutput = NTLbasisPotential (ppOriginalBasis, m, n);
      		totalPot += PotOutput;
      		SSOutput = NTLsquaredSum (ppOriginalBasis, m, n);
      		totalSS += SSOutput;
	
		sprintf (filename, "%s/%s/%s_%d_%06d.txt", _PARENT_OUTPUT_DIR_, dir, _OUTPUT_FILE_PREFIX_, m, i);
		cout << filename << endl;

		ofstream outFile (filename);

		outFile << nSwapCount << "\t";
 		outFile << nReductionCount << "\t";
   		outFile << time_taken << "\t";
   		outFile << V1NormOutput << "\t";
   		outFile << SVPOutputFactor << "\t";
        	outFile << RHFOutput << "\t";
   		outFile << ODOutput << "\t";
   		outFile << PotOutput << "\t";
   		outFile << SSOutput << "\n";

		outFile.close();

	}

	ppOriginalBasis.kill();

	exit(0);
}
