/***************************************************
 * Program:     preprocessed_reduction
 
 * Copyright (C) 2023; Jack Moyler, Sanjay Bhattacherjee.

 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.

 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.

 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*************************************************/
#include "common.h"

using namespace std;

#define EXECUTABLE_NAME "preprocessed_reduction"

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
		cout << "1: " << _SS_LLL_ << endl;
		cout << "2: " << _SS_GGLLL_ << endl;
		cout << "3: " << _Pot_LLL_ << endl;
		cout << "4: " << _Pot_GGLLL_ << endl;

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

	int ** ppBasisInt;
   	ppBasisInt = (int **) calloc (n, sizeof(int *));
   	for (int i=0; i<n; i++) {
      		ppBasisInt[i] = (int *) calloc (m, sizeof(int));
	}

        mat_ZZ ppBasisZZ;
	ppBasisZZ.SetDims(m,m);
	ZZ volume;

	char filename[20000];
	char dir[10000];
	struct stat stats;

	long double delta_threshold = 0.999;

   	// The value of eta used in S^2: 1 - 10^(-6)
	long double eta = 0.999999;

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
        
	sprintf (dir,"%s/%s",_PARENT_INPUT_DIR_, _INPUT_DIR_FPLLL_RED_);
	cout << dir << endl;
	stat(dir, &stats);
	if (!S_ISDIR(stats.st_mode)){
		cout << "Input directory doesn't exist." << endl;
		exit(4);
	}

	switch (algo) {
		case 1:
			sprintf (dir, "%s/%s",_PARENT_OUTPUT_DIR_PREPROCESSED_, _SS_LLL_);
			break;
		case 2:
			sprintf (dir, "%s/%s",_PARENT_OUTPUT_DIR_PREPROCESSED_, _SS_GGLLL_);
			break;
		case 3:
			sprintf (dir, "%s/%s",_PARENT_OUTPUT_DIR_PREPROCESSED_, _Pot_LLL_);
			break;
		case 4:
			sprintf (dir, "%s/%s",_PARENT_OUTPUT_DIR_PREPROCESSED_, _Pot_GGLLL_);
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
		
		sprintf (filename, "%s/%s/%s_%d_%06d.txt", _PARENT_INPUT_DIR_ , _INPUT_DIR_FPLLL_RED_, _INPUT_FILE_PREFIX_, m, i);

		char c;
		ifstream inFile (filename);
		inFile >> c;
		for (int i=0; i<m; i++) {
			inFile >> c;
			for (int j=0; j<m; j++) {
				inFile >> ppBasisInt[i][j];
			}
			inFile >> c;
		}
		inFile >> c;
		inFile.close();	
		
		ppBasisZZ = copyBasisToNTL (ppBasisInt, m, n);

		volume = abs(determinant(ppBasisZZ, 1));
		RRSVPConstant = SVP_Challenge_Factor(n, volume);

      		long long int nSwapCount = 0;
      		long long int nReductionCount = 0;
      		long long int nMinDeltaReductions = 0;
      		long double time_taken = 0.0;
      		V1NormInput = NTL_vector_norm (ppBasisZZ[0], m);
   		SVPInputFactor = V1NormInput / RRSVPConstant;
      		RHFInput = NTLrootHermiteFactor (ppBasisZZ, n, m);
      		ODInput = NTLorthogonalityDefect (ppBasisZZ, n, m);
      		PotInput = NTLbasisPotential (ppBasisZZ, n, m);
      		SSInput = NTLsquaredSum (ppBasisZZ, n, m);

      		clock_t begin_time;
      		clock_t end_time;
		switch (algo) {
			case 1: // _SS_LLL_
      				begin_time = clock();
      				ppBasisInt = SS_LLL_std (ppBasisInt, m, n, eta, &nSwapCount, &nReductionCount);
      				end_time = clock();
				break;
			case 2: // _SS_GGLLL_
      				begin_time = clock();
      				ppBasisInt = SSGGLLL_std (ppBasisInt, m, n, eta, &nSwapCount, &nReductionCount);
      				end_time = clock();
				break;
			case 3: // _Pot_LLL_
      				begin_time = clock();
      				ppBasisInt = Pot_LLL_std (ppBasisInt, m, n, delta_threshold, &nSwapCount, &nReductionCount);
      				end_time = clock();
				break;
			case 4: // _Pot_GGLLL_
      				begin_time = clock();
      				ppBasisInt = PotGGLLL_std (ppBasisInt, m, n, delta_threshold, &nSwapCount, &nReductionCount);
      				end_time = clock();
				break;
			default:
				cout << "No such basis reduction algorithm." << endl;
				exit(algo);
		}
      		
		time_taken = (long double) (end_time - begin_time) / CLOCKS_PER_SEC;
		
		ppBasisZZ = copyBasisToNTL (ppBasisInt, m, n);

		totalSwaps += nSwapCount;
      		totalReductions += nReductionCount;
      		totalTime += time_taken;

      		V1NormOutput = NTL_vector_norm (ppBasisZZ[0], m);
      		totalV1Norm += V1NormOutput;
   		SVPOutputFactor = V1NormOutput / RRSVPConstant;
		totalSVPFactor += SVPOutputFactor;
      		RHFOutput = NTLrootHermiteFactor (ppBasisZZ, m, n);
      		totalRHF += RHFOutput;
      		ODOutput = NTLorthogonalityDefect (ppBasisZZ, m, n);
      		totalOD += ODOutput;
      		PotOutput = NTLbasisPotential (ppBasisZZ, m, n);
      		totalPot += PotOutput;
      		SSOutput = NTLsquaredSum (ppBasisZZ, m, n);
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

	ppBasisZZ.kill();
	deleteBasis (ppBasisInt, n);

	exit(0);
}
