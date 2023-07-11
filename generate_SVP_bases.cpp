/*******************************************************************************
 * Program:     generate_SVP_bases

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
*******************************************************************************/
#include "common.h"

#include <iostream>
#include <fstream>
#include <string>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

using namespace std;
using namespace NTL;

int main (int argc, char * argv[]) {

	// argv[0]: executable name
	// argv[1]: m
	// argv[2]: start basis number
	// argv[3]: end basis number

	if (argc != 4) {
		cout << "Usage: <executable> <m> <start_basis_#> <end_basis_#>" << endl;
		exit(-1);
	}

	int m = atoi(argv[1]);
	int basis_number_start = atoi(argv[2]);
	int basis_number_end = atoi(argv[3]);
       
    	mat_ZZ ppOriginalBasis;
	ppOriginalBasis.SetDims(m,m);

	char filename[200];

	struct stat stats;
	stat(_INPUT_DIR_, &stats);
	if (!S_ISDIR(stats.st_mode)){
		if (!mkdir (_INPUT_DIR_, S_IRWXO|S_IRWXU|S_IRGRP)) {
			cout << "Directory created successfully." << endl;
		} else {
			cout << "Directory could not be created." << endl;
		}
	} else {
		cout << "Directory already exists." << endl;
	}

	for (int i=basis_number_start; i<=basis_number_end; i++) {

		sprintf (filename, "%s/%s_%d_%06d.txt", _INPUT_DIR_, _INPUT_FILE_PREFIX_, m, i);

          	ppOriginalBasis = generateNTLChallengeBasis (m);

		ofstream outFile (filename);
		outFile << m << endl; 

		for (int i=0; i<m; i++) {
			for (int j=0; j<m; j++) {
	      			outFile << ppOriginalBasis[i][j] << " ";
			}

			outFile << endl;
		}
		outFile.close();

	}
	return 0;

}


