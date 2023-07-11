/*******************************************************************************
 * Header File: common.cpp
 *              Contains definitions of all common functions required
 *              for lattice-based cryptography
*******************************************************************************/
#include "common.h"

// #define __DEBUG__

int setPrecision(int m) {
      	if (m <= 40) {         
        	RR::SetPrecision(1100);
   
         	#ifdef NTL_ZZ_NBITS 
         	#undef NTL_ZZ_NBITS
         	#endif
         	#define NTL_ZZ_NBITS (400)

      	} else if (m <= 45) {
         	RR::SetPrecision(1250);

         	#ifdef NTL_ZZ_NBITS
         	#undef NTL_ZZ_NBITS
         	#endif
         	#define NTL_ZZ_NBITS (450)

      	} else if (m <= 50) {
         	RR::SetPrecision(1325);
         	#ifdef NTL_ZZ_NBITS
         	#undef NTL_ZZ_NBITS
         	#endif
         	#define NTL_ZZ_NBITS (500)
      	} else if (m <= 55) {
         	RR::SetPrecision(1375);
 
         	#ifdef NTL_ZZ_NBITS
         	#undef NTL_ZZ_NBITS
         	#endif
         	#define NTL_ZZ_NBITS (550)
      	} else if (m <= 60) {
         	RR::SetPrecision(1550);

         	#ifdef NTL_ZZ_NBITS
         	#undef NTL_ZZ_NBITS
         	#endif
         	#define NTL_ZZ_NBITS (600)
      	} else if (m <= 65) {
         	RR::SetPrecision(1650);

         	#ifdef NTL_ZZ_NBITS
         	#undef NTL_ZZ_NBITS
         	#endif
         	#define NTL_ZZ_NBITS (650)
      	} else if (m <= 70) {
         	RR::SetPrecision(1750);

         	#ifdef NTL_ZZ_NBITS
         	#undef NTL_ZZ_NBITS
         	#endif
         	#define NTL_ZZ_NBITS (700)
      	} else if (m <= 75) {
         	RR::SetPrecision(1900);

         	#ifdef NTL_ZZ_NBITS
         	#undef NTL_ZZ_NBITS
         	#endif
         	#define NTL_ZZ_NBITS (750)
      	} else if (m <= 80) {
         	RR::SetPrecision(2000);

         	#ifdef NTL_ZZ_NBITS
         	#undef NTL_ZZ_NBITS
         	#endif
         	#define NTL_ZZ_NBITS (800)
      	} else if (m <= 85) {
         	RR::SetPrecision(2200);

         	#ifdef NTL_ZZ_NBITS
         	#undef NTL_ZZ_NBITS
         	#endif
         	#define NTL_ZZ_NBITS (850)
      	} else if (m <= 90) {
         	RR::SetPrecision(2300);

         	#ifdef NTL_ZZ_NBITS
         	#undef NTL_ZZ_NBITS
         	#endif
         	#define NTL_ZZ_NBITS (900)
      	} else if (m <= 100) {
         	RR::SetPrecision(3000);

         	#ifdef NTL_ZZ_NBITS
         	#undef NTL_ZZ_NBITS
         	#endif
         	#define NTL_ZZ_NBITS (1000)
      	} else if (m <= 105) {
         	RR::SetPrecision(3200);
   
         	#ifdef NTL_ZZ_NBITS
         	#undef NTL_ZZ_NBITS
         	#endif
         	#define NTL_ZZ_NBITS (1100)
      	} else if (m <= 110) {
         	RR::SetPrecision(3400);

         	#ifdef NTL_ZZ_NBITS
         	#undef NTL_ZZ_NBITS
         	#endif
         	#define NTL_ZZ_NBITS (1100)
      	} else if (m <= 115) {
         	RR::SetPrecision(3600);
   
         	#ifdef NTL_ZZ_NBITS
         	#undef NTL_ZZ_NBITS
         	#endif
         	#define NTL_ZZ_NBITS (1150)
      	} else if (m <= 120) {
         	RR::SetPrecision(3800);

         	#ifdef NTL_ZZ_NBITS
         	#undef NTL_ZZ_NBITS
         	#endif
         	#define NTL_ZZ_NBITS (1200)
      	} else if (m <= 125) {
         	RR::SetPrecision(4000);
 
         	#ifdef NTL_ZZ_NBITS
         	#undef NTL_ZZ_NBITS
         	#endif
         	#define NTL_ZZ_NBITS (1250)
      	} else if (m <= 130) {
         	RR::SetPrecision(4500);
 
         	#ifdef NTL_ZZ_NBITS
         	#undef NTL_ZZ_NBITS
         	#endif
         	#define NTL_ZZ_NBITS (1300)
      	} else if (m <= 140) {
         	RR::SetPrecision(5000);
 
         	#ifdef NTL_ZZ_NBITS
         	#undef NTL_ZZ_NBITS
         	#endif
         	#define NTL_ZZ_NBITS (1400)
      	} else if (m <= 150) {
         	RR::SetPrecision(5500);
 
         	#ifdef NTL_ZZ_NBITS
         	#undef NTL_ZZ_NBITS
         	#endif
         	#define NTL_ZZ_NBITS (1500)
      	} else {
         	printf ("\nDimension not yet tested to find precisions. \n");
         	exit (3);
      	}
	return 0;
}

// Deallocate ppBasis
int deleteBasis (int ** ppBasis, int n) {
   for (int i=0; i<n; i++) {
      free (ppBasis[i]);
   }
   free (ppBasis);
   return 0;
}

// Deallocate ppBasisDouble
template <typename T1>
int deleteBasisDouble (T1 ** ppBasisDouble, int n) {
   for (int i=0; i<n; i++) {
      free (ppBasisDouble[i]);
   }
   free (ppBasisDouble);
   return 0;
}

// Compute SVP Challenge factor
RR SVP_Challenge_Factor(int n, ZZ volume) {
   	RR n_RR;
   	conv(n_RR, n);

   	RR power, volume_RR;
   	div(power, (RR) 1.0, n_RR);
   	conv(volume_RR, volume);
  
   	double gamma;
   	RR gamma_RR, nth_root_vol, nth_root_gamma; 
   	gamma = tgamma(n/2 + 1);
   	conv(gamma_RR, gamma);
  
   	nth_root_gamma = pow(gamma_RR, power); 
  
  	RR pi_RR, root_pi_RR;
   	ComputePi(pi_RR);
   	root_pi_RR = SqrRoot(pi_RR);
   
   	pow(nth_root_vol, volume_RR, power);

   	RR RRfactor;
   	RRfactor = (nth_root_gamma / root_pi_RR) * nth_root_vol;

   	return RRfactor;
}

// Comparator function to sort pairs according to second value
bool cmp_double(pair<string, double>& a, pair<string, double>& b){
    return a.second < b.second;
}

// Comparator function to sort pairs according to second value
bool cmp_long_double(pair<string, long double>& a, pair<string, long double>& b){
    return a.second < b.second;
}

// Function to sort the map according
// to value in a (key-value) pairs (fix)
void map_sort_long_double_fix(map<string, long double>& M){
  
    // Declare vector of pairs
    vector<pair<string, long double> > A;
  
    // Copy key-value pair from Map to vector of pairs
    for (auto& it : M) {
        A.push_back(it);
    }
  
    // Sort using comparator function
    sort(A.begin(), A.end(), cmp_long_double);
     
    long double min_elem = A[0].second;
    string min_label = A[0].first;
    
    cout << std::fixed;
    // Print the sorted value
    for (auto& it : A) {
       if (it.first == min_label) {
          cout << std::setprecision(15) << it.first << "\t" << it.second << "    =    X " << endl;
       } else {
	  long double mult = (it.second / min_elem);
          cout << std::setprecision(15) << it.first << "\t" << it.second << "    =    " << mult << "X " << endl;
       }
    }
}

// Function to sort the map according
// to value in a (key-value) pairs (fix)
void map_sort_double_fix(map<string, double>& M){
  
    // Declare vector of pairs
    vector<pair<string, double> > A;
  
    // Copy key-value pair from Map to vector of pairs
    for (auto& it : M) {
        A.push_back(it);
    }
  
    // Sort using comparator function
    sort(A.begin(), A.end(), cmp_double);
     
    double min_elem = A[0].second;
    string min_label = A[0].first;
    
    cout << std::fixed;
    // Print the sorted value
    for (auto& it : A) {
       if (it.first == min_label) {
          cout << std::setprecision(15) << it.first << "\t" << it.second << "    =    X " << endl;
       } else {
	  double mult = (it.second/min_elem);
          cout << std::setprecision(15) << it.first << "\t" << it.second << "    =    " << mult << "X " << endl;
       }
    }
}


// Compare elements for use in qsort algorithm
int compare (const void * a, const void * b) {
   if ( *(long double *) a < *(long double *) b ) return -1;
   if ( *(long double *) a == *(long double *) b ) return 0;
   if ( *(long double *) a > *(long double *) b ) return 1;
   return -2;
}

// Generate a random NTL basis (mat_ZZ) with n integer vectors, each of dimension m
mat_ZZ generateNTLChallengeBasis (int n) {
   mat_ZZ NTLBasis;
   ZZ prime;
   long length, err; 
   length = 10*n;
   GenPrime(prime, length, err = 80); 
   NTLBasis.SetDims(n,n);
   NTLBasis[0][0] = prime;
   for (int i=1; i<n; i++) {
      ZZ r;
      RandomBnd(r, prime);
      NTLBasis[i][0] = r;
      for (int j=1; j<n; j++){
         if (j == i) {
	    NTLBasis[i][j] = 1;
	 } else {
	    NTLBasis[i][j] = 0;
	 }	 
      }
   }
   return NTLBasis;   
}

// Find the inner products of vectors pV1 with pV2
template <typename T1, typename T2>
long double inner_product_template (T1 * pV1, T2 * pV2, int m) {
   long double result = 0;
   for (int i=0; i<m; i++) {
      result += (long double)pV1[i] * (long double)pV2[i];
   }
   return result;
}

// Compute norm of a vec_ZZ
RR NTL_vector_norm (vec_ZZ V1, int m) {
   ZZ norm_squared;
   InnerProduct(norm_squared, V1, V1);
   RR RR_norm_squared, norm;
   conv(RR_norm_squared, norm_squared);
   SqrRoot(norm, RR_norm_squared);

   return norm;
}

// Reduce vector pV2 with pV1
int reduce (int * pV1, int * pV2, int m, int closest_integer_to_mu) {
   for (int i=0; i<m; i++) {
      pV2[i] -= closest_integer_to_mu * pV1[i];
   }
   return 0;
}

template <typename T1>
int reduceDouble (T1 * pVDouble1, T1 * pVDouble2, int m, T1 mu) {
   for (int i=0; i<m; i++) {
      pVDouble2[i] -= mu * pVDouble1[i];
   }
   return 0;
}

// Reduce vec_ZZ V2 with vec_ZZ V1 (only NTL)
int NTL_reduce (vec_ZZ& V1, vec_ZZ& V2, int m, ZZ closest_ZZ_to_mu) {
   for (int i=0; i<m; i++) {
      ZZ V1_mu;
      mul(V1_mu, closest_ZZ_to_mu, V1[i]);
      sub(V2[i], V2[i], V1_mu);
   }
   return 0;
}

// Copy int **ppBasis to mat_ZZ BasisNTL and return
mat_ZZ copyBasisToNTL (int ** ppBasis, int m, int n) {
   mat_ZZ BasisNTL;

   BasisNTL.SetDims(n,m);
   for (int i=0; i<m; i++) {
      for (int j=0; j<n; j++) {
         BasisNTL[i][j] = ppBasis[i][j];
      }
   }
   return BasisNTL;
}


// Copy mat_ZZ ppBasis to mat_RR BasisNTL and return
mat_RR copyNTLBasisToRR (mat_ZZ ppBasis, int m, int n) {
   mat_RR BasisRR;

   BasisRR.SetDims(n,m);
   for (int i=0; i<m; i++) {
      for (int j=0; j<n; j++) {
         conv(BasisRR[i][j], ppBasis[i][j]);
      }
   }
   return BasisRR;
}

// Copy pBasis into ppBasisDouble and return
long double ** copyBasisToDouble (int ** ppBasis, int m, int n) {
   long double ** ppBasisDouble;
   // Create the 2-dimensional array
   ppBasisDouble = (long double **) calloc (n, sizeof(long double *));
   for (int i=0; i<n; i++) {
      ppBasisDouble[i] = (long double *) calloc (m, sizeof(long double));
      // Copy the values
      for (int j=0; j<m; j++) {
         ppBasisDouble[i][j] = (long double) ppBasis[i][j];
      }
   }
   return ppBasisDouble;
}

//  Rounds a/d to nearest integer, breaking ties
//  by rounding towards zero.  Assumes d > 0.
//  Stored in q
//  Copied from NTL Library (github) 
int BalDiv(ZZ& q, const ZZ& a, const ZZ& d) {

   
	NTL_ZZRegister(r);
	DivRem(q, r, a, d);
	add(r, r, r);
	
	long cmp = compare(r, d);
	if (cmp > 0 || (cmp == 0 && q < 0))
      	add(q, q, 1);
   	return 0;
}

// Exact division: q = a/b where b|a
// Copied from NTL Library (github) 
int ExactDiv(ZZ& qq, const ZZ& a, const ZZ& b) {
	NTL_ZZRegister(q);
	NTL_ZZRegister(r);

	DivRem(q, r, a, b);
	if (!IsZero(r)) {
		cerr << "a = " << a << "\n";
		cerr << "b = " << b << "\n";
	   	LogicError("ExactDiv: nonzero remainder");
	}
   	qq = q;
   	return 0;
}

// Gram-Schmidt Orthogonalisation
long double ** GSO (int ** pBasis, int m, int n) {
   long double ** ppBasisGSO;
   ppBasisGSO = copyBasisToDouble (pBasis, m, n);

   for (int i=1; i<n; i++) {
      for (int j=i-1; j>=0; j--) {

         // Find the value of mu_ij
         long double mu_ij =
                 inner_product_template (pBasis[i], ppBasisGSO[j], m) /
                 inner_product_template (ppBasisGSO[j], ppBasisGSO[j], m);

         // Reduce vector ppBasisGSO[i] with ppBasisGSO[j]
         // Compute b_i - mu_ij x b_j*
         reduceDouble (ppBasisGSO[j], ppBasisGSO[i], m, mu_ij);
      }
   }
   return ppBasisGSO;
}

// Gram-Schmidt Orthogonalisation (NTL datatypes)
mat_RR GSO_RR (mat_ZZ ppBasis, int m, int n) {
   mat_RR ppBasisGSO;
   mat_RR ppBasisRR;
   ppBasisGSO.SetDims(m, n);
   ppBasisRR.SetDims(m, n);
   for (int i=0; i<m; i++) {
      for (int j=0; j<n; j++) {
         conv(ppBasisGSO[i][j], ppBasis[i][j]);
         conv(ppBasisRR[i][j], ppBasis[i][j]);
      }
   }
   
  for (int i=0; i<n; i++) {
     for (int j=i-1; j>=0; j--) {
        RR mu_ij, num, den;
	InnerProduct(num, ppBasisRR[i], ppBasisGSO[j]);
	InnerProduct(den, ppBasisGSO[j], ppBasisGSO[j]);

	div(mu_ij, num, den);

	vec_RR product;
	mul(product, mu_ij, ppBasisGSO[j]);
	sub(ppBasisGSO[i], ppBasisGSO[i], product);
     }
  }
  return ppBasisGSO; 
}

// Update the GSO information in the way required for DeepLLL
// Algorithm 4 in Yamaguchi Yasuda NuTMic paper (NuTMic 2017)
template <typename T1>
int deep_LLL_GSO_update (T1 ** ppM, T1 * pB, int k , int i, int n) {
   // Initialise arrays pP, pS, pD as described in [YY18]
   T1 * pP;
   pP = (T1 *) calloc (n, sizeof(T1));
   T1 * pS;
   pS = (T1 *) calloc (n, sizeof(T1));
   T1 * pD;
   pD = (T1 *) calloc (n, sizeof(T1));

   for (int j=0; j<n; j++) {
      pP[j] = 0.0;
      pS[j] = 0.0;
      pD[j] = 0.0;
   }


   pP[k] = pB[k];
   pD[k] = pB[k];

   // Lines 2-4 of Algorithm 4 in YY18
   for (int j=k-1; j>=i; j--) {
      pP[j] = ppM[k][j] * pB[j];
      pD[j] = pD[j+1] + (ppM[k][j]*pP[j]);
   }

   // Line 5 of Algorithm 4 in YY18 - Intialising values of S
   for (int j=i; j<n; j++) {
      pS[j] = (T1) 0.0;
   }

   // Lines 6 - 18 of Algorithm 4 in YY18
   for (int j=k; j>i; j--) {
      T1 T = 0.0;
      T = ppM[k][j-1] / pD[j];
      for (int l=n-1; l>k;l--) {
         pS[l] = pS[l] + (ppM[l][j] * pP[j]);
         ppM[l][j] = ppM[l][j-1] - (T * pS[l]);
      }

      for (int l=k; l>=j+2;l--) {
         pS[l] = pS[l] + (ppM[l-1][j] * pP[j]);
         ppM[l][j] = ppM[l-1][j-1] - (T * pS[l]);
      }

      if (j != k) {
         pS[j+1] = pP[j];
         ppM[j+1][j] = ppM[j][j-1] - (T * pS[j+1]);
      }
   }

   // Lines 20 - 22 of Algorithm 4 in YY18
   T1 T = 0.0;
   T = (T1) 1.0 / pD[i];


   for (int l=n-1; l>k; l--) {
      ppM[l][i] = T * (pS[l] + ppM[l][i] * pP[i]);
   }
   for (int l=k; l>=i+2; l--) {
      ppM[l][i] = T * (pS[l] + ppM[l-1][i] * pP[i]);
   }

   // Lines 23 - 29 of Algorithm 4 in YY18
   ppM[i+1][i] = T * pP[i];

   for (int j=0; j<i; j++) {
      T1 epsilon = 0.0;
      epsilon = ppM[k][j];
      for (int l=k; l>i; l--) {
         ppM[l][j] = ppM[l-1][j];
      }
         ppM[i][j] = epsilon;
   }

   // Lines 31 - 32 of Algorithm 4 in YY18
   for (int j=k; j>i;j--) {
      pB[j] = (pD[j] * pB[j-1]) / pD[j-1];
   }
   pB[i] = pD[i];

   return 0;
}

// Update the GSO information in the way required for DeepLLL (RR datatypes)
// Algorithm 4 in Yamaguchi Yasuda NuTMic paper (NuTMic 2017)
int deep_LLL_GSO_update_RR (mat_RR &ppM, vec_RR &pB, int k , int i, int n) {
   // Initialise arrays pP, pS, pD as described in [YY18]
   vec_RR pP;
   pP.SetLength(n);
   vec_RR pS;
   pS.SetLength(n);
   vec_RR pD;
   pD.SetLength(n);

   for (int j=0; j<n; j++) {
      pP[j] = 0;
      pS[j] = 0;
      pD[j] = 0;
   }


   pP[k] = pB[k];
   pD[k] = pB[k];

   // Lines 2-4 of Algorithm 4 in YY18
   for (int j=k-1; j>=i; j--) {
      pP[j] = ppM[k][j] * pB[j];
      pD[j] = pD[j+1] + (ppM[k][j]*pP[j]);
   }

   // Line 5 of Algorithm 4 in YY18 - Intialising values of S
   for (int j=i; j<n; j++) {
      pS[j] = 0;
   }

   // Lines 6 - 18 of Algorithm 4 in YY18
   for (int j=k; j>i; j--) {
      RR T;
      T = 0;
      T = ppM[k][j-1] / pD[j];
      for (int l=n-1; l>k;l--) {
         pS[l] = pS[l] + (ppM[l][j] * pP[j]);
         ppM[l][j] = ppM[l][j-1] - (T * pS[l]);
      }

      for (int l=k; l>=j+2;l--) {
         pS[l] = pS[l] + (ppM[l-1][j] * pP[j]);
         ppM[l][j] = ppM[l-1][j-1] - (T * pS[l]);
      }

      if (j != k) {
         pS[j+1] = pP[j];
         ppM[j+1][j] = ppM[j][j-1] - (T * pS[j+1]);
      }
   }

   // Lines 20 - 22 of Algorithm 4 in YY18
   RR T;
   inv(T, pD[i]);


   for (int l=n-1; l>k; l--) {
      ppM[l][i] = T * (pS[l] + ppM[l][i] * pP[i]);
   }
   for (int l=k; l>=i+2; l--) {
      ppM[l][i] = T * (pS[l] + ppM[l-1][i] * pP[i]);
   }

   // Lines 23 - 29 of Algorithm 4 in YY18
   ppM[i+1][i] = T * pP[i];

   for (int j=0; j<i; j++) {
      RR epsilon;
      epsilon = ppM[k][j];
      for (int l=k; l>i; l--) {
         ppM[l][j] = ppM[l-1][j];
      }
         ppM[i][j] = epsilon;
   }

   // Lines 31 - 32 of Algorithm 4 in YY18
   for (int j=k; j>i;j--) {
      pB[j] = (pD[j] * pB[j-1]) / pD[j-1];
   }
   pB[i] = pD[i];

   return 0;
}

// Find the indices k,i where we perform a deep insertion for SS-GGLLL (standard C datatypes)
template <typename T1>
T1 find_insertion_indices_SSGGLLL (T1 ** ppM, int ** ppBasis, T1 * pB, int n, int m, int * pK, int * pJ) {

   int k_max = 1;
   int i_max = 0;

   T1 S_max = 0.0;
   T1 S_ik = 0.0;

   // For k = 1,...,n and j = k-1,...,1,0, we find the indices k, i where we insert b_k to reduce the value of SS the most
   for (int k=1; k<n; k++) {
      // Initialise to the k-1 case
      T1 delta_k_1 = 0.0;
      delta_k_1 = (pB[k] + (ppM[k][k-1] * ppM[k][k-1] * pB[k-1])) / pB[k-1];

      // Initialise D_j = D_{k-1}
      T1 D_j = 0.0;
      D_j = (ppM[k][k-1] * ppM[k][k-1] * pB[k-1]) + pB[k];

      // Initialise S_ik = S_{i-1,k}
      // After initialising S and D, we can update them in the 'for' loop below for general S_jk and D_j for 0 <= j < k-1
      T1 S_ik = 0.0;
      S_ik = ppM[k][k-1] * ppM[k][k-1] * (1 - delta_k_1);

      // Check if k, k-1 case improves upon current best
      if (S_ik > S_max) {
         S_max = S_ik;
         k_max = k;
         i_max = k-1;
      }

      for (int l=k-2; l>=0; l--) {
         D_j += ppM[k][l] * ppM[k][l] * pB[l];
         S_ik += ppM[k][l] * ppM[k][l] * pB[l] * ((pB[l]/D_j) - 1);

         if (S_ik > S_max) {
            k_max = k;
            i_max = l;
            S_max = S_ik;
         }
      }
   }
   *pK = k_max;
   *pJ = i_max;

   return S_max;
}


RR RR_find_insertion_indices_SSGGLLL (mat_RR ppM, mat_ZZ ppBasis, vec_RR pB, int n, int m, int * pK, int * pJ) {

   int k_max = 1;
   int i_max = 0;

   RR S_max, S_ik; 
   S_max = 0;
   S_ik = 0;

   // For k = 1,...,n and j = k-1,...,1,0, we find the indices k, i where we insert b_k to reduce the value of SS the most
   for (int k=1; k<n; k++) {
      // Initialise to the k-1 case
      RR delta_k_1;
      delta_k_1 = 0;
      delta_k_1 = (pB[k] + (ppM[k][k-1] * ppM[k][k-1] * pB[k-1])) / pB[k-1];

      // Initialise D_j = D_{k-1}
      RR D_j;
      D_j = 0;
      D_j = (ppM[k][k-1] * ppM[k][k-1] * pB[k-1]) + pB[k];

      // Initialise S_ik = S_{i-1,k}
      // After initialising S and D, we can update them in the 'for' loop below for general S_jk and D_j for 0 <= j < k-1
      RR S_ik;
      S_ik = 0;
      S_ik = ppM[k][k-1] * ppM[k][k-1] * (1 - delta_k_1);

      // Check if k, k-1 case improves upon current best
      if (S_ik > S_max) {
         S_max = S_ik;
         k_max = k;
         i_max = k-1;
      }

      for (int l=k-2; l>=0; l--) {
         D_j += ppM[k][l] * ppM[k][l] * pB[l];
         S_ik += ppM[k][l] * ppM[k][l] * pB[l] * ((pB[l]/D_j) - 1);

         if (S_ik > S_max) {
            k_max = k;
            i_max = l;
            S_max = S_ik;
         }
      }
   }
   *pK = k_max;
   *pJ = i_max;

   return S_max;
}

// Find insertion indices for Pot-GGLLL (standard C datatypes) 
long double find_insertion_indices_PotGGLLL (long double ** ppM, int ** ppBasis, long double * pB, int n, int m, int * pK, int * pJ) {

   int k_min = 0;
   int i_min = 0;

   // Initialise P = 1.0, P_min = 1.0
   // We reinitialise P to be 1.0 every time we increment k (in standard PotLLL, for each k, we start with P = 1)
   long double P = 1.0;
   long double P_min = 1.0;

   // Initialise P_min to be the k=1, i=0 case
   long double P_init = 1.0;
   P_init = (pB[1] + ppM[1][0]*ppM[1][0]*pB[0]) / pB[0];

   if (P_init < P_min) {
      P_min = P_init;
      k_min = 1;
      i_min = 0;
   }
   long double sum = 0.0;
   long double P_mult = 0.0;
   for (int k=2;k<n;k++) {
      P = 1;
      for (int j=k-1; j>=0; j--) {
	 sum = 0.0;
         for (int r=j; r<k; r++) {
            sum += ppM[k][r]*ppM[k][r]*pB[r];
         }

         P_mult = (pB[k] + sum) / pB[j];
         P = P * P_mult;
         if (P<P_min) {
            k_min = k;
            i_min = j;
            P_min = P;
         }
      }
   }

   *pK = k_min;
   *pJ = i_min;

   return P_min;
}

// Find insertion indices for Pot-GGLLL (using NTL datatypes)
RR RR_find_insertion_indices_PotGGLLL (mat_RR ppM, mat_RR ppDelta, mat_ZZ ppBasis, vec_RR pB, int n, int m, int * pK, int * pJ) {

   int k_min = 0;
   int i_min = 0;

   // Initialise P = 1.0, P_min = 1.0
   // We reinitialise P to be 1.0 every time we increment k (in standard PotLLL, for each k, we start with P = 1)
   RR P, P_min;
   P = 1;
   P_min = 1;

   // Initialise P_min to be the k=1, i=0 case
   RR P_init;
   P_init = (pB[1] + ppM[1][0]*ppM[1][0]*pB[0]) / pB[0];

   if (P_init < P_min) {
      P_min = P_init;
      k_min = 1;
      i_min = 0;
   }

   for (int k=2;k<n;k++) {
      P = 1;
      for (int j=k-1; j>=0; j--) {
         RR sum;
	 sum = 0;
         for (int r=j; r<k; r++) {
            sum += ppM[k][r]*ppM[k][r]*pB[r];
         }

      RR P_mult;
      P_mult = (pB[k] + sum) / pB[j];
      P = P * P_mult;
      if (P<P_min) {
         k_min = k;
         i_min = j;
         P_min = P;
      }
      }
   }
   *pK = k_min;
   *pJ = i_min;

   return P_min;
}



// RR Root Hermite Factor
RR NTLrootHermiteFactor (mat_ZZ ppBasis, int m, int n) {
   RR RR_n;
   conv(RR_n, n);
   vec_RR V1;
   V1.SetLength(m);
   mat_RR ppBasisGSO;
   ppBasisGSO = GSO_RR (ppBasis, m, n);
   
   for (int i=0; i<m; i++){
      conv(V1[i], ppBasis[0][i]);
   }
   RR volume;
   volume = 1;

   for (int i=0; i<n; i++) {
      RR gsoVectorNormSquared, gsoVectorNorm;
      InnerProduct(gsoVectorNormSquared, ppBasisGSO[i],ppBasisGSO[i]);
      SqrRoot(gsoVectorNorm, gsoVectorNormSquared);

      mul(volume, volume, gsoVectorNorm);
   }
   RR power;
   inv(power, RR_n);

   RR nthRootVolume;
   pow(nthRootVolume, volume, power);

   RR b1NormSquared;
   InnerProduct (b1NormSquared, V1, V1);

   RR b1Norm;
   SqrRoot(b1Norm, b1NormSquared);

   RR hermiteFactor, RHF;
   div(hermiteFactor, b1Norm, nthRootVolume);

   pow(RHF, hermiteFactor, power);

   return RHF;   
}

RR NTLorthogonalityDefect (mat_ZZ ppBasis, int m, int n) {
   mat_RR ppBasisGSO;
   ppBasisGSO = GSO_RR(ppBasis, m, n);

   RR volumeGSO;
   volumeGSO = 1;
   for (int i=0; i<n; i++) {
      RR gsoVectorNormSquared;
      InnerProduct(gsoVectorNormSquared, ppBasisGSO[i], ppBasisGSO[i]);
      mul(volumeGSO, volumeGSO, gsoVectorNormSquared); 
   }
   SqrRoot(volumeGSO, volumeGSO);

   RR volumeBasis;
   ZZ volumeBasisSquared;
   volumeBasisSquared = 1;
   for (int j=0; j<n; j++) {
      ZZ vectorNormSquared;
      InnerProduct(vectorNormSquared, ppBasis[j], ppBasis[j]);
      mul(volumeBasisSquared, volumeBasisSquared, vectorNormSquared);
   }
   conv (volumeBasis, volumeBasisSquared);
   SqrRoot(volumeBasis,volumeBasis);

   RR orthDefect;
   div(orthDefect, volumeBasis, volumeGSO);
   return orthDefect;
}

RR NTLsquaredSum (mat_ZZ ppBasis, int m, int n) {
   mat_RR ppBasisGSO;
   ppBasisGSO = GSO_RR(ppBasis, m, n);
   RR squaredSum;
   squaredSum = 0;

   for (int i=0; i<n; i++) {
      RR gsoNormSquared;
      InnerProduct(gsoNormSquared, ppBasisGSO[i], ppBasisGSO[i]);
      add(squaredSum, squaredSum, gsoNormSquared);
   }

   return squaredSum;
}

RR NTLbasisPotential (mat_ZZ ppBasis, int m, int n) {
   mat_RR ppBasisGSO;
   ppBasisGSO = GSO_RR(ppBasis, m, n);

   RR potential;
   potential = 0;
   for (int i=0; i<n; i++) {
      // The potential term for b_i is ||b_i*||^(2(n-i)|| (since indexing from 0)
      // This is the same as (||b_i*||^2)^(n-i)
      // Log of the potential term is (n-i) * log(||b_i*||^2)
      RR gsoVectorNormSquared;
      InnerProduct(gsoVectorNormSquared, ppBasisGSO[i], ppBasisGSO[i]);
      RR RR_power;
      int power = n-i;
      conv(RR_power, power);
      RR logNormSquared;
      log(logNormSquared, gsoVectorNormSquared);
      RR potentialTerm; 
      mul(potentialTerm, RR_power, logNormSquared);
      potential += potentialTerm;
   }
   return potential;
}

// Check if a Basis is LLL reduced
int NTL_LLL_check (mat_ZZ ppBasis, int m, int n, RR delta_threshold) {

   int isLLLReduced = TRUE;

   mat_RR ppBasisGSO;

   mat_RR ppBasisRR;
   ppBasisRR = copyNTLBasisToRR(ppBasis, m, n);

   // Find the GSO of the Basis
   ppBasisGSO = GSO_RR (ppBasis, m, n);

   #ifdef __DEBUG__
   cout << "\n\n\nThe GSO of the basis is:\n" << ppBasisGSO;
   #endif

   // This 1-dimensional array pB will store all the values of ||b_i*||^2
   vec_RR pB;
   pB.SetLength(n);
   #ifdef __DEBUG__
   printf ("\nThe values of ||b_i*||^2 are:\n");
   #endif
   for (int i=0; i<n; i++) {
      InnerProduct(pB[i], ppBasisGSO[i], ppBasisGSO[i]);
      #ifdef __DEBUG__
      cout << pB[i] << endl;
      #endif
   }

   // This 2-dimensional nxn array ppM will store all the values of mu_ij
   mat_RR ppM;
   // This 2-dimensional nxn array ppDelta will store all the values of delta_ij
   mat_RR ppDelta;
   ppM.SetDims(n,m);
   ppDelta.SetDims(n,m);
   #ifdef __DEBUG__
   cout << "The values of mu_ij are:\n" << ppM << endl;
   #endif
   for (int i=1; i<n; i++) {
      for (int j=0; j<i; j++) {
	 RR inner;
         InnerProduct(inner, ppBasisRR[i], ppBasisGSO[j]);	 
         div(ppM[i][j], inner, pB[j]);
         if (ppM[i][j]<-0.5 || ppM[i][j]>0.5) {
            isLLLReduced = FALSE;
            cout << "\nmu" << i << j << "  " << ppM[i][j] << endl;
            return isLLLReduced;
         }
	 RR frac, product;
	 div(frac, pB[i], pB[j]);
         mul(product, ppM[i][j], ppM[i][j]);
	 add(ppDelta[i][j], frac, product);
         #ifdef __DEBUG__
         cout << "\nMU_ij  " << ppM[i][j] << "\nDelta_ij  " << ppDelta[i][j] << endl;
         #endif
      }
      #ifdef __DEBUG__
      printf ("\n");
      #endif
      // Check the Lovaz Condition
      if (pB[i] < ((delta_threshold - (ppM[i][i-1] * ppM[i][i-1])) * pB[i-1])) {
         isLLLReduced = FALSE;
         printf ("\nLC[%d][%d]", i, i-1);
         return isLLLReduced;
      }
   }
    
   ppM.kill();
   ppDelta.kill();
   ppBasisGSO.kill();
   ppBasisRR.kill();

   return isLLLReduced;
}

// The standard LLL implementation using NTL mat_ZZ datatype for basis
mat_ZZ LLL_mat_ZZ (mat_ZZ ppBasis, int m, int n, RR delta, long long int * pSwaps, long long int * pReductions) {

   mat_RR ppBasisGSO;
   mat_RR ppBasisRR;
   // Delta needs to be type RR
   ppBasisRR = copyNTLBasisToRR(ppBasis, m, n);
   // Find the GSO of the Basis
   ppBasisGSO = GSO_RR (ppBasis, m, n);

   #ifdef __DEBUG__
   cout <<"\n\n\nThe GSO of the basis is:\n" << ppBasisGSO << endl;
   #endif

   // This 1-dimensional array pB will store all the values of ||b_i*||^2
   vec_RR pB;
   pB.SetLength(n);
   #ifdef __DEBUG__
   printf ("\nThe values of ||b_i*||^2 are:\n");
   #endif
   for (int i=0; i<n; i++) {
      InnerProduct(pB[i], ppBasisGSO[i], ppBasisGSO[i]);
      #ifdef __DEBUG__
      cout << pB[i];
      #endif
   }

   // This 2-dimensional nxn array ppM will store all the values of mu_ij
   mat_RR ppM, ppDelta;
   ppM.SetDims(n,m);
   ppDelta.SetDims(n,m);
   #ifdef __DEBUG__
   printf ("\nThe values of mu_ij are:");
   #endif
   for (int i=1; i<n; i++) {
      for (int j=0; j<i; j++) {
	 RR inner;
	 InnerProduct(inner, ppBasisRR[i], ppBasisGSO[j]);
         div(ppM[i][j], inner, pB[j]);
         #ifdef __DEBUG__
         cout << ppM[i][j];
         #endif
      }
      #ifdef __DEBUG__
      printf ("\n");
      #endif
   }

   // The LLL task
   *pSwaps = 0;
   *pReductions = 0;
   int fullReductions;
   fullReductions = 0;
   int k = 1;
   while (k<n) {


      // Reduce b_k with all the previous basis vectors
      for (int j=k-1; j>=0; j--) {
         //double half_plus_fudge = 0.55;
	 //If want fudge, change the lines below from 0.50 to half_plus_fudge
         RR closest_int_RR;
	 ZZ closest_integer_to_mu_kj;
	 round(closest_int_RR, ppM[k][j]);
	 conv(closest_integer_to_mu_kj, closest_int_RR);
         
	 if (closest_integer_to_mu_kj != 0) {
            // Count the number of reductions
            (*pReductions)++;
	    // Reduce b_k with b_j
            NTL_reduce (ppBasis[j], ppBasis[k], m, closest_integer_to_mu_kj);
   
            // mu_ki -= CLOSEST_INTEGER (mu_kj)
            for (int i=j-1; i>=0; i--) {
               // By Exercise 17.4.4 from Galbraith v2i
	       RR term;
	       mul (term, closest_int_RR, ppM[j][i]);
               sub(ppM[k][i], ppM[k][i], term);
            }

            // mu_kj = mu_kj - [mu_kj]
            sub(ppM[k][j], ppM[k][j], closest_int_RR);
         }
      }
      (fullReductions++) ;

      // Check Lovasz condition
      // if (LC (pB[k], pB[k-1], delta, ppM[k][k-1]))
      if (pB[k] >= ((delta - (ppM[k][k-1] * ppM[k][k-1])) * pB[k-1])) {

         #ifdef __DEBUG__
         cout << "\nLC succeeded for k="<< k << " delta=" << ((pB[k]/pB[k-1]) + (ppM[k][k-1] * ppM[k][k-1]));
         cout << "\nLC succeeded for k="<< k << " delta=" << ((pB[k]/pB[k-1]) + (ppM[k][k-1] * ppM[k][k-1]));
         #endif
         k++;

      } else {

         #ifdef __DEBUG__
         printf ("\nLC failed for k=" << k << " delta=%Lf" << ((pB[k]/pB[k-1]) + (ppM[k][k-1] * ppM[k][k-1]));
         #endif

         // Swap vectors k, k-1
         //swap(ppBasis[k], ppBasis[k-1]);
	 vec_ZZ pTemp;
         pTemp = ppBasis[k];
         ppBasis[k] = ppBasis[k-1];
         ppBasis[k-1] = pTemp;

         // Count the number of swaps
         (*pSwaps)++;

         // Swap the first k-2 elements in rows k, k-1 of matrix M
         for (int c=k-2; c>=0; c--) {
            // Swap rows k-1 and k: ppM[k-1][c] with ppM[k][c]
            RR temp;
	    temp = ppM[k-1][c];
            ppM[k-1][c] = ppM[k][c];
            ppM[k][c] = temp;
         }

         // Updating B_k and B_{k-1}
         RR mu_k_k1;
	 mu_k_k1 = ppM[k][k-1];
         RR B_k1_updated;
	 B_k1_updated = pB[k] + (mu_k_k1*mu_k_k1)*pB[k-1];
         ppM[k][k-1] = (mu_k_k1 * pB[k-1]) / B_k1_updated;
         pB[k] = (pB[k-1] * pB[k]) / B_k1_updated;
         pB[k-1] = B_k1_updated;

         // Update columns k-1 and k in M
         for (int r=k+1; r<n; r++) {
            RR T;
	    T = ppM[r][k];
            ppM[r][k] = ppM[r][k-1] - (mu_k_k1 * T);
            ppM[r][k-1] = T + (ppM[k][k-1] * ppM[r][k]);
         }

         #ifdef __DEBUG__
         printf ("\nThe Basis vectors after the swap are:\n");
         cout << ppBasis;
         printf ("\nThe updated values of mu_ij are:");
         for (int i=0; i<n; i++) {
            for (int j=0; j<i; j++) {
               cout << ppM[i][j];
            }
            printf ("\n");
         }
         #endif

         // Update k
         k = (1 > (k-1)) ? 1 : k-1;
      }

   }

   ppM.kill();
   ppBasisGSO.kill();
   ppBasisRR.kill();

   return ppBasis;
}

// The PotLLL implementation with NTL datatypes
// At every iteration, find insert the vector b_k position i such that the basis potential reduces by a sufficiently large amount
mat_ZZ NTL_PotLLL (mat_ZZ ppBasis, int m, int n, RR delta, long long int * pSwaps, long long int * pReductions) {

   mat_RR ppBasisGSO;
   ppBasisGSO.SetDims(n,m);
   // Find the GSO of the Basis
   ppBasisGSO = GSO_RR (ppBasis, m, n);

   #ifdef __DEBUG__
   printf ("\n\n\nThe GSO of the basis is:\n");
   cout << ppBasisGSO << endl;
   #endif

   // This 1-dimensional array pB will store all the values of ||b_i*||^2
   vec_RR pB;
   pB.SetLength(n);
   #ifdef __DEBUG__
   printf ("\nThe values of ||b_i*||^2 are:\n");
   #endif
   for (int i=0; i<n; i++) {
      InnerProduct(pB[i], ppBasisGSO[i], ppBasisGSO[i]);
      #ifdef __DEBUG__
      cout << pB[i] << endl;
      #endif
   }

   // This 2-dimensional nxn array ppM will store all the values of mu_ij
   mat_RR ppM;
   // This 2-dimensional nxn array ppDelta will store all the values of delta_ij
   mat_RR ppDelta;
   ppM.SetDims(n,m);
   ppDelta.SetDims(n,m);
   #ifdef __DEBUG__
   printf ("\n\nThe values of mu_ij are:");
   #endif

   mat_RR ppBasisRR;
   ppBasisRR = copyNTLBasisToRR(ppBasis, m, n);

   for (int i=1; i<n; i++) {
      for (int j=0; j<i; j++) {
	 RR inner;
	 InnerProduct(inner, ppBasisRR[i], ppBasisGSO[j]);
	 div(ppM[i][j], inner, pB[j]);
         ppDelta[i][j] = ((pB[i]/pB[j]) + (ppM[i][j] * ppM[i][j]));
         #ifdef __DEBUG__
         cout << ppM[i][j] << endl << ppDelta[i][j] << endl;
         #endif
      }
      #ifdef __DEBUG__
      printf ("\n");
      #endif
   }

   // The PotLLL task
   *pSwaps = 0;
   *pReductions = 0;
   int k = 1;

   while (k<n) {

      // Reduce b_k with all previous basis vectors
      for (int j=k-1; j>=0; j--) {
         RR closest_integer_RR;
	 round(closest_integer_RR, ppM[k][j]);
         ZZ closest_integer_to_mu_kj;
	 conv(closest_integer_to_mu_kj, closest_integer_RR);
         
	 if (closest_integer_to_mu_kj != 0) {
            // Count the number of reductions
            (*pReductions)++;
            // Reduce b_k with b_j
            NTL_reduce (ppBasis[j], ppBasis[k], m, closest_integer_to_mu_kj);

            // mu_ki -= CLOSEST_INTEGER (mu_kj);
            for (int i = j-1; i>=0; i--) {
               // By Exercise 17.4.4 from Galbraith v2
               ppM[k][i] -= closest_integer_RR * ppM[j][i];
            }

            // mu_kj = mu_kj - [mu_kj]
            ppM[k][j] -= closest_integer_RR;
         }
      }
      // Initialise P = 1, P_min = 1, k = 0
      // For j = k-1,...,1,0, we find the index where an insertion of vector b_k in position j reduces the potential the most
      RR P, P_min;
      P = 1;
      P_min = 1;
      int l = 1;

      for (int j=k-1; j>=0; j--) {
         RR sum;
	 sum = 0;
         for (int r=j; r<k; r++) {
            sum += ppM[k][r]*ppM[k][r]*pB[r];
         }

         RR P_mult;
         P_mult = (pB[k] + sum) / pB[j];
         P = P * P_mult;
         if (P<P_min) {
            l = j;
            P_min = P;
         }
      }

      #ifdef __DEBUG__
      printf ("\nP_min = %Lf", P_min);
      #endif

      if (delta > P_min) {
         // b_l <- b_k, and vectors b_l, ... , b_k-1 are shifted up one position
         vec_ZZ pTemp;
         pTemp = ppBasis[k];

         for (int j=k-1; j>=l; j--) {
            ppBasis[j+1] = ppBasis[j];
         }
         ppBasis[l] = pTemp;

         // Update values of ppM, pB
         (*pSwaps)++;

         deep_LLL_GSO_update_RR (ppM, pB, k , l, n);
         k = (1 > (l)) ? 1 : l;
         } else {
         k++;
      }
   }
   ppM.kill();
   ppDelta.kill();
   ppBasisGSO.kill();
   
   return ppBasis;
}

// The SS-LLL implementation using NTL Datatypes
// At every iteration, find insert the vector b_k position i such that the value of S^2 reduces as much as possible
mat_ZZ NTL_SS_LLL (mat_ZZ ppBasis, int m, int n, RR eta, long long int * pSwaps, long long int * pReductions) {

   mat_RR ppBasisGSO, ppBasisRR;
   ppBasisGSO.SetDims(n,m);
   ppBasisRR.SetDims(n,m);

   ppBasisRR = copyNTLBasisToRR(ppBasis, m, n);
   // Find the GSO of the Basis
   ppBasisGSO = GSO_RR (ppBasis, m, n);

   #ifdef __DEBUG__
   printf ("\n\n\nThe GSO of the basis is:\n");
   cout << ppBasisGSO;
   #endif

   // This 1-dimensional array pB will store all the values of ||b_i*||^2
   vec_RR pB;
   pB.SetLength(n);
   #ifdef __DEBUG__
   printf ("\nThe values of ||b_i*||^2 are:\n");
   #endif
   for (int i=0; i<n; i++) {
      InnerProduct(pB[i], ppBasisGSO[i], ppBasisGSO[i]);
      #ifdef __DEBUG__
      cout << pB[i];
      #endif
   }

   // This 2-dimensional nxn array ppM will store all the values of mu_ij
   mat_RR ppM;
   // This 2-dimensional nxn array ppDelta will store all the values of delta_ij
   mat_RR ppDelta;
   ppM.SetDims(n,m);
   ppDelta.SetDims(n,m);
   #ifdef __DEBUG__
   printf ("\n\nThe values of mu_ij are:");
   #endif
   for (int i=1; i<n; i++) {
      for (int j=0; j<i; j++) {
	 RR inner;
	 InnerProduct(inner, ppBasisRR[i], ppBasisGSO[j]);
         div(ppM[i][j], inner, pB[j]);
         ppDelta[i][j] = ((pB[i]/pB[j]) + (ppM[i][j] * ppM[i][j]));
         #ifdef __DEBUG__
         cout << ppM[i][j] << endl << ppDelta[i][j] << endl;
         #endif
      }
      #ifdef __DEBUG__
      printf ("\n");
      #endif
   }

   // The S^2LLL task
   *pSwaps = 0;
   *pReductions = 0;
   int k = 1;
   RR eta_p;
   eta_p = 0.000001;

   while (k<n) {

      RR squared_sum;
      squared_sum = 0;
      for (int i=0; i<n; i++) {
         squared_sum += pB[i];
      }

      // Reduce b_k with all previous basis vectors
      for (int j=k-1; j>=0; j--) {
         RR closest_integer_RR;
	 round(closest_integer_RR, ppM[k][j]);  
         ZZ closest_integer_to_mu_kj;
	 conv(closest_integer_to_mu_kj, closest_integer_RR);
         
	 if (closest_integer_to_mu_kj != 0) {
            // Count the number of reductions
            (*pReductions)++;
            // Reduce b_k with b_j
            NTL_reduce (ppBasis[j], ppBasis[k], m, closest_integer_to_mu_kj);

            // mu_ki -= CLOSEST_INTEGER (mu_kj);
            for (int i = j-1; i>=0; i--) {
               // By Exercise 17.4.4 from Galbraith v2
               ppM[k][i] -= closest_integer_RR * ppM[j][i];
            }

            // mu_kj = mu_kj - [mu_kj]
            sub(ppM[k][j], ppM[k][j], closest_integer_RR);
         }
      }

      // For j = k-1,...,1,0, we find the index i where we insert b_k to reduce the value of SS the most
      RR delta_k_1;
      delta_k_1 = 0.0;
      delta_k_1 = (pB[k] + (ppM[k][k-1] * ppM[k][k-1] * pB[k-1])) / pB[k-1];

      // Initialise D_j = D_{k-1}
      RR D_j;
      D_j = (ppM[k][k-1] * ppM[k][k-1] * pB[k-1]) + pB[k];

      // Initialise S_ik = S_{i-1,k}
      // After initialising S and D, we can update them in the 'for' loop below for general S_jk and D_j for 0 <= j < k-1
      RR S_ik;
      S_ik = ppM[k][k-1] * ppM[k][k-1] * (1 - delta_k_1);

      int i = k-1;

      RR S_lk;
      S_lk = S_ik;

      for (int l=k-2; l>=0; l--) {
         D_j += ppM[k][l] * ppM[k][l] * pB[l];
         S_lk += ppM[k][l] * ppM[k][l] * pB[l] * ((pB[l]/D_j) - 1);

         if (S_lk > S_ik) {
            i = l;
            S_ik = S_lk;
         }
      }

      if (S_ik <= eta_p * squared_sum) {
         k++;
      } else {
         // b_l <- b_k, and vectors b_l, ... , b_k-1 are shifted up one position
         vec_ZZ pTemp;
         pTemp = ppBasis[k];

         for (int j=k-1; j>=i; j--) {
            ppBasis[j+1] = ppBasis[j];
         }
         ppBasis[i] = pTemp;

         // Update values of ppM, pB
         (*pSwaps)++;

         deep_LLL_GSO_update_RR (ppM, pB, k , i, n);
         k = (1 >(i)) ? 1 : i;
      }
   }
   ppM.kill();
   ppDelta.kill();
   ppBasisGSO.kill();
   ppBasisRR.kill();
   
   return ppBasis;
}

// The Pot-GGLLL implementation with NTL datatypes
// At every iteration, find indices k, i such that the Potential of the basis is reduced by the largest amount.
// When this happens we insert b_k in position i
mat_ZZ NTL_PotGGLLL (mat_ZZ ppBasis, int m, int n, RR delta_threshold, long long int * pSwaps, long long int * pReductions) {

   mat_RR ppBasisGSO, ppBasisRR;
   ppBasisGSO.SetDims(n,m);
   ppBasisRR.SetDims(n,m);

   // Find the GSO of the Basis
   ppBasisGSO = GSO_RR (ppBasis, m, n);
   ppBasisRR = copyNTLBasisToRR (ppBasis, m, n);
   #ifdef __DEBUG__
   printf ("\n\n\nThe GSO of the basis is:\n");
   cout << ppBasisGSO;
   #endif

   // This 1-dimensional array pB will store all the values of ||b_i*||^2
   vec_RR pB;
   pB.SetLength(n);
   
   #ifdef __DEBUG__
   printf ("\nThe values of ||b_i*||^2 are:\n");
   #endif
   for (int i=0; i<n; i++) {
      InnerProduct(pB[i], ppBasisGSO[i], ppBasisGSO[i]);
      #ifdef __DEBUG__
      cout << pB[i] << endl;
      #endif
   }

   // This 2-dimensional nxn array ppM will store all the values of mu_ij
   mat_RR ppM;
   // This 2-dimensional nxn array ppDelta will store all the values of delta_ij
   mat_RR ppDelta;
   ppM.SetDims(n,m);
   ppDelta.SetDims(n,m);
   #ifdef __DEBUG__
   printf ("\n\nThe values of mu_ij are:");
   #endif

   for (int i=1; i<n; i++) {
      for (int j=0; j<i; j++) {
	 RR inner;
         InnerProduct(inner, ppBasisRR[i], ppBasisGSO[j]);
	 div(ppM[i][j], inner, pB[j]);
         ppDelta[i][j] = ((pB[i]/pB[j]) + (ppM[i][j] * ppM[i][j]));
         #ifdef __DEBUG__
         cout << ppM[i][j] << endl << ppDelta[i][j] << endl;
         #endif
      }
      #ifdef __DEBUG__
      printf ("\n");
      #endif
   }

   // The greedy global LLL task
   *pSwaps = 0;
   *pReductions = 0;
   RR P_min;

   // First reduce every vector by all its previous vectors
   for (int k=1; k<n; k++) {

      // Reduce b_k with all the previous basis vectors
      for (int j=k-1; j>=0; j--) {
         RR closest_integer_RR;
	 round(closest_integer_RR, ppM[k][j]);
	 ZZ closest_integer_to_mu_kj;
	 conv(closest_integer_to_mu_kj, closest_integer_RR);

	 if (closest_integer_to_mu_kj != 0) {
            // Count the number of reductions
            (*pReductions)++;
            // Reduce b_k with b_j
            NTL_reduce (ppBasis[j], ppBasis[k], m, closest_integer_to_mu_kj);

            // mu_ki -= CLOSEST_INTEGER (mu_kj)
            for (int i=j-1; i>=0; i--) {
               // By Exercise 17.4.4 from Galbraith v2
               ppM[k][i] -= closest_integer_RR * ppM[j][i];
               ppDelta[k][i] = ((pB[k]/pB[i]) + (ppM[k][i] * ppM[k][i]));
            }

            // mu_kj = mu_kj - [mu_kj]
            ppM[k][j] -= closest_integer_RR;
            ppDelta[k][j] = ((pB[k]/pB[j]) + (ppM[k][j] * ppM[k][j]));
	 }
      }
   }

   // Find the index where there is a potential swap
   int k = 0;
   int i = 0;

   P_min = RR_find_insertion_indices_PotGGLLL (ppM, ppDelta, ppBasis, pB, n, m, &k, &i);
   #ifdef __DEBUG__
   printf ("\n\n\nPotGGLLL: P_min=%Lf\n\n\n", P_min);
   #endif

   while (P_min < delta_threshold) {

      #ifdef __DEBUG__
      printf ("Insertion at k=%d, i=%d for delta=%Lf\n", k, i, P_min);
      #endif

      vec_ZZ pTemp;
      pTemp = ppBasis[k];

      for (int j=k-1; j>=i; j--) {
         ppBasis[j+1] = ppBasis[j];
      }
      ppBasis[i] = pTemp;

         // Update values of ppM, pB
      deep_LLL_GSO_update_RR (ppM, pB, k , i, n);


      // Count the number of swaps
      (*pSwaps)++;
      // if (*pSwaps > 500) {
            //break;
       //}


      // Reduce every vector b_{l}, i <= l < n, by all its previous vectors b_{k-1} ... b_{l-1}
      // Full size reduction or partial size reduction here?
      //for (int l=1; l<n; l++)
      for (int l=i; l<n; l++)
      {
         // Reduce b_l with all the previous basis vectors
         // for (int j=k-1; j<=l-1; j++)
         // for (int j=l-1; j>=k-1; j--)
         for (int j=l-1; j>=0; j--)
         {
            RR closest_int_mu_lj_RR;
	    round(closest_int_mu_lj_RR, ppM[l][j]);
            ZZ closest_integer_to_mu_lj;
	    conv(closest_integer_to_mu_lj, closest_int_mu_lj_RR);

	    if (closest_integer_to_mu_lj != 0) {
               // Count the number of reductions
               (*pReductions)++;
               // Reduce b_l with b_j
               NTL_reduce (ppBasis[j], ppBasis[l], m, closest_integer_to_mu_lj);

               // mu_li -= CLOSEST_INTEGER (mu_lj)
               for (int i=j-1; i>=0; i--) {
                  // By Exercise 17.4.4 from Galbraith v2
                  ppM[l][i] -= closest_int_mu_lj_RR * ppM[j][i];
                  ppDelta[l][i] = ((pB[l]/pB[i]) + (ppM[l][i] * ppM[l][i]));
               }

               // mu_lj = mu_lj - [mu_lj]
               ppM[l][j] -= closest_int_mu_lj_RR;
               ppDelta[l][j] = ((pB[l]/pB[j]) + (ppM[l][j] * ppM[l][j]));
	    }
         }
      }

      #ifdef __DEBUG__
      printf ("\nThe Basis vectors are:\n");
      cout << ppBasis;
      printf ("\nThe values of mu_ij are:");
      for (int i=0; i<n; i++) {
         for (int j=0; j<i; j++) {
            cout << ppM[i][j] << endl << ppDelta[i][j] << endl;
         }
         printf ("\n");
      }
      #endif
      // Update k
      P_min = RR_find_insertion_indices_PotGGLLL (ppM, ppDelta, ppBasis, pB, n, m, &k, &i);
   }

   ppM.kill();
   ppDelta.kill();
   ppBasisGSO.kill();
   
   return ppBasis;
}

// The SS-GGLLL implementation with NTL datatypes
// At every iteration, find indices k, i such that SS(B) - SS(C) is maximised (i.e. max reduction of SS over all indices).
// When this happens we insert b_k in position i
mat_ZZ NTL_SSGGLLL (mat_ZZ ppBasis, int m, int n, RR eta, long long int * pSwaps, long long int * pReductions) {
   RR eta_p;
   eta_p = 1 - eta;

   mat_RR ppBasisGSO, ppBasisRR;
   ppBasisGSO.SetDims(n,m);
   ppBasisRR.SetDims(n,m);
   // Find the GSO of the Basis
   ppBasisGSO = GSO_RR (ppBasis, m, n);
   ppBasisRR = copyNTLBasisToRR(ppBasis, m, n);

   #ifdef __DEBUG__
   printf ("\n\n\nThe GSO of the basis is:\n");
   cout << ppBasisGSO;
   #endif

   // This 1-dimensional array pB will store all the values of ||b_i*||^2
   vec_RR pB;
   pB.SetLength(n);
   #ifdef __DEBUG__
   printf ("\nThe values of ||b_i*||^2 are:\n");
   #endif

   for (int i=0; i<n; i++) {
      InnerProduct(pB[i], ppBasisGSO[i], ppBasisGSO[i]);
      #ifdef __DEBUG__
      cout << pB[i] << endl;
      #endif
   }

   // This 2-dimensional nxn array ppM will store all the values of mu_ij
   mat_RR ppM;
   ppM.SetDims(n,m);
   #ifdef __DEBUG__
   printf ("\n\nThe values of mu_ij are:");
   #endif
   
   for (int i=1; i<n; i++) {
      for (int j=0; j<i; j++) {
	 RR inner;
	 InnerProduct(inner, ppBasisRR[i], ppBasisGSO[j]);
         div(ppM[i][j], inner, pB[j]);
         #ifdef __DEBUG__
         cout << ppM[i][j] << endl << ppDelta[i][j] << endl;
         #endif
      }
      #ifdef __DEBUG__
      printf ("\n");
      #endif
   }

   // The greedy global LLL task
   RR S_max, SS_sub_S_max;
   S_max = 0;
   RR squared_sum;
   squared_sum = 0;
   for (int i=0; i<n; i++) {
      squared_sum += pB[i];
   }
   SS_sub_S_max = squared_sum;
  
   // First reduce every vector by all its previous vectors
   for (int k=1; k<n; k++) {

      // Reduce b_k with all the previous basis vectors
      for (int j=k-1; j>=0; j--) {
         RR closest_integer_RR;
         round(closest_integer_RR, ppM[k][j]);
	 ZZ closest_integer_to_mu_kj;
         conv(closest_integer_to_mu_kj, closest_integer_RR);	 

	 if (closest_integer_to_mu_kj != 0) {
            // Count the number of reductions
            (*pReductions)++;
            // Reduce b_k with b_j
            NTL_reduce (ppBasis[j], ppBasis[k], m, closest_integer_to_mu_kj);
		
            // mu_ki -= CLOSEST_INTEGER (mu_kj)
            for (int i=j-1; i>=0; i--) {
               // By Exercise 17.4.4 from Galbraith v2
               ppM[k][i] -= closest_integer_RR * ppM[j][i];
            }

            // mu_kj = mu_kj - [mu_kj]
            ppM[k][j] -= closest_integer_RR;
         }
      }
   }
   // Find the index where there is a potential swap
   int k = 0;
   int i = 0;
   //cout << "Basis After Size-Reduction: " << endl << ppBasis << endl;
   S_max = RR_find_insertion_indices_SSGGLLL (ppM, ppBasis, pB, n, m, &k, &i);
   //#ifdef __DEBUG__
   // cout << "SS(B) - SS(C) = " << S_max << endl;
   //#endif

   RR SS_recomputed;
   while (S_max > eta_p * squared_sum) {

      //#ifdef __DEBUG__
      //cout << "Insertion at k = " << k << "  i = " << i << "  SS(B)-SS(C) = " << S_max << endl;
      //#endif

      vec_ZZ pTemp;
      pTemp = ppBasis[k];

      for (int j=k-1; j>=i; j--) {
         ppBasis[j+1] = ppBasis[j];
      }
      ppBasis[i] = pTemp;

      // Update values of ppM, pB
      deep_LLL_GSO_update_RR (ppM, pB, k , i, n);

      // Update squared_sum
      SS_sub_S_max -= S_max;
      SS_recomputed = 0;
      for (int i=0; i<n; i++) {
         SS_recomputed += pB[i];
      }
      squared_sum = SS_recomputed;
      //cout << "Updated SS (SS(B) - S_max) = " << SS_sub_S_max << endl;
      //cout << "Updated SS (recomputed from GSO) = " << SS_recomputed << endl;
      //cout << "Updated Basis after Insertion: " << endl << ppBasis << endl;
      // Count the number of swaps
      (*pSwaps)++;


      // Reduce every vector b_{l}, i <= l < n, by all its previous vectors b_{k-1} ... b_{l-1}
      // Full size reduction or partial size reduction here?
      //for (int l=1; l<n; l++)
      for (int l=i; l<n; l++)
      {
         // Reduce b_l with all the previous basis vectors
         // for (int j=k-1; j<=l-1; j++)
         // for (int j=l-1; j>=k-1; j--)
         for (int j=l-1; j>=0; j--) {
            RR closest_integer_mu_lj_RR;
	    round(closest_integer_mu_lj_RR, ppM[l][j]);
            ZZ closest_integer_to_mu_lj;
            conv(closest_integer_to_mu_lj, closest_integer_mu_lj_RR);

	    if (closest_integer_to_mu_lj != 0) {
               // Count the number of reductions
               (*pReductions)++;
               // Reduce b_l with b_j
               NTL_reduce (ppBasis[j], ppBasis[l], m, closest_integer_to_mu_lj);

               // mu_li -= CLOSEST_INTEGER (mu_lj)
               for (int i=j-1; i>=0; i--) {
                  // By Exercise 17.4.4 from Galbraith v2
                  ppM[l][i] -= closest_integer_mu_lj_RR * ppM[j][i];
               }

               // mu_lj = mu_lj - [mu_lj]
               ppM[l][j] -= closest_integer_mu_lj_RR;
	    }
         }
      }
      //cout << "Basis after Size-Reduction: " << endl << ppBasis << endl;
      #ifdef __DEBUG__
      printf ("\nThe Basis vectors are:\n");
      cout << ppBasis << endl;
      printf ("\nThe values of mu_ij are:");
      for (int i=0; i<n; i++) {
         for (int j=0; j<i; j++) {
            cout << ppM[i][j] << endl;
         }
         printf ("\n");
      }
      #endif

      // Update k
     // ppBasisInt = copyNTLBasisToInt(ppBasis, m, n);
      S_max = RR_find_insertion_indices_SSGGLLL (ppM, ppBasis, pB, n, m, &k, &i);
   }

   ppM.kill();
   ppBasisGSO.kill();
   ppBasisRR.kill();
   
   return ppBasis;
}

// The SS-GGLLL implementation with standard C datatypes
// At every iteration, find indices k, i such that SS(B) - SS(C) is maximised (i.e. max reduction of SS over all indices).
// When this happens we insert b_k in position i
int ** SSGGLLL_std (int ** ppBasis, int m, int n, long double eta, long long int * pSwaps, long long int * pReductions) {
   	long double eta_p = 0;
   	eta_p = 1 - eta;

	long double ** ppBasisGSO;
	ppBasisGSO = GSO (ppBasis, m, n);
	long double ** ppBasisDouble;
	ppBasisDouble = copyBasisToDouble (ppBasis, m, n);
	long double ** ppM;
	ppM = (long double **) calloc (n, sizeof(long double *));
   	for (int i=0; i<n; i++) {
   		ppM[i] = (long double *) calloc (m, sizeof(long double));
	}

   	#ifdef __DEBUG__
   	printf ("\n\n\nThe GSO of the basis is:\n");
	printBasisDoubleAsColumns(ppBasisGSO, m, n);
   	#endif

   	// This 1-dimensional array pB will store all the values of ||b_i*||^2
   	long double * pB;
   	pB = (long double *) calloc (m, sizeof(long double));

   	#ifdef __DEBUG__
   	printf ("\nThe values of ||b_i*||^2 are:\n");
   	#endif

   	for (int i=0; i<n; i++) {
     		pB[i] = inner_product_template(ppBasisGSO[i], ppBasisGSO[i], m);
      		#ifdef __DEBUG__
      		cout << pB[i] << endl;
      		#endif
   	}

   	#ifdef __DEBUG__
   	printf ("\n\nThe values of mu_ij are:");
   	#endif

	long double inner = 0;

   	for (int i=1; i<n; i++) {
      		for (int j=0; j<i; j++) {
	 		inner = inner_product_template(ppBasisDouble[i], ppBasisGSO[j], m);
         		ppM[i][j] = inner / pB[j];

         		#ifdef __DEBUG__
         		cout << ppM[i][j] << " ";
         		#endif
      		}
      	#ifdef __DEBUG__
     	printf ("\n");
      	#endif
   	}

   	// The greedy global LLL task
   	long double S_max = 0;
        long double SS_sub_S_max = 0;
   	long double squared_sum = 0;

	long double bound = 0;

	for (int i=0; i<n; i++) {
      		squared_sum += pB[i];
   	}

	SS_sub_S_max = squared_sum;

	int closest_integer = 0;
	long double closest_integer_Lf = 0.0;
	// First reduce every vector by all its previous vectors
   	for (int k=1; k<n; k++) {

      		// Reduce b_k with all the previous basis vectors
      		for (int j=k-1; j>=0; j--) {
         		closest_integer = CLOSEST_INTEGER(ppM[k][j]);
			closest_integer_Lf = (long double) closest_integer;
	 		if (closest_integer != 0) {
            			// Count the number of reductions
            			(*pReductions)++;
            			// Reduce b_k with b_j
            			reduce (ppBasis[j], ppBasis[k], m, closest_integer);

            			// mu_ki -= CLOSEST_INTEGER (mu_kj)
            			for (int i=j-1; i>=0; i--) {
               				// By Exercise 17.4.4 from Galbraith v2
               				ppM[k][i] -= closest_integer_Lf * ppM[j][i];
            			}

            			// mu_kj = mu_kj - [mu_kj]
            			ppM[k][j] -= closest_integer_Lf;
         		}
      		}
   	}

	// Find the index where there is a potential swap
 	int k = 0;
   	int i = 0;

	S_max = find_insertion_indices_SSGGLLL (ppM, ppBasis, pB, n, m, &k, &i);

	bound = eta_p * squared_sum;
	long double SS_recomputed = 0;

      	int * pTemp;
   	pTemp = (int *) calloc (m, sizeof(int));


   	while (S_max > bound) {

      		pTemp = ppBasis[k];

      		for (int j=k-1; j>=i; j--) {
         		ppBasis[j+1] = ppBasis[j];
      		}

		ppBasis[i] = pTemp;

      		// Update values of ppM, pB
      		deep_LLL_GSO_update (ppM, pB, k , i, n);

      		// Update squared_sum
      		SS_sub_S_max -= S_max;
      		SS_recomputed = 0;
      		for (int i=0; i<n; i++) {
        		SS_recomputed += pB[i];
      		}
      		squared_sum = SS_recomputed;

		// Count the number of swaps
      		(*pSwaps)++;

      		// Reduce every vector b_{l}, i <= l < n, by all its previous vectors b_{k-1} ... b_{l-1}
      		// Full size reduction or partial size reduction here?
      		//for (int l=1; l<n; l++)
      		for (int l=i; l<n; l++) {
        	 // Reduce b_l with all the previous basis vectors
        		for (int j=l-1; j>=0; j--) {
	    			closest_integer = CLOSEST_INTEGER(ppM[l][j]);
				closest_integer_Lf = (long double) closest_integer;
	    			if (closest_integer != 0) {
               				// Count the number of reductions
               				(*pReductions)++;
               				// Reduce b_l with b_j
               				reduce (ppBasis[j], ppBasis[l], m, closest_integer);

               				// mu_li -= CLOSEST_INTEGER (mu_lj)
               				for (int i=j-1; i>=0; i--) {
                  				// By Exercise 17.4.4 from Galbraith v2
                  				ppM[l][i] -= closest_integer_Lf * ppM[j][i];
               				}

        	       			// mu_lj = mu_lj - [mu_lj]
               				ppM[l][j] -= closest_integer_Lf;
	    			}
         		}
      		}

		//cout << "Basis after Size-Reduction: " << endl << ppBasis << endl;

		#ifdef __DEBUG__
      		printf ("\nThe Basis vectors are:\n");
      		cout << ppBasis << endl;
      		printf ("\nThe values of mu_ij are:");
      		for (int i=0; i<n; i++) {
         		for (int j=0; j<i; j++) {
            			cout << ppM[i][j] << endl;
         		}
         		printf ("\n");
      		}
      		#endif

      		// Update k
      		S_max = find_insertion_indices_SSGGLLL (ppM, ppBasis, pB, n, m, &k, &i);
		bound = eta_p * squared_sum;
   	}
   	deleteBasisDouble(ppM, n);
   	deleteBasisDouble(ppBasisGSO, n);
   	deleteBasisDouble(ppBasisDouble, n);

   	return ppBasis;
}

// The SS-LLL implementation with standard C datatypes
// At every iteration, find insert the vector b_k position i such that the value of S^2 reduces as much as possible
int ** SS_LLL_std (int ** ppBasis, int m, int n, long double eta, long long int * pSwaps, long long int * pReductions) {

   	long double eta_p = 0;
   	eta_p = 1 - eta;

	long double ** ppBasisGSO;
	ppBasisGSO = GSO (ppBasis, m, n);
	long double ** ppBasisDouble;
	ppBasisDouble = copyBasisToDouble (ppBasis, m, n);
	long double ** ppM;
	ppM = (long double **) calloc (n, sizeof(long double *));
   	for (int i=0; i<n; i++) {
   		ppM[i] = (long double *) calloc (m, sizeof(long double));
	}

   	#ifdef __DEBUG__
   	printf ("\n\n\nThe GSO of the basis is:\n");
	printBasisDoubleAsColumns(ppBasisGSO, m, n);
   	#endif

   	// This 1-dimensional array pB will store all the values of ||b_i*||^2
   	long double * pB;
   	pB = (long double *) calloc (m, sizeof(long double));

   	#ifdef __DEBUG__
   	printf ("\nThe values of ||b_i*||^2 are:\n");
   	#endif

   	for (int i=0; i<n; i++) {
     		pB[i] = inner_product_template(ppBasisGSO[i], ppBasisGSO[i], m);
      		#ifdef __DEBUG__
      		cout << pB[i] << endl;
      		#endif
   	}

   	#ifdef __DEBUG__
   	printf ("\n\nThe values of mu_ij are:");
   	#endif

	long double inner = 0;

   	for (int i=1; i<n; i++) {
      		for (int j=0; j<i; j++) {
	 		inner = inner_product_template(ppBasisDouble[i], ppBasisGSO[j], m);
         		ppM[i][j] = inner / pB[j];

         		#ifdef __DEBUG__
         		cout << ppM[i][j] << " ";
         		#endif
      		}
      	#ifdef __DEBUG__
     	printf ("\n");
      	#endif
   	}

	// The S^2LLL task
   	*pSwaps = 0;
   	*pReductions = 0;
   	int k = 1;
	int closest_integer = 0;

	long double delta_k_1 = 0.0;
	long double D_j = 0.0;
	long double S_ik = 0.0;
	long double S_lk = 0.0;

	int * pTemp;
   	pTemp = (int *) calloc (m, sizeof(int));

   	while (k<n) {
      		long double squared_sum = 0;

		for (int i=0; i<n; i++) {
         		squared_sum += pB[i];
      		}

      		// Reduce b_k with all previous basis vectors
      		for (int j=k-1; j>=0; j--) {
	 		closest_integer = CLOSEST_INTEGER(ppM[k][j]);

	 		if (closest_integer != 0) {
            			// Count the number of reductions
            			(*pReductions)++;
            			// Reduce b_k with b_j
            			reduce (ppBasis[j], ppBasis[k], m, closest_integer);

            			// mu_ki -= CLOSEST_INTEGER (mu_kj);
            			for (int i = j-1; i>=0; i--) {
               				// By Exercise 17.4.4 from Galbraith v2
               				ppM[k][i] -= closest_integer * ppM[j][i];
            			}

            			// mu_kj = mu_kj - [mu_kj]
            			ppM[k][j] -= closest_integer;
         		}
     	 	}

      		// For j = k-1,...,1,0, we find the index i where we insert b_k to reduce the value of SS the most
      		delta_k_1 = (pB[k] + (ppM[k][k-1] * ppM[k][k-1] * pB[k-1])) / pB[k-1];

      		// Initialise D_j = D_{k-1}
      		D_j = (ppM[k][k-1] * ppM[k][k-1] * pB[k-1]) + pB[k];

      		// Initialise S_ik = S_{i-1,k}
      		// After initialising S and D, we can update them in the 'for' loop below for general S_jk and D_j for 0 <= j < k-1
      		S_ik = ppM[k][k-1] * ppM[k][k-1] * (1 - delta_k_1);

      		int i = k-1;

      		S_lk = S_ik;

      		for (int l=k-2; l>=0; l--) {
        		D_j += ppM[k][l] * ppM[k][l] * pB[l];
         		S_lk += ppM[k][l] * ppM[k][l] * pB[l] * ((pB[l]/D_j) - 1);

         		if (S_lk > S_ik) {
            			i = l;
            			S_ik = S_lk;
         		}
      		}

      		if (S_ik <= eta_p * squared_sum) {
         		k++;
      		} else {
         		// b_l <- b_k, and vectors b_l, ... , b_k-1 are shifted up one position
         		pTemp = ppBasis[k];

         		for (int j=k-1; j>=i; j--) {
            			ppBasis[j+1] = ppBasis[j];
         		}

			ppBasis[i] = pTemp;

         		// Update values of ppM, pB
         		(*pSwaps)++;

         		deep_LLL_GSO_update (ppM, pB, k , i, n);
         		k = (1 >(i)) ? 1 : i;
      		}
   	}
  	deleteBasisDouble(ppM, n);
  	deleteBasisDouble(ppBasisGSO, n);
   	deleteBasisDouble(ppBasisDouble, n);

   	return ppBasis;
}

// The PotLLL implementation with standard C datatypes
// At every iteration, find insert the vector b_k position i such that the basis potential reduces by a sufficiently large amount
int ** Pot_LLL_std (int ** ppBasis, int m, int n, long double delta, long long int * pSwaps, long long int * pReductions) {

	long double ** ppBasisGSO;
	ppBasisGSO = GSO (ppBasis, m, n);
	long double ** ppBasisDouble;
	ppBasisDouble = copyBasisToDouble (ppBasis, m, n);
	long double ** ppM;
	ppM = (long double **) calloc (n, sizeof(long double *));
   	for (int i=0; i<n; i++) {
   		ppM[i] = (long double *) calloc (m, sizeof(long double));
	}

   	#ifdef __DEBUG__
   	printf ("\n\n\nThe GSO of the basis is:\n");
	printBasisDoubleAsColumns(ppBasisGSO, m, n);
   	#endif

   	// This 1-dimensional array pB will store all the values of ||b_i*||^2
   	long double * pB;
   	pB = (long double *) calloc (m, sizeof(long double));

   	#ifdef __DEBUG__
   	printf ("\nThe values of ||b_i*||^2 are:\n");
   	#endif

   	for (int i=0; i<n; i++) {
     		pB[i] = inner_product_template(ppBasisGSO[i], ppBasisGSO[i], m);
      		#ifdef __DEBUG__
      		cout << pB[i] << endl;
      		#endif
   	}

   	#ifdef __DEBUG__
   	printf ("\n\nThe values of mu_ij are:");
   	#endif

	long double inner = 0;

   	for (int i=1; i<n; i++) {
      		for (int j=0; j<i; j++) {
	 		inner = inner_product_template(ppBasisDouble[i], ppBasisGSO[j], m);
         		ppM[i][j] = inner / pB[j];

         		#ifdef __DEBUG__
         		cout << ppM[i][j] << " ";
         		#endif
      		}
      	#ifdef __DEBUG__
     	printf ("\n");
      	#endif
   	}

   	// The PotLLL task
   	*pSwaps = 0;
   	*pReductions = 0;
   	int k = 1;
   	long double C = 0.0;
   	int temp = 0;
	int closest_integer = 0;
	long double closest_integer_Lf = 0.0;

	long double P = 1.0;
	long double P_min = 1.0;
	int l = 1;
	long double sum = 0.0;
	long double P_mult = 0.0;
	int * pTemp;
   	pTemp = (int *) calloc (m, sizeof(int));

   	while (k<n) {
      		// Reduce b_k with all previous basis vectors
      		for (int j=k-1; j>=0; j--) {
	 		closest_integer = CLOSEST_INTEGER(ppM[k][j]);
        		closest_integer_Lf = (long double) closest_integer;
	 		if (closest_integer != 0) {
            			// Count the number of reductions
            			(*pReductions)++;
            			// Reduce b_k with b_j
            			reduce (ppBasis[j], ppBasis[k], m, closest_integer);

            			// mu_ki -= CLOSEST_INTEGER (mu_kj);
            			for (int i = j-1; i>=0; i--) {
               				// By Exercise 17.4.4 from Galbraith v2
               				ppM[k][i] -= closest_integer_Lf * ppM[j][i];
            			}

            			// mu_kj = mu_kj - [mu_kj]
            			ppM[k][j] -= closest_integer_Lf;
         		}
      		}

      		// Initialise P = 1, P_min = 1, k = 0
      		// For j = k-1,...,1,0, we find the index where an insertion of vector b_k in position j reduces the potential the most
      		P = 1;
      		P_min = 1;
      		l = 1;

      		for (int j=k-1; j>=0; j--) {
	 		sum = 0.0;
         		for (int r=j; r<k; r++) {
            			sum += ppM[k][r]*ppM[k][r]*pB[r];
         		}

         		P_mult = (pB[k] + sum) / pB[j];
         		P = P * P_mult;
         		if (P<P_min) {
            			l = j;
            			P_min = P;
         		}
      		}

      		#ifdef __DEBUG__
      		printf ("\nP_min = %Lf", P_min);
      		#endif

      		if (delta > P_min) {
         		// b_l <- b_k, and vectors b_l, ... , b_k-1 are shifted up one position
         		pTemp = ppBasis[k];

         		for (int j=k-1; j>=l; j--) {
            			ppBasis[j+1] = ppBasis[j];
         		}
         		ppBasis[l] = pTemp;

         		// Update values of ppM[i][j] for j<i, pB[i]
         		(*pSwaps)++;

         		deep_LLL_GSO_update (ppM, pB, k , l, n);

         		k = (1 > (l)) ? 1 : l;
         	} else {
         	k++;
      		}
   	}
   	deleteBasisDouble(ppM, n);
   	deleteBasisDouble(ppBasisGSO, n);
   	deleteBasisDouble(ppBasisDouble, n);

   	return ppBasis;
}

// The Pot-GGLLL implementation with standard C datatypes
// At every iteration, find indices k, i such that the Potential of the basis is reduced by the largest amount.
// When this happens we insert b_k in position i
int ** PotGGLLL_std (int ** ppBasis, int m, int n, long double delta_threshold, long long int * pSwaps, long long int * pReductions) {

	long double ** ppBasisGSO;
	ppBasisGSO = GSO (ppBasis, m, n);
	long double ** ppBasisDouble;
	ppBasisDouble = copyBasisToDouble (ppBasis, m, n);
	long double ** ppM;
	ppM = (long double **) calloc (n, sizeof(long double *));
   	for (int i=0; i<n; i++) {
   		ppM[i] = (long double *) calloc (m, sizeof(long double));
	}


   	#ifdef __DEBUG__
 	printf ("\n\n\nThe GSO of the basis is:\n");
 	cout << ppBasisGSO;
 	#endif

   	// This 1-dimensional array pB will store all the values of ||b_i*||^2
   	long double * pB;
   	pB = (long double *) calloc (m, sizeof(long double));

   	#ifdef __DEBUG__
   	printf ("\nThe values of ||b_i*||^2 are:\n");
   	#endif

	for (int i=0; i<n; i++) {
     		pB[i] = inner_product_template(ppBasisGSO[i], ppBasisGSO[i], m);
      		#ifdef __DEBUG__
      		cout << pB[i] << endl;
      		#endif
	}

   	#ifdef __DEBUG__
  	 printf ("\n\nThe values of mu_ij are:");
   	#endif

  	for (int i=1; i<n; i++) {
      		for (int j=0; j<i; j++) {
	 		long double inner;
         		inner = inner_product_template(ppBasisDouble[i], ppBasisGSO[j], m);
	 		ppM[i][j] = inner / pB[j];
         		#ifdef __DEBUG__
         		cout << ppM[i][j] << endl;
         		#endif
      		}
      		#ifdef __DEBUG__
      		printf ("\n");
      		#endif
   	}

	// The greedy global LLL task
   	*pSwaps = 0;
   	*pReductions = 0;
   	long double P_min = 0;
        int closest_integer = 0;
	long double closest_integer_Lf = 0.0;

   	// First reduce every vector by all its previous vectors
   	for (int k=1; k<n; k++) {

   	   // Reduce b_k with all the previous basis vectors
   		for (int j=k-1; j>=0; j--) {
	 		closest_integer = CLOSEST_INTEGER(ppM[k][j]);
        		closest_integer_Lf = (long double) closest_integer;

		 	if (closest_integer != 0) {
         			// Count the number of reductions
            			(*pReductions)++;
            			// Reduce b_k with b_j
            			reduce (ppBasis[j], ppBasis[k], m, closest_integer);

            			// mu_ki -= CLOSEST_INTEGER (mu_kj)
            			for (int i=j-1; i>=0; i--) {
               				// By Exercise 17.4.4 from Galbraith v2
               				ppM[k][i] -= closest_integer_Lf * ppM[j][i];
            			}

            			// mu_kj = mu_kj - [mu_kj]
            			ppM[k][j] -= closest_integer_Lf;
	 		}
      		}
   	}

   	// Find the index where there is a potential swap
   	int k = 0;
   	int i = 0;

	int * pTemp;
   	pTemp = (int *) calloc (m, sizeof(int));


   	P_min = find_insertion_indices_PotGGLLL (ppM, ppBasis, pB, n, m, &k, &i);
   	#ifdef __DEBUG__
   	printf ("\n\n\nPotGGLLL: P_min=%Lf\n\n\n", P_min);
   	#endif

   	while (P_min < delta_threshold) {

      		#ifdef __DEBUG__
      		printf ("Insertion at k=%d, i=%d for delta=%Lf\n", k, i, P_min);
      		#endif

      		pTemp = ppBasis[k];

      		for (int j=k-1; j>=i; j--) {
         		ppBasis[j+1] = ppBasis[j];
      		}

		ppBasis[i] = pTemp;

         	// Update values of ppM, pB
      		deep_LLL_GSO_update (ppM, pB, k , i, n);


      		// Count the number of swaps
      		(*pSwaps)++;


      		// Reduce every vector b_{l}, i <= l < n, by all its previous vectors b_{k-1} ... b_{l-1}
      		// Full size reduction or partial size reduction here?

		for (int l=i; l<n; l++) {

			// Reduce b_l with all the previous basis vectors
         		// for (int j=k-1; j<=l-1; j++)
         		// for (int j=l-1; j>=k-1; j--)

			for (int j=l-1; j>=0; j--) {
	    			closest_integer = CLOSEST_INTEGER(ppM[l][j]);
        			closest_integer_Lf = (long double) closest_integer;

	    			if (closest_integer != 0) {
               				// Count the number of reductions
               				(*pReductions)++;
               				// Reduce b_l with b_j
               				reduce (ppBasis[j], ppBasis[l], m, closest_integer);

               				// mu_li -= CLOSEST_INTEGER (mu_lj)
               				for (int i=j-1; i>=0; i--) {
                  				// By Exercise 17.4.4 from Galbraith v2
                  				ppM[l][i] -= closest_integer_Lf * ppM[j][i];
               				}

               				// mu_lj = mu_lj - [mu_lj]
               				ppM[l][j] -= closest_integer_Lf;
	    			}
         		}
      		}

      		#ifdef __DEBUG__
      		printf ("\nThe Basis vectors are:\n");
      		cout << ppBasis;
      		printf ("\nThe values of mu_ij are:");
      		for (int i=0; i<n; i++) {
         		for (int j=0; j<i; j++) {
            			cout << ppM[i][j] << endl;
         		}

			printf ("\n");
      		}
      		#endif

		// Update k
      		P_min = find_insertion_indices_PotGGLLL (ppM, ppBasis, pB, n, m, &k, &i);
   	}
	deleteBasisDouble (ppM, n);
	deleteBasisDouble (ppBasisGSO, n);
	deleteBasisDouble (ppBasisDouble, n);

   	return ppBasis;
}
