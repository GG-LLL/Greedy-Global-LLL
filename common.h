/*******************************************************************************
 * Header File: common.h
*******************************************************************************/
#include <fstream>
#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <gmp.h>
#include <mpfr.h>
#include <NTL/mat_ZZ.h>
#include <NTL/LLL.h>
#include <NTL/fileio.h>
#include <NTL/vec_double.h>
#include <NTL/ZZ.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>
#include <NTL/xdouble.h>
#include <NTL/quad_float.h>
#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <bits/stdc++.h>

#define _PARENT_DIR_ "../generate_SVP_bases"
#define _PARENT_INPUT_DIR_ "../generate_SVP_bases"
#define _PARENT_OUTPUT_DIR_ "../SVP_Outputs"
#define _PARENT_OUTPUT_DIR_PREPROCESSED_ "../SVP_Outputs_Preprocessed"
#define _INPUT_DIR_ "SVP_Bases"
#define _INPUT_DIR_FPLLL_RED_ "SVP_Bases_LLLReduced"
#define _INPUT_FILE_PREFIX_ "SVP_Basis"
#define _OUTPUT_FILE_PREFIX_ "out"

#define _LLL_         "LLL"
#define _SS_LLL_      "SS_LLL"
#define _SS_GGLLL_    "SS_GGLLL"
#define _Pot_LLL_     "Pot_LLL"
#define _Pot_GGLLL_   "Pot_GGLLL"

using namespace std;
using namespace NTL;

#define TRUE 1
#define FALSE 0

#define FILE_ERROR -1

#define NO_BITS_INT (sizeof(unsigned int)*8)
#define NO_BITS_DOUBLE (sizeof(long double)*8)
#define MASK_MSB_INT ((unsigned int)(1<<(NO_BITS_INT-1)))
#define MASK_MSB_DOUBLE ((unsigned int)(1<<(NO_BITS_DOUBLE-1)))
#define ALL_ONE_INT ((unsigned int)(((MASK_MSB_INT-1)<<1)+1))
#define ALL_ONE_DOUBLE ((unsigned int)(((MASK_MSB_DOUBLE-1)<<1)+1))o

#define CLOSEST_INTEGER(X) ( (X>0) ?  ((long double)((long double)X - (int)X) > ((long double)0.51))? (int)X + 1: (int)X : ((long double)((long double)X - (int)X) < ((long double)-0.51))? ((int)X - 1): (int)X )

#define LC (Bi, Bj, delta, mu_ij) ((Bi >= ((delta - (mu_ij*mu_ij)) * Bj)) ? TRUE : FALSE)

// Set precision of ZZ and RR
int setPrecision(int m);

// Deallocate ppBasis
int deleteBasis (int ** ppBasis, int n);

// Deallocate ppBasisDouble
template <typename T1>
int deleteBasisDouble (T1 ** ppBasisDouble, int n);

// Compute SVP Challenge factor
RR SVP_Challenge_Factor(int n, ZZ volume);

// Comparator function to sort pairs according to second value
bool cmp_double(pair<string, double>& a, pair<string, double>& b);

// Comparator function to sort pairs (RR) according to second value
bool cmp_long_double(pair<string, long double>& a, pair<string, long double>& b);

// Function to sort the map according
// to value in a (key-value) pairs (fix)
void map_sort_long_double_fix(map<string, long double>& M);

// Function to sort the map according
// to value in a (key-value) pairs (fix)
void map_sort_double_fix(map<string, double>& M);

// Compare elements for use in qsort algorithm
int compare (const void * a, const void * b);

// Generate a random NTL basis (mat_ZZ) with n integer vectors, each of dimension m
mat_ZZ generateNTLChallengeBasis (int n);

// Find the inner products of vectors pV1 with pV2
template <typename T1, typename T2>
long double inner_product_template (T1 * pV1, T2 * pV2, int m);

// Compute norm of a vec_ZZ
RR NTL_vector_norm (vec_ZZ V1, int m);

// Reduce vector pV2 with pV1
int reduce (int * pV1, int * pV2, int m, int closest_integer_to_mu);
template <typename T1>
int reduceDouble (T1 * pVDouble1, T1 * pVDouble2, int m, T1 mu);

// Reduce vec_ZZ V2 with vec_ZZ V1 (only NTL)
int NTL_reduce (vec_ZZ& V1, vec_ZZ& V2, int m, ZZ closest_ZZ_to_mu);

// Copy int **ppBasis to mat_ZZ BasisNTL and return
mat_ZZ copyBasisToNTL (int ** ppBasis, int m, int n);

// Copy mat_ZZ ppBasis to mat_RR BasisNTL and return
mat_RR copyNTLBasisToRR (mat_ZZ ppBasis, int m, int n);

// Copy pBasis into ppBasisDouble and return
long double ** copyBasisToDouble (int ** ppBasis, int m, int n);

//  Rounds a/d to nearest integer, breaking ties
//  by rounding towards zero.  Assumes d > 0.
//  Stored in q
//  Copied from NTL Library (github) 
int BalDiv(ZZ& q, const ZZ& a, const ZZ& d);

// Exact division: q = a/b where b|a
// Copied from NTL Library (github) 
int ExactDiv(ZZ& qq, const ZZ& a, const ZZ& b);

// Gram-Schmidt Orthogonalisation
long double ** GSO (int ** pBasis, int m, int );

// Gram-Schmidt Orthogonalisation using RR datatype
mat_RR GSO_RR (mat_ZZ ppBasis, int m, int n);

// Update the GSO information in the way required for DeepLLL
// Algorithm 4 in Yamaguchi Yasuda NuTMic paper (NuTMic 2017)
template <typename T1>
int deep_LLL_GSO_update (T1 ** ppM, T1 * pB, int k , int i, int n);

// Update the GSO information in the way required for DeepLLL (RR datatypes)
// Algorithm 4 in Yamaguchi Yasuda NuTMic paper (NuTMic 2017)
int deep_LLL_GSO_update_RR (mat_RR &ppM, vec_RR &pB, int k , int i, int n);

// Find the insertion indices for SS-GGLLL
RR RR_find_insertion_indices_SSGGLLL (mat_RR ppM, mat_ZZ ppBasis, vec_RR pB, int n, int m, int * pK, int * pJ);
template <typename T1>
T1 find_insertion_indices_SSGGLLL (T1 ** ppM, int ** ppBasis, T1 * pB, int n, int m, int * pK, int * pJ);

//Find the insertion indices for Pot-GGLLL
RR RR_find_insertion_indices_PotGGLLL (mat_RR ppM, mat_RR ppDelta, mat_ZZ ppBasis, vec_RR pB, int n, int m, int * pK, int * pJ);
long double find_insertion_indices_PotGGLLL (long double ** ppM, int ** ppBasis, long double * pB, int n, int m, int * pK, int * pJ);

// Compute Root Hermite Factor (RR)
RR NTLrootHermiteFactor (mat_ZZ ppBasis, int m, int n); 

// Compute orthogonality defect (RR)
RR NTLorthogonalityDefect (mat_ZZ ppBasis, int m, int n);

// Compute squared sum (RR)
RR NTLsquaredSum (mat_ZZ ppBasis, int m, int n);

// Compute log of the potential (RR)
RR NTLbasisPotential (mat_ZZ ppBasis, int m, int n);

// Check if a Basis is LLL reduced
int NTL_LLL_check (mat_ZZ ppBasis, int m, int n, RR delta_threshold);

// The standard LLL algorithm using NTL mat_ZZ datatype for basis
mat_ZZ LLL_mat_ZZ (mat_ZZ ppBasis, int m, int n, RR delta, long long int * pSwaps, long long int * pReductions);

// The PotLLL implementation with NTL datatypes / standard C datatypes
// At every iteration, find and insert the vector b_k in position i such that the basis potential reduces by a sufficiently large amount
mat_ZZ NTL_PotLLL (mat_ZZ ppBasis, int m, int n, RR delta, long long int * pSwaps, long long int * pReductions);
int ** Pot_LLL_std (int ** ppBasis, int m, int n, long double delta, long long int * pSwaps, long long int * pReductions);

// The S^2LLL implementation using NTL datatypes / standard C datatypes
// At every iteration, find and insert the vector b_k in position i such that the value of S^2 reduces as much as possible
mat_ZZ NTL_SS_LLL (mat_ZZ ppBasis, int m, int n, RR eta, long long int * pSwaps, long long int * pReductions);
int ** SS_LLL_std (int ** ppBasis, int m, int n, long double eta, long long int * pSwaps, long long int * pReductions);

// The Pot-GGLLL implementations with NTL datatypes / standard C datatypes
// At every iteration, find indices k, i such that the Potential of the basis is reduced by the largest amount.
// When this happens we insert b_k in position i
mat_ZZ NTL_PotGGLLL (mat_ZZ ppBasis, int m, int n, RR delta_threshold, long long int * pSwaps, long long int * pReductions);
int ** PotGGLLL_std (int ** ppBasis, int m, int n, long double delta_threshold, long long int * pSwaps, long long int * pReductions);

// The SS-GGLLL implementation with NTL datatypes / standard C datatypes
// At every iteration, find indices k, i such that SS(B) - SS(C) is maximised (i.e. max reduction of SS over all indices).
// When this happens we insert b_k in position i
mat_ZZ NTL_SSGGLLL (mat_ZZ ppBasis, int m, int n, RR eta, long long int * pSwaps, long long int * pReductions);
int ** SSGGLLL_std (int ** ppBasis, int m, int n, long double eta, long long int * pSwaps, long long int * pReductions);
