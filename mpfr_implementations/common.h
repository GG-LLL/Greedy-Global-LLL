/* 
 * *****************************************************************************
 * Filename: common.h
 * Authors : Sanjay Bhattacherjee and Jack Moyler
 * *****************************************************************************
*/

#ifndef COMMON_H
#define COMMON_H

#include <cstring>
#include "test_utils.h"
#include <fplll.h>

using namespace std;
using namespace fplll;

#define __DEBUG__

#define DEFAULT_PRECISION 53
#define DEFAULT_SUBLATTICE_SIZE 0
#define DEFAULT_SEARCH_SWITCH 0
#define DEFAULT_PRECISION_CORRECTION_LOOPS_ALLOWED 2

#define SEARCH_SWITCH_MIN_PRECISION_BINARY 1
#define SEARCH_SWITCH_MIN_PRECISION_EXHAUSTIVE 2
#define SEARCH_SWITCH_MIN_RHF_RUNTIME_EXHAUSTIVE 3

#define PRECISION_STEP_DOWN 5
#define PRECISION_STEP_UP 5
#define PRECISION_JUMP_FACTOR 2

#define ETA 0.51

// -----------------------------------------------------------------------------
// Errors
// Convert these into an enum
#define INDEX_SEARCH_FAIL 0
// Note that the reduction logic depends on INDEX_SEARCH_FAIL being set to 0
#define PRECISION_FAIL 1
#define ERROR_COMMANDLINE_ARGUMENT_NUMBER 2
#define ERROR_ALGORITHM_NUMBER 3
#define ERROR_FILE_READ 4
#define ERROR_INVALID_INTEGER_TYPE_IN_BASIS 5
#define ERROR_INVALID_PRECISION_RANGE 6
// -----------------------------------------------------------------------------

#define TYPE_MPZ_T 1
#define TYPE_LONG_LONG_INT 2
#define TYPE_LONG_INT 3
#define TYPE_INT 4

#define SIZE_LONG_LONG_INT ((4 * sizeof(long long int)) - 1)
#define SIZE_LONG_INT ((4 * sizeof(long int)) - 1)
#define SIZE_INT ((4 * sizeof(int)) - 1)

#define MAX_LONG_LONG_INT ((long long int) ((1 << SIZE_LONG_LONG_INT) - 1))
#define MAX_LONG_INT ((long int) ((1 << SIZE_LONG_INT) - 1))
#define MAX_INT ((int) ((1 << SIZE_INT) - 1))

template<typename Z_T>
inline bool basesEqual (ZZ_mat<Z_T> A, ZZ_mat<Z_T> A_previous, int d) 
{
    for (int i=0; i<d; i++)
    {
        for (int j=0; j<d; j++)
        {
            if (A[i][j] != A_previous[i][j])
            {
                return false;
            }
        }
    }
    return true;
}

template<typename Z_T>
inline void printBasis (Z_T ** &ppA, int d) 
{
    cout << "["; 
    for (int i=0; i<d; i++)
    {
        cout << "[";
        for (int j=0; j<d; j++)
        {
            cout << ppA[i][j] << " " ;
        }
        cout << "]" << endl;
    }
    cout << "]" << endl;

    return;
}

template <typename Z_T> 
inline void size_reduce (ZZ_mat<Z_T> &A, int d, int i, int j, Z_NR<Z_T> &nint_mu)
{
    nint_mu.neg(nint_mu); 
    for (int l=0; l<d; l++)
    {
        A[i][l].addmul(nint_mu, A[j][l]);
    }
    return;
}

template <typename T>
inline void size_reduce (T ** &ppA, int d, int i, int j, long int nint_mu)
{
    for (int l=0; l<d; l++)
    {
        ppA[i][l] -= ((T) nint_mu * ppA[j][l]);
    }
    return;
}

template <typename Z_T> 
inline void size_reduce_wrapper (
                ZZ_mat<Z_T> &A,
                FP_mat<mpfr_t> &mu,
                int d,
                int i,
                int j,
                mpfr_t &nint_mu_f, 
                FP_NR<mpfr_t> &nint_mu_FP,
                mpz_t &nint_mu_z,
                Z_NR<mpz_t> &max_nint_mu,
                Z_NR<mpz_t> &nint_mu,
                int start_index = 0
                )

{        
    mpfr_round(nint_mu_f, mu[i][j].get_data());
    mpfr_get_z (nint_mu_z, nint_mu_f, MPFR_RNDN);
    nint_mu = nint_mu_z;
    max_nint_mu = (nint_mu > max_nint_mu) ? nint_mu : max_nint_mu;

    size_reduce (A, d, i, j, nint_mu);

    // Update mu_ij
    nint_mu_FP = nint_mu_f;
    for (int l=j-1; l>=start_index; l--) 
    {
        mu[i][l] = mu[i][l] - nint_mu_FP * mu[j][l];
    }
    mu[i][j] = mu[i][j] - nint_mu_FP;

    return;
}

template <typename T> 
inline void size_reduce_wrapper (
                T ** &ppA,
                FP_mat<mpfr_t> &mu,
                int d,
                int i,
                int j,
                mpfr_t &nint_mu_f, 
                FP_NR<mpfr_t> &nint_mu_FP,
                long int &max_nint_mu,
                long int &nint_mu,
                int start_index = 0
                )
{
    mpfr_round(nint_mu_f, mu[i][j].get_data());
    nint_mu = mpfr_get_si (nint_mu_f, MPFR_RNDN);
    max_nint_mu = (nint_mu > max_nint_mu) ? nint_mu : max_nint_mu;

    size_reduce (ppA, d, i, j, nint_mu);
    nint_mu_FP = nint_mu_f;
            
    // Update mu_ij
    for (int l=j-1; l>=start_index; l--) 
    {
        mu[i][l] = mu[i][l] - nint_mu_FP * mu[j][l];
    }
    mu[i][j] = mu[i][j] - nint_mu_FP;

    return;
}


inline int getIntegerType (ZZ_mat<mpz_t> A, int d)
{
    int int_type = TYPE_MPZ_T;
    Z_NR<mpz_t> max_element;
    max_element.abs(A[0][0]);
    Z_NR<mpz_t> mpz_temp;

    if (d > 450)
        return int_type;

    for (int i=0; i<d; i++)
    {
        for (int j=0; j<d; j++)
        {
            mpz_temp.abs(A[i][j]);
            if (mpz_temp > max_element)
            {
                max_element = mpz_temp; 
            }
        }
    }

    if (max_element <= MAX_LONG_LONG_INT)
    {
        int_type = TYPE_LONG_LONG_INT;
    }

    if (max_element <= MAX_LONG_INT)
    {
        int_type = TYPE_LONG_INT;
    }

    if (max_element <= MAX_INT)
    {
        int_type = TYPE_INT;
    }
    cout << "int_type = " << int_type << endl;
    return int_type;
}

inline void copy_basis_mpz_mpfr (ZZ_mat<mpz_t> &A, FP_mat<mpfr_t> &A_mpfr, FP_mat<mpfr_t> &gso, FP_mat<mpfr_t> &mu, int d) 
{
    for (int i=0; i<d; i++)
    {
        for (int j=0; j<d; j++) 
        {
            A_mpfr(i,j).set_z(A(i,j));
            gso(i,j).set_z(A(i,j));
            mu[i][j] = 0.0;
        }
    }
    return; 
}

template <typename Z_T>
inline void copy_basis_Z_T_mpfr (Z_T ** &ppA, FP_mat<mpfr_t> &A_mpfr, FP_mat<mpfr_t> &gso, FP_mat<mpfr_t> &mu, int d) 
{
    Z_NR<mpz_t> temp;
    for (int i=0; i<d; i++)
    {
        for (int j=0; j<d; j++) 
        {
            temp = ppA[i][j];
            A_mpfr(i,j).set_z(temp);
            gso(i,j).set_z(temp);
            mu[i][j] = 0.0;
        }
    }
    return; 
}

template<typename T>
inline void copy_basis (T &A_destination, T &A_source, int d) 
{
    for (int i=0; i<d; i++)
    {
        for (int j=0; j<d; j++) 
        {
            A_destination[i][j] = A_source[i][j];
        }
    }
    return; 
}

template <typename Z_T>
void copy_basis_mpz_Z_T (Z_T ** &ppA, ZZ_mat<mpz_t> &A_source, int d) 
{
    long temp;
    ppA = (Z_T **) calloc (d, sizeof(Z_T *));
    for (int i=0; i<d; i++)
    {
        ppA[i] = (Z_T *) calloc (d, sizeof(Z_T));
    }
    for (int i=0; i<d; i++)
    {
        for (int j=0; j<d; j++) 
        {
            temp = A_source[i][j].get_si();
            ppA[i][j] = (Z_T) temp;
        }
    }
    return; 
}

// Computing the initial GSO, mu matrices, and B vector
template <typename FP_T>
inline void compute_GSO (
                ZZ_mat<mpz_t> &A,
                FP_mat<FP_T> &gso,
                FP_mat<FP_T> &mu,
                vector<FP_NR<FP_T>> &B ,
                int d,
                int start_index = 0,
                int end_index = 0
                )
{
    end_index = (!end_index) ? d-1 : end_index;

    // Copy mpz_t basis to mpfr_t basis A_mpfr, GSO
    // Also set each mu[i][j]=0
    FP_mat<FP_T> A_FP(d,d);
    copy_basis_mpz_mpfr (A, A_FP, gso, mu, d);

    FP_NR<FP_T> mu_num, mu_neg;
    for (int i=start_index; i<=end_index; i++) 
    {
        for (int j=i-1; j>=start_index; j--) 
        {
            // Compute numerator and denominator of mu
            A_FP[i].dot_product (mu_num, gso[j]);

            // Compute and store mu_ij in mu matrix
            mu[i][j] = mu_num / B[j];         

            // Update GSO as gso[i] - mu*gso[j]
            mu_neg = -mu[i][j];
            for (int l=0; l<d; l++) 
            {
                gso[i][l].addmul(mu_neg, gso[j][l]);
            }
        }
        // Once GSO updated for all j, store <b_j*,b_j*> in B
        // Then use B[j] instead of denominator computed above? 
        gso[i].dot_product (B[i], gso[i]);
    }
    return;
}

// Computing the initial GSO, mu matrices, and B vector
template <typename Z_T, typename FP_T>
inline void compute_GSO (
                Z_T ** &ppA,
                FP_mat<FP_T> &gso,
                FP_mat<FP_T> &mu,
                vector<FP_NR<FP_T>> &B ,
                int d,
                int start_index = 0,
                int end_index = 0
                )
{
    end_index = (!end_index) ? d-1 : end_index;

    // Copy mpz_t basis to mpfr_t bases A_mpfr and gso
    // Set all mu[i][j] = 0.0
    FP_mat<FP_T> A_FP(d,d);
    copy_basis_Z_T_mpfr (ppA, A_FP, gso, mu, d);

    FP_NR<FP_T> mu_num, mu_neg;
    // vector<FP_NR<FP_T>> fpVtmp;
    // fpVtmp.resize(d,0.0);
    for (int i=start_index; i<=end_index; i++) 
    {
        for (int j=i-1; j>=start_index; j--) 
        {
            // Compute numerator and denominator of mu
            A_FP[i].dot_product (mu_num, gso[j]);

            // Compute and store mu_ij in mu matrix
            mu[i][j] = mu_num / B[j];         

            // Update GSO as gso[i] - mu*gso[j]
            mu_neg = -mu[i][j];
            for (int l=0; l<d; l++) 
            {
                gso[i][l].addmul(mu_neg, gso[j][l]);
            }
        }
        // Once GSO updated for all j, store <b_j*,b_j*> in B
        // Then use B[j] instead of denominator computed above? 
        gso[i].dot_product (B[i], gso[i]);
    }
    return;
}

// Using the algorithm from YY18 for GSO update
template <typename FP_T>
inline int deep_LLL_GSO_update (
                FP_mat<FP_T> &mu,
                vector<FP_NR<FP_T>> &B,
                int k,
                int i,
                int d,
                int start_index = 0,
                int end_index = 0
                ) 
{
    end_index = (!end_index) ? d-1 : end_index;

    // Initialise and set P, S, D vectors
    vector<FP_NR<FP_T>> P, S, D;
    P.resize (d, 0.0);
    S.resize (d, 0.0);
    D.resize (d, 0.0);

    // Initialise and set constants
    FP_NR<FP_T> T = 0.0;
    FP_NR<FP_T> fptmp = 0.0;

    P[k] = B[k];
    D[k] = B[k];

    for (int j=k-1; j>=i; j--) 
    {
	// P[j] = r[k][j]
        P[j] = mu[k][j] * B[j];
        D[j] = D[j+1] + (mu[k][j] * P[j]);
    }

    for (int j=k; j>i; j--) 
    {
        T = mu[k][j-1] / D[j];

        for (int l=end_index; l>k; l--)
        {
            S[l] = S[l] + (mu[l][j] * P[j]);
            mu[l][j] = mu[l][j-1] - (T * S[l]);            
        }

        for (int l=k; l>=j+2; l--)
        {
            S[l] = S[l] + (mu[l-1][j] * P[j]);
            mu[l][j] = mu[l-1][j-1] - (T * S[l]);    
        }

        if (j != k)
        {
            S[j+1] = P[j];
            mu[j+1][j] = mu[j][j-1] - (T * S[j+1]);
        }
        
    }

    T = 1.0 / D[i];

    for (int l=end_index; l>k; l--)
    {
        mu[l][i] = T * (S[l] + mu[l][i] * P[i]);
    }

    for (int l=k; l>=i+2; l--) 
    {
        mu[l][i] = T * (S[l] + mu[l-1][i] * P[i]);
    }

    mu[i+1][i] = T * P[i];

    for (int j=start_index; j<i; j++) 
    {
        fptmp = mu[k][j];

        for (int l=k; l>i; l--)
        {
            mu[l][j] = mu[l-1][j];
        }

        mu[i][j] = fptmp;
    }

    for (int j=k; j>i; j--)
    {
        B[j] = (D[j] * B[j-1]) / D[j-1];
    }

    B[i] = D[i];

    // Free P,S,D vectors
    
    return 0;
}

/*
// Using the algorithm from YY18 for GSO update
// Using r[i][j] as well as mu[i][j]
template <typename FP_T>
inline int deep_LLL_GSO_update_COPY (
                FP_mat<FP_T> &mu,
                FP_mat<FP_T> &r,
                vector<FP_NR<FP_T>> &B,
                int k,
                int i,
                int d,
                int start_index = 0,
                int end_index = 0
                ) 
{
    end_index = (!end_index) ? d-1 : end_index;

    // Initialise and set P, S, D vectors
    vector<FP_NR<FP_T>> P, S, D;
    P.resize (d, 0.0);
    S.resize (d, 0.0);
    D.resize (d, 0.0);

    // Initialise and set constants
    FP_NR<FP_T> T = 0.0;
    FP_NR<FP_T> fptmp = 0.0;

    P[k] = B[k];
    D[k] = B[k];

    for (int j=k-1; j>=i; j--) 
    {
	P[j] = r[k][j];
        // P[j] = mu[k][j] * B[j];
        D[j] = D[j+1] + (mu[k][j] * P[j]);
    }

    for (int j=k; j>i; j--) 
    {
        T = mu[k][j-1] / D[j];

        for (int l=end_index; l>k; l--)
        {
            S[l] = S[l] + (mu[l][j] * P[j]);
            mu[l][j] = mu[l][j-1] - (T * S[l]);            
	    //r[l][j] = (mu[l][j-1] - (T * S[l])) * B[j];
        }

        for (int l=k; l>=j+2; l--)
        {
            S[l] = S[l] + (mu[l-1][j] * P[j]);
            mu[l][j] = mu[l-1][j-1] - (T * S[l]);    
	    // r[l][j] = ...
        }

        if (j != k)
        {
            S[j+1] = P[j];
            mu[j+1][j] = mu[j][j-1] - (T * S[j+1]);
	    // r[j+1][j] = ...
        }
        
    }

    T = 1.0 / D[i];

    for (int l=end_index; l>k; l--)
    {
        mu[l][i] = T * (S[l] + mu[l][i] * P[i]);
        // r[l][i] = ...
    }

    for (int l=k; l>=i+2; l--) 
    {
        mu[l][i] = T * (S[l] + mu[l-1][i] * P[i]);
        // r[l][i] = ...
    }

    mu[i+1][i] = T * P[i];
    // r[i+1][i] = ...

    for (int j=start_index; j<i; j++) 
    {
        fptmp = mu[k][j];

        for (int l=k; l>i; l--)
        {
            mu[l][j] = mu[l-1][j];
            // r[l][j] = ...
        }

        mu[i][j] = fptmp;
        // r[i][j] = ...
    }

    for (int j=k; j>i; j--)
    {
        B[j] = (D[j] * B[j-1]) / D[j-1];
    }

    B[i] = D[i];

    // Free P,S,D vectors
    
    return 0;
}
*/

// Compute the Root Hermite Factor of a Lattice Basis
template <typename FP_T>
FP_NR<FP_T> compute_RHF (
                vector<FP_NR<FP_T>> &B, 
                int d,
                int start_index = 0,
                int end_index = 0
                )
{
    end_index = (!end_index) ? d-1 : end_index;
   
    // Compute volume of lattice as product of GSO lengths
    // Do this by taking inner product of GSO vector with itself and square root
    FP_NR<FP_T> volumeSquared = 1.0;
    FP_NR<FP_T> volume;
   
    for (int i=start_index; i<=end_index; i++) {
       //RR_vol = RR_vol * gsoVectorNorm;
       volumeSquared *= B[i];
    }

    volume.sqrt(volumeSquared);

    // Compute the dth root of the volume using root function (from fplll)
    FP_NR<FP_T> dthRootVolume;
    dthRootVolume.root(volume, (end_index - start_index + 1));

    // Compute the norm of b1 using the square root of inner product
    FP_NR<FP_T> b1Norm; 
    b1Norm.sqrt(B[start_index]);

    // Compute the Hermite Factor as gamma = ||b_1|| / Vol(L)^(1/n)
    //RR hermiteFactor;
    FP_NR<FP_T> hermiteFactor;
    hermiteFactor = (b1Norm / dthRootVolume);

    // Compute the Root Hermite Factor (RHF) as gamma^(1/n)
    FP_NR<FP_T> RHF;
    RHF.root(hermiteFactor, (end_index - start_index + 1));

    return RHF;
}

#endif
