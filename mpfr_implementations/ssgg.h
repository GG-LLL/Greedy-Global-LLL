/* 
 * *****************************************************************************
 * Filename: ssgg.h
 * Authors : Sanjay Bhattacherjee and Jack Moyler
 * *****************************************************************************
*/

#include "common.h"

#ifndef SSGG
#define SSGG

#define SS_GG_DELTA 0.999999

template <typename FP_T>
inline FP_NR<FP_T> compute_SS (
                FP_NR<FP_T> &squared_sum,
                vector<FP_NR<FP_T>> &B,
                int d,
                int start_index = 0,
                int end_index = 0
                )
{
    end_index = (!end_index) ? d-1 : end_index;

    squared_sum = 0.0; 
    for (int i=start_index; i<=end_index; i++) 
    {
        squared_sum += B[i];
    } 
    return squared_sum;
}

template <typename FP_T>
inline bool index_search_SS (
                FP_mat<FP_T> &mu,
                vector<FP_NR<FP_T>> &B,
                int d,
                FP_NR<FP_T> squared_sum,
                FP_NR<FP_T> delta_p,
                int &k_max,
                int &i_max,
                int start_index = 0,
                int end_index = 0
                ) 
{
    end_index = (!end_index) ? d-1 : end_index;

    k_max = INDEX_SEARCH_FAIL;
    i_max = start_index; 

    FP_NR<FP_T> SS_change = 0.0;
    FP_NR<FP_T> SS_change_max = delta_p * squared_sum;
    FP_NR<FP_T> projection_l_k = 0.0;
    FP_NR<FP_T> ratio_projection_l_k = 0.0;
    FP_NR<FP_T> inv_ratio_projection_l_k = 0.0;

    for (int k=start_index+1; k<=end_index; k++) 
    {
        // projection_l_k = \pi_{k-1}(\vb_{k}) (initially)
        // D^(k)_{k-1} = (projection_l_k)^2
        projection_l_k = B[k] + (mu[k][k-1] * mu[k][k-1] * B[k-1]);
        // projection_l_l = \pi_{l}(\vb_{l}) = pB[k-1]
        // ratio_projection_l_k = projection_l_k / projection_l_l
        ratio_projection_l_k = projection_l_k / B[k-1];
        // SS_change = Old_SS - New_SS
        SS_change = mu[k][k-1] * mu[k][k-1] * (1.0 - ratio_projection_l_k);  

        if (SS_change > SS_change_max) 
        {
            SS_change_max = SS_change;
            k_max = k;
            i_max = k-1;
        }

        for (int i = k-2; i >= start_index; i--) 
        {
            // projection_l_k = \pi_{i}(\vb_{k}) 
            // D^(k)_{i} = (projection_l_k)^2
            projection_l_k += (mu[k][i] * mu[k][i] * B[i]);
            // projection_l_l = \pi_{l}(\vb_{l}) = pB[i]
            // ratio_projection_l_k = projection_l_k / projection_l_l
            inv_ratio_projection_l_k = B[i] / projection_l_k;
            // SS_change = Old_SS - New_SS
            SS_change += mu[k][i] * mu[k][i] * B[i] * (inv_ratio_projection_l_k - 1.0);  
            
            if (SS_change > SS_change_max) 
            {
                SS_change_max = SS_change;
                k_max = k;
                i_max = i;
            }
        }
    }
    
    if (k_max != INDEX_SEARCH_FAIL)
    {
        return true;
    } 
    return false;
}

template <typename FP_T>
void ssgg_reduced_check_internal (
                FP_mat<FP_T> &mu,
                vector<FP_NR<FP_T>> &B,
                FP_NR<FP_T> delta_p,
                int precision,
                int d,
                FP_NR<FP_T> eta,
                int &k,
                int &i,
                FP_NR<FP_T> &squared_sum,
                bool &size_reduction_required,
                bool &deep_insertion_required,
                bool &index_search_required,
                int start_index = 0,
                int end_index = 0
                )
{
    end_index = (!end_index) ? d-1 : end_index;

    size_reduction_required = true;
    deep_insertion_required = true;
    index_search_required = true;
    
    for (k=start_index+1; k<=end_index; k++) 
    {
        for (int j=k-1; j>=start_index; j--)
        {
            if (mu[k][j] < -eta || mu[k][j] > eta)
            {
                cout << "ERROR: BASIS NOT SIZE REDUCED" << endl;
                cout << "mu[k][j] = " << mu[k][j] << endl;
                return;
            }
        }
    }
    
    size_reduction_required = false;
    k = INDEX_SEARCH_FAIL;
    i = -1;

    squared_sum = compute_SS (squared_sum, B, d, start_index, end_index);
    deep_insertion_required = index_search_SS (mu, B, d, squared_sum, delta_p, k, i, start_index, end_index);
    index_search_required = false;

    if (deep_insertion_required) 
    {
        cout << "ERROR: BASIS SIZE REDUCED BUT NOT SS REDUCED" << endl;
    }

    return;
}

// SS-GG rediction for ZZ_mat<mpz_t> type basis
int ssgg (
                ZZ_mat<mpz_t> &A,
                FP_NR<mpfr_t> delta,
                int d,
                int precision,
                int precision_correction_loops_allowed,
                FP_NR<mpfr_t> eta,
                FP_mat<mpfr_t> &gso,
                FP_mat<mpfr_t> &mu,
                vector<FP_NR<mpfr_t>> &B,
                int start_index = 0,
                int end_index = 0
                )
{
    end_index = (!end_index) ? d-1 : end_index;

    // FP_NR<mpfr_t>::set_prec(precision);

    FP_NR<mpfr_t> delta_p;
    delta_p = 1.0 - delta;
    
    FP_NR<mpfr_t> squared_sum;
    FP_NR<mpfr_t> previous_squared_sum;
    
    mpz_t nint_mu_z;
    mpz_init (nint_mu_z);
    mpfr_t nint_mu_f;
    mpfr_init2 (nint_mu_f, mpfr_prec_t (precision));
    FP_NR<mpfr_t> nint_mu_FP = 0.0;
    Z_NR<mpz_t> max_nint_mu, nint_mu, MAX_nint_mu;
    nint_mu = 0;
    max_nint_mu = 0; 
    MAX_nint_mu = MAX_INT; 

    // GSO - initial
    compute_GSO (A, gso, mu, B, d, start_index, end_index);

    bool size_reduction_required = true;
    bool deep_insertion_required = true;
    bool index_search_required = true;

    int k = start_index+1;
    int i = start_index;
    int nInsertions = 0;
    int nLoops = 1;

    while (size_reduction_required || deep_insertion_required)
    {
        if (nLoops > precision_correction_loops_allowed)
        {
            cout << "nLoops " << nLoops << endl;
            cout << "Precision " << precision << " insufficient." << endl;
            return PRECISION_FAIL;
        }

        if (size_reduction_required)
        {
            // Size reduction - initial
            for (int l=k; l<=end_index; l++)
            {
                for (int j=l-1; j>=start_index; j--)
                {
                    // Takes mu[i][j]
                    // Computes nearest mpz_t to mu_[i][j]
                    if (mu[l][j] < -eta || mu[l][j] > eta)
                    {
                        size_reduce_wrapper (A, mu, d, l, j, 
                        nint_mu_f, 
                        nint_mu_FP,
                        nint_mu_z,
                        max_nint_mu,
                        nint_mu,
                        start_index
                        );
                    }
                }
            }

            size_reduction_required = false;
            deep_insertion_required = true;
            index_search_required = true;
        }

        if (index_search_required)
        {
            // Computing the squared sum on input
            squared_sum = compute_SS (squared_sum, B, d, start_index, end_index);
            previous_squared_sum = squared_sum;

            // Index search - initial
            deep_insertion_required = index_search_SS (mu, B, d, squared_sum, delta_p, k, i, start_index, end_index);
        }
        cout << "SS on input = " << squared_sum << endl;

        // Loop
        while (deep_insertion_required) 
        {            
            nInsertions++;
            // Deep insert b_k in position i
            // rotate_right function from fplll
            // (v[first],...,v[last]) becomes (v[last],v[first],...,v[last-1])
            A.rotate_right(i, k);

            // Update the GSO     
            deep_LLL_GSO_update (mu, B, k, i, d, start_index, end_index); 

            // Size reduce from i to end_index
            for (int l=i+1; l<=end_index; l++)
            {
                for (int j=l-1; j>=start_index; j--)
                {
                    // Takes mu[i][j]
                    // Computes nearest mpz_t to mu_[i][j]
                    if (mu[l][j] < -eta || mu[l][j] > eta)
                    {
                        size_reduce_wrapper(A, mu, d, l, j, 
                        nint_mu_f, 
                        nint_mu_FP,
                        nint_mu_z,
                        max_nint_mu,
                        nint_mu,
                        start_index
                        );
                    }
                }
            }

            // Index search
            squared_sum = compute_SS (squared_sum, B, d, start_index, end_index);
            deep_insertion_required = index_search_SS (mu, B, d, squared_sum, delta_p, k, i, start_index, end_index);
        }

        // Check the SS_GG_Reducedness
        compute_GSO (A, gso, mu, B, d, start_index, end_index);

        ssgg_reduced_check_internal(
                        mu,
                        B,
                        delta_p,
                        precision,
                        d,
                        eta,
                        k,
                        i,
                        squared_sum,
                        size_reduction_required,
                        deep_insertion_required,
                        index_search_required,
                        start_index,
                        end_index
                        );

        if (squared_sum >= previous_squared_sum)
        {
            cout << "nLoops " << nLoops << endl;
            cout << "previous_squared_sum " << previous_squared_sum << endl;
            cout << "squared_sum " << squared_sum << endl;
            cout << "Precision " << precision << " insufficient." << endl;
            return PRECISION_FAIL;
        }

        nLoops++;
    }
    cout << "SS on output = " << squared_sum << endl;

    mpfr_clear(nint_mu_f);
    cout << "nInsertions = " << nInsertions << endl;
    Z_NR<mpz_t> nint_mu_difference;
    nint_mu_difference.sub(MAX_nint_mu, max_nint_mu);
    cout << "MAX_nint_mu - max_nint_mu = " << nint_mu_difference << endl;

    return RED_SUCCESS;
}

// SS-GG rediction for Z_T ** type basis
template <typename Z_T>
int ssgg (
                Z_T ** &ppA,
                FP_NR<mpfr_t> delta,
                int d,
                int precision,
                int precision_correction_loops_allowed,
                FP_NR<mpfr_t> eta,
                FP_mat<mpfr_t> &gso,
                FP_mat<mpfr_t> &mu,
                vector<FP_NR<mpfr_t>> &B,
                int start_index = 0,
                int end_index = 0
                )
{
    end_index = (!end_index) ? d-1 : end_index;

    // FP_NR<mpfr_t>::set_prec(precision);

    FP_NR<mpfr_t> delta_p;
    delta_p = 1.0 - delta;
    
    FP_NR<mpfr_t> squared_sum;
    FP_NR<mpfr_t> previous_squared_sum;
    
    FP_NR<mpfr_t> nint_mu_FP = 0.0;
    mpfr_t nint_mu_f;
    mpfr_init2(nint_mu_f, mpfr_prec_t (precision));
    long int max_nint_mu = 0;
    long int nint_mu = 0;

    Z_T * pTemp;

    // GSO - initial
    compute_GSO (ppA, gso, mu, B, d, start_index, end_index);

    bool size_reduction_required = true; 
    bool deep_insertion_required = true; 
    bool index_search_required = true; 

    int k = start_index+1;
    int i = start_index;
    int nInsertions = 0;
    int nLoops = 1;

    while (size_reduction_required || deep_insertion_required)
    {
        if (nLoops > precision_correction_loops_allowed)
        {
            cout << "nLoops " << nLoops << endl;
            cout << "Precision " << precision << " insufficient." << endl;
            return PRECISION_FAIL;
        }

        if (size_reduction_required)
        {
            // Size reduction - initial
            for (int l=k; l<=end_index; l++)
            {
                for (int j=l-1; j>=start_index; j--)
                {
                    // Takes mu[i][j]
                    // Computes nearest mpz_t to mu_[i][j]
                    if (mu[l][j] < -eta || mu[l][j] > eta)
                    {
                        size_reduce_wrapper(ppA, mu, d, l, j, 
                        nint_mu_f, 
                        nint_mu_FP,
                        max_nint_mu,
                        nint_mu,
                        start_index
                        );
                    }
                }
            }

            size_reduction_required = false;
            deep_insertion_required = true;
            index_search_required = true;
        }

        if (index_search_required)
        {
            // Computing the squared sum on input
            squared_sum = compute_SS (squared_sum, B, d, start_index, end_index);
            previous_squared_sum = squared_sum;

            // Index search - initial
            deep_insertion_required = index_search_SS (mu, B, d, squared_sum, delta_p, k, i, start_index, end_index);
        }

        cout << "SS on input = " << squared_sum << endl;

        // Loop
        while (deep_insertion_required) 
        {            
            nInsertions++;
            // Deep insert b_k in position i
            // (v[first],...,v[last]) becomes (v[last],v[first],...,v[last-1])
            pTemp = ppA[k];
            for (int j=k-1; j>=i; j--)
            {
                ppA[j+1] = ppA[j];
            }
            ppA[i] = pTemp;

            // Update the GSO     
            deep_LLL_GSO_update (mu, B, k, i, d, start_index, end_index); 

            // Size reduce from i to n
            for (int l=i+1; l<=end_index; l++)
            {
                for (int j=l-1; j>=start_index; j--)
                {
                    // Takes mu[i][j]
                    // Computes nearest mpz_t to mu_[i][j]
                    if (mu[l][j] < -eta || mu[l][j] > eta)
                    {
                        size_reduce_wrapper(ppA, mu, d, l, j, 
                        nint_mu_f, 
                        nint_mu_FP,
                        max_nint_mu,
                        nint_mu,
                        start_index
                        );
                    }
                }
            }

            // Index search
            squared_sum = compute_SS (squared_sum, B, d, start_index, end_index);
            deep_insertion_required = index_search_SS (mu, B, d, squared_sum, delta_p, k, i, start_index, end_index);
        }

        // Check the SS_GG_Reducedness
        compute_GSO (ppA, gso, mu, B, d, start_index, end_index);

        ssgg_reduced_check_internal (
                        mu,
                        B,
                        delta_p,
                        precision,
                        d,
                        eta,
                        k,
                        i,
                        squared_sum,
                        size_reduction_required,
                        deep_insertion_required,
                        index_search_required,
                        start_index,
                        end_index
                        );

        if (squared_sum >= previous_squared_sum)
        {
            cout << "nLoops " << nLoops << endl;
            cout << "Precision " << precision << " insufficient." << endl;
            return PRECISION_FAIL;
        }

        nLoops++;
    }
    cout << "SS on output = " << squared_sum << endl;

    mpfr_clear(nint_mu_f);
    cout << "nInsertions = " << nInsertions << endl;
    cout << "max_nint_mu = " << max_nint_mu << endl;
    cout << "MAX_INT - max_nint_mu = " << MAX_INT - max_nint_mu << endl;
    cout << "MAX_LONG_INT - max_nint_mu = " << MAX_LONG_INT - max_nint_mu << endl;

    return RED_SUCCESS;
}

int ssgg_wrapper (
                ZZ_mat<mpz_t> &A,
                FP_NR<mpfr_t> delta_SS,
                int d,
                int precision,
                int precision_correction_loops_allowed,
                FP_NR<mpfr_t> eta,
                long double &time_taken,
                FP_NR<mpfr_t> &RHF,
                int start_index = 0,
                int end_index = 0,
                int sublattice_size = DEFAULT_SUBLATTICE_SIZE
                )
{
    end_index = (!end_index) ? d-1 : end_index;
    sublattice_size = (sublattice_size <= 0) ? d : sublattice_size;

    FP_NR<mpfr_t>::set_prec(precision);
    // gso - The GSO basis
    // mu - The matrix containing all mu_{i,j} for 0 <= j < i < d
    // B - Vector of ||b_i*||^2 for 0 <= i < d
    FP_mat<mpfr_t> gso(d,d);
    FP_mat<mpfr_t> mu(d,d);
    vector<FP_NR<mpfr_t>> B;
    B.resize(d,0.0);

    int prec_status; 

    int ** ppIntA = NULL;
    int ** ppLongIntA = NULL;
    int ** ppLongLongIntA = NULL;

    clock_t begin_time;
    clock_t end_time;

    int int_type=0;
    int int_type_new=0;

    // *************************************************************************
    begin_time = clock();

    int_type = getIntegerType(A, d);
    switch (int_type)
    {
        case TYPE_MPZ_T:
            break;
        case TYPE_LONG_LONG_INT:
            copy_basis_mpz_Z_T (ppLongLongIntA, A, d); 
            break;
        case TYPE_LONG_INT:
            copy_basis_mpz_Z_T (ppLongIntA, A, d); 
            break;
        case TYPE_INT:
            copy_basis_mpz_Z_T (ppIntA, A, d); 
            break;
        default:
            cout << "ERROR: BASIS INTEGER DATATYPE NOT SUPPORTED" << endl;
            exit(ERROR_INVALID_INTEGER_TYPE_IN_BASIS);
    };

    int sublattice_start, sublattice_end;
    sublattice_start = start_index;
    sublattice_end = 0;
    while (sublattice_start < end_index)
    {
        sublattice_end += sublattice_size;
        sublattice_end = (sublattice_end > end_index) ? end_index : sublattice_end;
        cout << endl;
        // cout << "sublattice_start = " << sublattice_start << endl;
        // cout << "sublattice_end = " << sublattice_end << endl;

        switch (int_type)
        {
            case TYPE_MPZ_T:
                prec_status = ssgg (A, delta_SS, d, precision, precision_correction_loops_allowed, eta, gso, mu, B, sublattice_start, sublattice_end);
                break;
            case TYPE_LONG_LONG_INT:
                prec_status = ssgg (ppLongLongIntA, delta_SS, d, precision, precision_correction_loops_allowed, eta, gso, mu, B, sublattice_start, sublattice_end);
                break;
            case TYPE_LONG_INT:
                prec_status = ssgg (ppLongIntA, delta_SS, d, precision, precision_correction_loops_allowed, eta, gso, mu, B, sublattice_start, sublattice_end);
                break;
            case TYPE_INT:
                prec_status = ssgg (ppIntA, delta_SS, d, precision, precision_correction_loops_allowed, eta, gso, mu, B, sublattice_start, sublattice_end);
                break;
            default:
                cout << "ERROR: BASIS INTEGER DATATYPE NOT SUPPORTED" << endl;
                exit(ERROR_INVALID_INTEGER_TYPE_IN_BASIS);
        };
        sublattice_start += sublattice_size;
    }

    int_type_new = getIntegerType(A, d);

    if (sublattice_size != d)
    {
        cout << "Preprocessing complete. Now running SSGG globally." << endl;
        switch (int_type)
        {
            case TYPE_MPZ_T:
                prec_status = ssgg (A, delta_SS, d, precision, precision_correction_loops_allowed, eta, gso, mu, B, start_index, end_index);
                break;
            case TYPE_LONG_LONG_INT:
                prec_status = ssgg (ppLongLongIntA, delta_SS, d, precision, precision_correction_loops_allowed, eta, gso, mu, B, start_index, end_index);
                break;
            case TYPE_LONG_INT:
                prec_status = ssgg (ppLongIntA, delta_SS, d, precision, precision_correction_loops_allowed, eta, gso, mu, B, start_index, end_index);
                break;
            case TYPE_INT:
                prec_status = ssgg (ppIntA, delta_SS, d, precision, precision_correction_loops_allowed, eta, gso, mu, B, start_index, end_index);
                break;
            default:
                cout << "ERROR: BASIS INTEGER DATATYPE NOT SUPPORTED" << endl;
                exit(ERROR_INVALID_INTEGER_TYPE_IN_BASIS);
        };
    }
    end_time = clock();
    // *************************************************************************

    time_taken = (long double) (end_time - begin_time) / CLOCKS_PER_SEC;
    cout << "Time=" << time_taken;
    RHF = compute_RHF (B, d, start_index, end_index);
    cout << "\tRHF=" << RHF << endl << endl;
            
    return prec_status;
}

#endif
