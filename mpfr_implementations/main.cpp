/* 
 * *****************************************************************************
 * Filename: main.cpp
 * Authors : Sanjay Bhattacherjee and Jack Moyler
 * *****************************************************************************
*/

#include "ssgg.h"
#include "potgg.h"

int main (int argc, char * argv[])
{
    if (argc<3 || argc>7) 
    {
        cout << "Usage: <executable> <input-file-name> <algorithm> [precision = 53] [search-switch = 0] [sublattice_size = 0] [precision-correction-loops-allowed = 2] [eta = 0.51]" << endl;
        cout << "Algorithm:" << endl;
        cout << " 1 :  SS-GGLLL" << endl;
        cout << " 2 :  Pot-GGLLL" << endl;
        cout << "precision:" << endl;
        cout << " the number of bits of significand" << endl;
        cout << "search-switch: " << endl;
        cout << " 0: no search, using only the specified precision" << endl;
        cout << " 1: (min precision) binary search" << endl;
        cout << " 2: (min precision) exhaustive search starting from the specified precision" << endl;
        cout << " 3: (precision for min RHF + min time) exhaustive search starting from the specified precision" << endl;
        cout << "sublattice_size:" << endl;
        cout << " non-zero if X-GG sublattice preprocessing is desired" << endl;
        cout << "precision-correction-loops-allowed: " << endl;
        cout << " k: GSO to be recomputed k-1 times" << endl;
        cout << "eta:" << endl;
        cout << " floating point size reduction relaxation" << endl;
        exit (ERROR_COMMANDLINE_ARGUMENT_NUMBER);
    }

    int precision;
    FP_NR<mpfr_t> RHF;
    long double time_taken;
    int sublattice_size;
    int search_switch;
    int precision_correction_loops_allowed;
    FP_NR<mpfr_t> eta, delta_SS, delta_Pot;

    int precision_min_RHF, precision_min_runtime;
    FP_NR<mpfr_t> RHF_min, RHF_min_runtime;
    long double runtime_min, runtime_min_RHF;

    int precision_fail, precision_mid, precision_success;
    int precision_jump;
    FP_NR<mpfr_t> prev_RHF;

    string input_filename (argv[1]);
    cout << "argv[1] input_filename: " << input_filename << endl;
    int algo_num = atoi(argv[2]); 
    cout << "argv[2] algo_num: " << algo_num << endl;
    precision =
            (argc>=4) ? atoi (argv[3]) : DEFAULT_PRECISION;
    cout << "argv[3] precision: " << precision << endl;
    search_switch = 
            (argc>=5) ? atoi (argv[4]) : DEFAULT_SEARCH_SWITCH;
    cout << "argv[4] search_switch: " << search_switch << endl;
    sublattice_size =
            (argc>=6) ? atoi (argv[5]) : DEFAULT_SUBLATTICE_SIZE;
    cout << "argv[5] sublattice_size: " << sublattice_size << endl;
    precision_correction_loops_allowed = 
            (argc>=7) ? atoi (argv[6]) : DEFAULT_PRECISION_CORRECTION_LOOPS_ALLOWED;
    cout << "argv[6] precision_correction_loops_allowed: " << precision_correction_loops_allowed << endl;
    eta = 
            (argc>=8) ? atoi (argv[7]) : ETA;
    cout << "argv[7] eta: " << eta << endl;
    delta_SS = SS_GG_DELTA;
    delta_Pot = POT_GG_DELTA;

    ZZ_mat<mpz_t> A;

    int status = 0;
    status |= read_file(A, input_filename.c_str());
    if (status) exit (ERROR_FILE_READ);

    int d = A.get_cols();

    int prec_status;

    if (algo_num == 1)
    {
        cout << endl;
        cout << endl << "Starting with precision = " << precision << endl;

        prec_status = ssgg_wrapper (A, delta_SS, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, 0, d-1, sublattice_size);

        cout << endl;

        if (search_switch == SEARCH_SWITCH_MIN_PRECISION_BINARY)
        {
            // Find an insufficient precison
            precision_fail = precision;
            precision_jump = PRECISION_STEP_DOWN;
            bool min_check = false;

            while (prec_status == RED_SUCCESS)
            {
                min_check = true;
		precision_success = precision_fail;
                cout << endl;
                cout << "Succeeded with precision_fail = " << precision_fail << endl;
                precision_fail -= precision_jump;
		precision_fail = (precision_fail <= 0) ? 1 : precision_fail;
                precision_jump *= PRECISION_JUMP_FACTOR;
                cout << "Checking with precision_fail = " << precision_fail << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
                prec_status = ssgg_wrapper (A, delta_SS, d, precision_fail, precision_correction_loops_allowed, eta, time_taken, RHF, 0, d-1, sublattice_size);
            }

            // Find a sufficient precison
            if (!min_check)
            {
                precision_jump = PRECISION_STEP_UP;
                precision_success = precision_fail + precision_jump;

                cout << endl;
                cout << "Checking with precision_success = " << precision_success << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
                prec_status = ssgg_wrapper (A, delta_SS, d, precision_success, precision_correction_loops_allowed, eta, time_taken, RHF, 0, d-1, sublattice_size);

                while (prec_status != RED_SUCCESS)
                {
                    cout << "Failed with precision_success = " << precision_success << endl;
                    precision_fail = precision_success;
                    precision_jump *= PRECISION_JUMP_FACTOR;
                    precision_success += precision_jump;
                    cout << "Checking with precision_success = " << precision_success << endl;

                    status |= read_file(A, input_filename.c_str());
                    if (status) exit (ERROR_FILE_READ);
                    prec_status = ssgg_wrapper (A, delta_SS, d, precision_success, precision_correction_loops_allowed, eta, time_taken, RHF, 0, d-1, sublattice_size);
                }
            }

            if (precision_fail > precision_success)
            {
                cout << "Error: In finding precision range; precision_fail = " << precision_fail << " precision_success = " << precision_success << endl;
                exit (ERROR_INVALID_PRECISION_RANGE);
            }

            // Conduct binary search for the minimum precision
            cout << endl << "Starting binary search with" << endl;
            cout << "Inital precision_fail = " << precision_fail << endl;
            cout << "Inital precision_success = " << precision_success << endl << endl;
            while (precision_fail+1 < precision_success)
            {
                precision_mid = precision_fail + precision_success;
                precision_mid = (precision_mid % 2) ? (precision_mid+1)/2 : precision_mid/2; 
                cout << endl << "precision_fail = " << precision_fail << endl;
                cout << "precision_mid = " << precision_mid << endl;
                cout << "precision_success = " << precision_success << endl << endl;
                cout << "Checking with precision_mid = " << precision_mid << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
                prec_status = 
                    ssgg_wrapper (A, delta_SS, d, precision_mid, precision_correction_loops_allowed, eta, time_taken, RHF, 0, d-1, sublattice_size);
                if (prec_status != RED_SUCCESS)
                    precision_fail = precision_mid;
                else
                    precision_success = precision_mid;
            }

            precision = precision_success;
            cout << "Running with Precision " << precision << endl;

            status |= read_file(A, input_filename.c_str());
            if (status) exit (ERROR_FILE_READ);
            prec_status = 
                ssgg_wrapper (A, delta_SS, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, 0, d-1, sublattice_size);
        }
        else if (search_switch == SEARCH_SWITCH_MIN_PRECISION_EXHAUSTIVE)
        {
            while (prec_status != RED_SUCCESS)
            {
                precision++;
                cout << "Running with Precision " << precision << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
                prec_status = ssgg_wrapper (A, delta_SS, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, 0, d-1, sublattice_size);
            }
        }
        else if (search_switch == SEARCH_SWITCH_MIN_RHF_RUNTIME_EXHAUSTIVE)
        {
            while (prec_status != RED_SUCCESS)
            {
                precision++;
                cout << "Running with Precision " << precision << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
                prec_status = ssgg_wrapper (A, delta_SS, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, 0, d-1, sublattice_size);
            }
            precision_min_RHF = precision;
            precision_min_runtime = precision;
            RHF_min = RHF;
            RHF_min_runtime = RHF;
            runtime_min = time_taken;
            runtime_min_RHF = time_taken;
            do
            {
                precision++;
                cout << "Running with Precision " << precision << endl;

                prev_RHF = RHF;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
                prec_status = ssgg_wrapper (A, delta_SS, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, 0, d-1, sublattice_size);

                if (prec_status == RED_SUCCESS)
                {
                    if (RHF < RHF_min)
                    {
                        RHF_min = RHF;
                        runtime_min_RHF = time_taken;
                        precision_min_RHF = precision;
                    }

                    if (time_taken < runtime_min)
                    {
                        RHF_min_runtime = RHF;
                        runtime_min = time_taken;
                        precision_min_runtime = precision;
                    }
                }
            }
            while (prev_RHF != RHF || prec_status != RED_SUCCESS);

            cout << "Min time:\ttime=" << runtime_min;
            cout << "\tRHF=" << RHF_min_runtime;
            cout << "\td=" << d;
            cout << "\talgo_num=" << algo_num;
            cout << "\tprec=" << precision_min_runtime;
            cout << "\tsublattice=" << sublattice_size;
            cout << "\tswitch=" << search_switch;
            cout << "\tloops=" << precision_correction_loops_allowed;
            cout << "\teta=" << eta << endl;

            cout << "Min RHF:\ttime=" << runtime_min_RHF;
            cout << "\tRHF=" << RHF_min;
            cout << "\td=" << d;
            cout << "\talgo_num=" << algo_num;
            cout << "\tprec=" << precision_min_RHF;
            cout << "\tsublattice=" << sublattice_size;
            cout << "\tswitch=" << search_switch;
            cout << "\tloops=" << precision_correction_loops_allowed;
            cout << "\teta=" << eta << endl;

        }
        else
        {
            // The default case has been handled already
        }
    } 
    else if (algo_num == 2) 
    {
        cout << endl;
        cout << endl << "Starting with precision = " << precision << endl;

        prec_status = potgg_wrapper (A, delta_Pot, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, 0, d-1, sublattice_size);

        cout << endl;

        if (search_switch == SEARCH_SWITCH_MIN_PRECISION_BINARY)
        {
            // Find an insufficient precison
            precision_fail = precision;
            precision_jump = PRECISION_STEP_DOWN;
            bool min_check = false;

            while (prec_status == RED_SUCCESS)
            {
                min_check = true;
		precision_success = precision_fail;
                cout << endl;
                cout << "Succeeded with precision_fail = " << precision_fail << endl;
                precision_fail -= precision_jump;
                precision_jump *= PRECISION_JUMP_FACTOR;
                cout << "Checking with precision_fail = " << precision_fail << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
                prec_status = potgg_wrapper (A, delta_Pot, d, precision_fail, precision_correction_loops_allowed, eta, time_taken, RHF, 0, d-1, sublattice_size);
            }

            // Find a sufficient precison
            if (!min_check)
            {
                precision_jump = PRECISION_STEP_UP;
                precision_success = precision_fail + precision_jump;

                cout << endl;
                cout << "Checking with precision_success = " << precision_success << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
                prec_status = potgg_wrapper (A, delta_Pot, d, precision_success, precision_correction_loops_allowed, eta, time_taken, RHF, 0, d-1, sublattice_size);

                while (prec_status != RED_SUCCESS)
                {
                    cout << "Failed with precision_success = " << precision_success << endl;
                    precision_fail = precision_success;
                    precision_jump *= PRECISION_JUMP_FACTOR;
                    precision_success += precision_jump;
                    cout << "Checking with precision_success = " << precision_success << endl;

                    status |= read_file(A, input_filename.c_str());
                    if (status) exit (ERROR_FILE_READ);
                    prec_status = potgg_wrapper (A, delta_Pot, d, precision_success, precision_correction_loops_allowed, eta, time_taken, RHF, 0, d-1, sublattice_size);
                }
            }

            if (precision_fail > precision_success)
            {
                cout << "Error: In finding precision range; precision_fail = " << precision_fail << " precision_success = " << precision_success << endl;
                exit (ERROR_INVALID_PRECISION_RANGE);
            }

            // Conduct binary search for the minimum precision
            cout << endl << "Starting binary search with" << endl;
            cout << "Inital precision_fail = " << precision_fail << endl;
            cout << "Inital precision_success = " << precision_success << endl << endl;
            while (precision_fail+1 < precision_success)
            {
                precision_mid = precision_fail + precision_success;
                precision_mid = (precision_mid % 2) ? (precision_mid+1)/2 : precision_mid/2; 
                cout << endl << "precision_fail = " << precision_fail << endl;
                cout << "precision_mid = " << precision_mid << endl;
                cout << "precision_success = " << precision_success << endl << endl;
                cout << "Checking with precision_mid = " << precision_mid << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
                prec_status = 
                    potgg_wrapper (A, delta_Pot, d, precision_mid, precision_correction_loops_allowed, eta, time_taken, RHF, 0, d-1, sublattice_size);
                if (prec_status != RED_SUCCESS)
                    precision_fail = precision_mid;
                else
                    precision_success = precision_mid;
            }

            precision = precision_success;
            cout << "Running with Precision " << precision << endl;

            status |= read_file(A, input_filename.c_str());
            if (status) exit (ERROR_FILE_READ);
            prec_status = 
                potgg_wrapper (A, delta_Pot, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, 0, d-1, sublattice_size);
        }
        else if (search_switch == SEARCH_SWITCH_MIN_PRECISION_EXHAUSTIVE)
        {
            while (prec_status != RED_SUCCESS)
            {
                precision++;
                cout << "Running with Precision " << precision << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
                prec_status = potgg_wrapper (A, delta_Pot, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, 0, d-1, sublattice_size);
            }
        }
        else if (search_switch == SEARCH_SWITCH_MIN_RHF_RUNTIME_EXHAUSTIVE)
        {
            while (prec_status != RED_SUCCESS)
            {
                precision++;
                cout << "Running with Precision " << precision << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
                prec_status = potgg_wrapper (A, delta_Pot, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, 0, d-1, sublattice_size);
            }
            precision_min_RHF = precision;
            precision_min_runtime = precision;
            RHF_min = RHF;
            RHF_min_runtime = RHF;
            runtime_min = time_taken;
            runtime_min_RHF = time_taken;
            do
            {
                precision++;
                cout << "Running with Precision " << precision << endl;

                prev_RHF = RHF;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
                prec_status = potgg_wrapper (A, delta_Pot, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, 0, d-1, sublattice_size);

                if (prec_status == RED_SUCCESS)
                {
                    if (RHF < RHF_min)
                    {
                        RHF_min = RHF;
                        runtime_min_RHF = time_taken;
                        precision_min_RHF = precision;
                    }

                    if (time_taken < runtime_min)
                    {
                        RHF_min_runtime = RHF;
                        runtime_min = time_taken;
                        precision_min_runtime = precision;
                    }
                }
            }
            while (prev_RHF != RHF || prec_status != RED_SUCCESS);

            cout << "Min time:\ttime=" << runtime_min;
            cout << "\tRHF=" << RHF_min_runtime;
            cout << "\td=" << d;
            cout << "\talgo_num=" << algo_num;
            cout << "\tprec=" << precision_min_runtime;
            cout << "\tsublattice=" << sublattice_size;
            cout << "\tswitch=" << search_switch;
            cout << "\tloops=" << precision_correction_loops_allowed;
            cout << "\teta=" << eta << endl;

            cout << "Min RHF:\ttime=" << runtime_min_RHF;
            cout << "\tRHF=" << RHF_min;
            cout << "\td=" << d;
            cout << "\talgo_num=" << algo_num;
            cout << "\tprec=" << precision_min_RHF;
            cout << "\tsublattice=" << sublattice_size;
            cout << "\tswitch=" << search_switch;
            cout << "\tloops=" << precision_correction_loops_allowed;
            cout << "\teta=" << eta << endl;

        }
        else
        {
            // The default case has been handled already
        }
    } 
    else 
    {
        cout << "Algorithm number should be 1 or 2" << endl;
        cout << "Algorithm 1 :  SS-GGLLL" << endl;
        cout << "Algorithm 2 :  Pot-GGLLL" << endl;
        exit (ERROR_ALGORITHM_NUMBER);
    }

    cout << endl << endl;
    cout << "Final   :\ttime=" << time_taken;
    cout << "\tRHF=" << RHF;
    cout << "\td=" << d;
    cout << "\talgo_num=" << algo_num;
    cout << "\tprec=" << precision;
    cout << "\tsublattice=" << sublattice_size;
    cout << "\tswitch=" << search_switch;
    cout << "\tloops=" << precision_correction_loops_allowed;
    cout << "\teta=" << eta << endl;

    return 0;
}
