#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include "common.h"
#include "random.h"
#include "sampling.h"
#include "protocolparties.h"
#include "cpucycles.h"
#include "arithmetic.h"

#define CPU_CYCLES (2.3 * 1000000000.0)
unsigned long long timing_overhead1;
unsigned long long timing_commit = 0;
unsigned long long timing_proofGen = 0;




void prover(poly_matrix A_com , poly_matrix A_hide , poly_matrix u_chall , poly_matrix z_chall , poly_matrix commitment , poly_matrix Inner_wu , poly_matrix Inner_wz)
{
    //scalar c_overal_coeffs[PARAM_N *(PARAM_M - PARAM_D)];
	//zero_poly(commitment , PARAM_N *(PARAM_M - PARAM_D));
	//poly_matrix c_overal = c_overal_coeffs;

	//scalar buOverall[PARAM_N];
	//zero_poly(Inner_wu , PARAM_N);
	//scalar bzOverall[PARAM_N];
	//zero_poly(Inner_wz , PARAM_N);
	unsigned long long begin_timing = 0;
	unsigned long long end_timing = 0;
	timing_overhead1 = cpucycles_overhead();
	for (int subs = 0; subs < PARAM_NumSubvectors; subs++)
	{
		scalar *b_coefs = (scalar *) malloc(PARAM_N*PARAM_D*PARAM_NumTrapdoors / 2 * sizeof(scalar));
		scalar *mu_coefs = (scalar *) malloc(PARAM_N*PARAM_D*PARAM_NumTrapdoors / 2 * sizeof(scalar));
		poly_matrix b_vec = b_coefs;
		poly_matrix mu_vec = mu_coefs;
		random_poly(b_vec , PARAM_N*PARAM_D*PARAM_NumTrapdoors / 2);

		random_poly(mu_vec , PARAM_N*PARAM_D*PARAM_NumTrapdoors / 2);
		
		scalar c_com_coefs [PARAM_N * (PARAM_M -PARAM_D)];
		scalar c_hide_coefs [PARAM_N * (PARAM_M -PARAM_D)];
		poly_matrix c_com = c_com_coefs;
		poly_matrix c_hide = c_hide_coefs;
		begin_timing = cpucycles_start();
		mul_normal_poly_matrix(c_com,A_com , b_vec , PARAM_M - PARAM_D , PARAM_NumTrapdoors*PARAM_D/2 , 1 );
		mul_normal_poly_matrix(c_hide,A_hide , mu_vec , PARAM_M - PARAM_D , PARAM_NumTrapdoors*PARAM_D/2 , 1);
		add_to_poly_matrix(c_com , c_hide , PARAM_M - PARAM_D , 1);
		freeze_poly(c_com,  PARAM_N*(PARAM_M - PARAM_D));
		add_to_poly_matrix(commitment , c_com , PARAM_M - PARAM_D, 1);
		freeze_poly(commitment, (PARAM_M - PARAM_D)*PARAM_N);
		end_timing = cpucycles_stop();
		timing_commit = timing_commit + (end_timing - begin_timing);
		
		scalar BTimesU[PARAM_N];
		scalar MuTimesU[PARAM_N];
		scalar BTimesZ[PARAM_N];
		scalar MuTimesZ[PARAM_N];
		begin_timing = cpucycles_start();
		poly_matrix u_chall_com = poly_matrix_element(u_chall, PARAM_D*PARAM_NumSubvectors*PARAM_NumTrapdoors , 0 , subs*PARAM_NumTrapdoors*PARAM_D/2);
		mul_normal_poly_matrix(BTimesU , u_chall_com , b_vec , 1,PARAM_D/2*PARAM_NumTrapdoors , 1 );
		poly_matrix u_chall_hide = poly_matrix_element(u_chall ,PARAM_NumSubvectors*PARAM_NumTrapdoors * PARAM_D , 0 , PARAM_NumSubvectors*PARAM_NumTrapdoors*PARAM_D / 2 +subs*PARAM_NumTrapdoors*PARAM_D/2 );
		mul_normal_poly_matrix(MuTimesU ,u_chall_hide , mu_vec , 1 , PARAM_D*PARAM_NumTrapdoors/2 , 1 );
		add_to_poly_matrix(BTimesU , MuTimesU , 1,1);
		freeze_poly(BTimesU , PARAM_N);
		add_to_poly_matrix(Inner_wu , BTimesU , 1,1);
		freeze_poly(Inner_wu , PARAM_N);

		
		poly_matrix z_chall_com = poly_matrix_element(z_chall , PARAM_NumSubvectors*PARAM_NumTrapdoors*PARAM_D,0 ,subs*PARAM_NumTrapdoors*PARAM_D/2 );
		mul_normal_poly_matrix(BTimesZ , z_chall_com , b_vec , 1,PARAM_D*PARAM_NumTrapdoors/2 , 1 );
		poly_matrix z_chall_hide = poly_matrix_element(z_chall ,PARAM_NumSubvectors*PARAM_NumTrapdoors * PARAM_D , 0 , PARAM_NumSubvectors*PARAM_NumTrapdoors*PARAM_D / 2 +subs*PARAM_NumTrapdoors*PARAM_D/2);
		mul_normal_poly_matrix(MuTimesZ ,z_chall_hide , mu_vec , 1 , PARAM_D*PARAM_NumTrapdoors/2 , 1 );
		add_to_poly_matrix(BTimesZ , MuTimesZ , 1,1);
		freeze_poly(BTimesZ , PARAM_N);
		add_to_poly_matrix(Inner_wz , BTimesZ , 1,1);
		freeze_poly(Inner_wz , PARAM_N);
		end_timing = cpucycles_stop();
		timing_proofGen= timing_proofGen + (end_timing - begin_timing);
	}
	timing_commit = timing_commit - timing_overhead1;
	timing_proofGen = (timing_proofGen - timing_overhead1);
	printf("[*] commitment timing: %lld cycles (%.2lf ms)\n", timing_commit, (timing_commit*1000)/CPU_CYCLES);
	printf("[*] proof generation timing: %lld cycles (%.2lf ms)\n", timing_proofGen, (timing_proofGen*1000)/CPU_CYCLES);

}

