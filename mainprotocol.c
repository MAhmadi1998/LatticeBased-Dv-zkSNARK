#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include "common.h"
#include "random.h"
#include "sampling.h"

#include "arithmetic.h"
#include "preprocessing.h"
#include "protocolparties.h"
#include "cpucycles.h"

#define CPU_CYCLES (2.3 * 1000000000.0)
unsigned long long timing_overhead;
unsigned long long timing_crsgen = 0;
unsigned long long timing_verify = 0;
unsigned long long prover_time = 0;
unsigned long long overall_proof = 0;
unsigned long long verification_time = 0;

unsigned long long sampling_timing = 0;


int main(void) {
	
	random_bytes_init();
	init_crt_trees();
	init_D_lattice_coeffs();
	init_cplx_roots_of_unity();
	timing_overhead = cpucycles_overhead();
	unsigned long long begin_timing = 0;
	unsigned long long end_timing = 0;

	printf("[*] Lattice-Based Designated Verifier zkSNARKs From Standard Assumptions\n");
	printf("[*] Note that this implementation is a proof-of-concept code\n");
	printf("[*] This code simulate commitent phase (Commitmen Key Gen + Commit algorithm)\n");
	printf("[*] Pok protocol for just lambda = 1\n");

	



 	//A_com is (PARAM_M - PARAM_D) * (NumOfTrapdoors* PARAM_D / 2) It would be in CRT Representation
	scalar *A_com_coefs = (scalar *)malloc(PARAM_N * (PARAM_M - PARAM_D)  *PARAM_NumTrapdoors * PARAM_D / 2 * sizeof(scalar));
	
	// A_hide is (PARAM_M - PARAM_D) * (NumOfTrapdoors* PARAM_D / 2) It would be in CRT Representation
	scalar *A_hide_coefs = (scalar *)malloc(PARAM_N * (PARAM_M - PARAM_D) * PARAM_NumTrapdoors * PARAM_D /2*sizeof(scalar));;
	
	poly_matrix A_com = A_com_coefs;
	poly_matrix A_hide = A_hide_coefs;

	// u_chall is in CRT domain
	scalar *u_chall = (scalar *)malloc(PARAM_N * (PARAM_NumTrapdoors * PARAM_NumSubvectors * PARAM_D) * sizeof(scalar));
	


	// z_chall is in CRT domin (The second round of protocol challenge)
	scalar *z_chall = (scalar *)malloc(PARAM_N * (PARAM_NumTrapdoors * PARAM_NumSubvectors* PARAM_D) * sizeof(scalar));
	
	signed_scalar *p_coeffs = (signed_scalar *)malloc(PARAM_N * (PARAM_M - PARAM_D) * sizeof(signed_scalar));
	signed_poly_matrix p = p_coeffs;



	printf("[*] Starting CRS Generation Algorithm for witness of lenght: %d ...\n" , PARAM_NumTrapdoors);
	

	begin_timing = cpucycles_start();
	CRS_generation_algorithm(A_com , A_hide , u_chall , z_chall , p);
	end_timing = cpucycles_stop();
	timing_crsgen = (end_timing -  begin_timing) - timing_overhead;
	printf("[*] CRS generation timing: %lld cycles (%.2lf ms)\n", timing_crsgen, (timing_crsgen*1000)/CPU_CYCLES);




	scalar *c_overal_coeffs = (scalar *)malloc(PARAM_N *(PARAM_M - PARAM_D)*sizeof(scalar));
	zero_poly(c_overal_coeffs , PARAM_N *(PARAM_M - PARAM_D));
	poly_matrix commitment = c_overal_coeffs;
	scalar *Inner_wu = (scalar *) malloc(PARAM_N *  sizeof(scalar)) ;
	zero_poly(Inner_wu , PARAM_N);
	scalar *Inner_wz = (scalar *) malloc(PARAM_N *  sizeof(scalar)) ;
	zero_poly(Inner_wz , PARAM_N);

	scalar z_sec2[PARAM_N ];
	scalar z_sec3[PARAM_N ];
	scalar z_sec4[PARAM_N ];
	random_poly(z_sec2 , PARAM_N);
	random_poly(z_sec3 , PARAM_N);
	printf("[*] Running Prover algorithm: Commitment + Response to 2 PoK Challenges ...\n");
	prover(A_com, A_hide , u_chall , z_chall , commitment , Inner_wu , Inner_wz);
	printf("\n [*] The value of u_chall.[b] is equal to:\n                          [mu]\n");
	print_poly_matrix(Inner_wu , 1,1);
  
	printf("[*] Starting Verification of one inner product challenge ...\n");
	printf("\n[*] The value of <p,c> + z_agg.[b] is equal to:\n                               [mu]\n");
	scalar *PtimesC = (scalar *) malloc(PARAM_N *  sizeof(scalar)) ; 
	//scalar PtimesC [PARAM_N];
	begin_timing = cpucycles_start();
	mul_normal_poly_matrix(PtimesC, (poly_matrix)p , commitment , 1 , PARAM_M - PARAM_D , 1 );
	end_timing = cpucycles_stop();
	timing_verify = (end_timing -  begin_timing) - timing_overhead;
	verification_time = verification_time + timing_verify;


	begin_timing = cpucycles_start();
	add_to_poly_matrix(Inner_wz , PtimesC , 1,1 );
	end_timing = cpucycles_stop();
	timing_verify = ((end_timing -  begin_timing) - timing_overhead);
	verification_time = verification_time + timing_verify;
	freeze_poly(Inner_wz , PARAM_N);
	
	//printf("\n ----------------------------------------- <b , z> + <p, c> is --------------------------------------------\n");

	print_poly_matrix(Inner_wz , 1,1);
	// Norm Infty Calc
	printf("[*] Verification timing %lld cycles (%.2lf ms)\n", verification_time, (verification_time*1000)/CPU_CYCLES);
	


	
	
	
	

	free(A_com_coefs);
	free(A_hide_coefs);
	free(u_chall);
	free(z_chall);
	free(p);
	free(Inner_wu);
	free(Inner_wz);
	free(PtimesC);
	//free(p_coeffs);

	return 0;
}

