#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include "common.h"
#include "random.h"
#include "sampling.h"
#include "preprocessing.h"

#include "arithmetic.h"



void CRS_generation_algorithm(poly_matrix A_com, poly_matrix A_hide, poly_matrix u_chall,  poly_matrix z_chall, signed_poly_matrix p)
{
    
	//scalar *A_com_coefs = (scalar *)malloc(PARAM_N * (PARAM_M - PARAM_D)  *PARAM_NumTrapdoors * PARAM_D / 2 * sizeof(scalar));
	
    scalar *T_coeffs_tmp = (scalar *)malloc(PARAM_N * PARAM_D * PARAM_D * PARAM_K* sizeof(scalar));
	poly_matrix T_tmp = T_coeffs_tmp;
	
	cplx *sch_comp_coeffs_tmp = (cplx *) malloc(PARAM_N * PARAM_D * (2 * PARAM_D + 1) * sizeof(cplx));
	cplx *cplx_T_coeffs_tmp = (cplx *) malloc(PARAM_N * PARAM_D * PARAM_D * PARAM_K* sizeof(cplx)) ;
	cplx_poly_matrix cplx_T_tmp = cplx_T_coeffs_tmp;
	cplx_poly_matrix sch_comp_tmp = sch_comp_coeffs_tmp;
	
	SampleR_matrix_centered((signed_poly_matrix) T_tmp, PARAM_D, PARAM_D * PARAM_K, PARAM_SIGMA);
	
	construct_complex_private_key(cplx_T_tmp,sch_comp_tmp , T_tmp);
	
	
	
	sample_perturb(p, cplx_T_tmp, sch_comp_tmp);	
	
	
	// Add q to p's coeffs so that they are positive, and put p in the CRT domain
	for(int i = 0 ; i < PARAM_N * (PARAM_M - PARAM_D) ; ++i)
		{
		p[i] += PARAM_Q;
		}

	freeze_poly((poly_matrix)p , PARAM_N *(PARAM_M - PARAM_D));
	
	//matrix_crt_representation((poly_matrix) p, PARAM_M - PARAM_D, 1, LOG_R);
	
	
	// Defining Gadget Matrix
	scalar *G_coefs = (scalar *)malloc(PARAM_N * PARAM_D * PARAM_D * PARAM_K * sizeof(scalar));

	poly_matrix G_orginal = G_coefs;
	int r = 0;
	for (size_t i = 0; i < PARAM_D; i++)
	{
		r = 0 ;
		for (size_t j = 0; j < PARAM_D*PARAM_K; j++)
		{	
			if (j < i*PARAM_K || j > (i+1)*PARAM_K - 1)
			{
				for (size_t k = 0; k < PARAM_N; k++)
				{
					G_orginal[i*PARAM_D*PARAM_K*PARAM_N + PARAM_N*j +k] = 0;
				
				}
			}

			else
			{
				G_orginal[i*PARAM_D*PARAM_K*PARAM_N + PARAM_N*j] = reduce_double_naive((double_scalar)pow(PARAM_B ,r));
				r++;
				for (size_t k = 1; k < PARAM_N; k++)
				{
					G_orginal[i*PARAM_D*PARAM_K*PARAM_N + PARAM_N*j +k] = 0;
				
				}
			}
			

			
		}
		
	}


	
	scalar *G_coefs_me = (scalar *)malloc(PARAM_N * PARAM_D * PARAM_D * PARAM_K * sizeof(scalar));
	
	poly_matrix G_me= G_coefs_me;
	
	scalar *A_coeffs = (scalar *)malloc(PARAM_N * PARAM_D * (PARAM_M - PARAM_D)* sizeof(scalar));
	scalar  *T_coeffs = (scalar *)malloc(PARAM_N * PARAM_D * PARAM_D * PARAM_K * sizeof(scalar));
	poly_matrix A = A_coeffs, T = T_coeffs;
	cplx *sch_comp_coeffs = (cplx *) malloc(PARAM_N * PARAM_D * (2 * PARAM_D + 1) * sizeof(cplx));
	cplx *cplx_T_coeffs = (cplx *) malloc(PARAM_N * PARAM_D * PARAM_D * PARAM_K * sizeof(cplx));
	cplx_poly_matrix cplx_T = cplx_T_coeffs;
	cplx_poly_matrix sch_comp = sch_comp_coeffs;

	scalar *synd_coefs = (scalar *) malloc(PARAM_N * (PARAM_D)* sizeof(scalar));
	poly_matrix Synd_Vec = synd_coefs;

	for (int counter = 0; counter < PARAM_NumTrapdoors; counter++)
	{

	
	
	memmove(G_me , G_orginal , PARAM_N*PARAM_D*PARAM_D*PARAM_K*(sizeof(scalar)));
	
	
	// Sampling A_i and corresponding R_i


	SampleR_matrix_centered((signed_poly_matrix) T, PARAM_D, PARAM_D * PARAM_K, PARAM_SIGMA);
	
	construct_complex_private_key(cplx_T,sch_comp , T);
	
	
	TrapGen(A, T , G_me);

	

	
	for (size_t nrows = 0; nrows < PARAM_D / 2; nrows++)
	{
		poly_matrix A_com_cols =  poly_matrix_element(A_com , PARAM_D/2 *PARAM_NumTrapdoors , 0 , counter*PARAM_D /2 + nrows);
		poly_matrix A_hide_cols =  poly_matrix_element(A_hide , PARAM_D/2 *PARAM_NumTrapdoors , 0 , counter*PARAM_D /2 + nrows);
		poly_matrix A_i_row = poly_matrix_element(A , (PARAM_M - PARAM_D) , nrows , 0);
		poly_matrix A_i_row_hide = poly_matrix_element(A , (PARAM_M - PARAM_D) , nrows + PARAM_D/2, 0);
		for (size_t ind = 0; ind < (PARAM_M - PARAM_D); ind++)
		{
			for (size_t degs = 0; degs < PARAM_N; degs++)
			{
				A_com_cols[PARAM_N*((PARAM_D/2)*PARAM_NumTrapdoors)*ind + degs] = A_i_row[PARAM_N*ind + degs];
				A_hide_cols[PARAM_N*((PARAM_D/2)*PARAM_NumTrapdoors)*ind + degs] = A_i_row_hide[PARAM_N*ind + degs];
			}
			
		}
		
	}
	
	
	

	// Sampling Challenges
	for (size_t Counter_Subvec = 0; Counter_Subvec < PARAM_NumSubvectors; Counter_Subvec++)
	{
		
		
		random_poly(Synd_Vec , PARAM_N * ( PARAM_D));
		
		
		
		poly_matrix u_com = poly_matrix_element(u_chall , PARAM_D , 0 , counter*(PARAM_D / 2)+ Counter_Subvec*PARAM_NumTrapdoors*PARAM_D/2);

		
		
		poly_matrix u_hide = poly_matrix_element(u_chall , PARAM_D , 0 , counter*(PARAM_D / 2) + Counter_Subvec*PARAM_NumTrapdoors*PARAM_D/2 + PARAM_NumSubvectors*PARAM_NumTrapdoors*PARAM_D/2);
	
		memmove(u_com , Synd_Vec , (PARAM_D / 2)*PARAM_N * sizeof(scalar));
	
		
		poly_matrix secondhalf_Vec  = poly_matrix_element(Synd_Vec , PARAM_D , 0 , PARAM_D / 2);
		
		memmove(u_hide , secondhalf_Vec , (PARAM_D / 2)*PARAM_N * sizeof(scalar));
		
		scalar *x_coefs = (scalar *)malloc(PARAM_N * (PARAM_M - PARAM_D)* sizeof(scalar));
		poly_matrix x = x_coefs;

		signed_scalar *z_coefs = (signed_scalar *)malloc(PARAM_D*PARAM_K*PARAM_N * sizeof(signed_scalar));
		signed_poly_matrix z = z_coefs;
	
		sample_pre(x,A,T,cplx_T,sch_comp , Synd_Vec, (poly_matrix) p , z);


		scalar *z_agg_com_coefs = (scalar *)malloc(PARAM_N * (PARAM_D/2) * sizeof(scalar));
		poly_matrix z_agg_com = z_agg_com_coefs;
		mul_normal_poly_matrix(z_agg_com , G_orginal , (poly_matrix)z , PARAM_D /2 , PARAM_D * PARAM_K, 1);
		scalar *z_agg_hide_coefs = (scalar *)malloc(PARAM_N * (PARAM_D/2)*sizeof(scalar)) ;
		poly_matrix z_agg_hide = z_agg_hide_coefs;
		poly_matrix G_down_half = poly_matrix_element(G_orginal , PARAM_D * PARAM_K , PARAM_D / 2 , 0);
		mul_normal_poly_matrix(z_agg_hide ,G_down_half , (poly_matrix) z , PARAM_D / 2 , PARAM_D * PARAM_K , 1 );

		poly_matrix z_com = poly_matrix_element(z_chall , PARAM_D , 0 , counter*(PARAM_D / 2) + Counter_Subvec*PARAM_NumTrapdoors*PARAM_D/2);
		poly_matrix z_hide = poly_matrix_element(z_chall , PARAM_D , 0 , counter*(PARAM_D / 2) + Counter_Subvec*PARAM_NumTrapdoors*PARAM_D/2 +  PARAM_NumSubvectors*PARAM_NumTrapdoors*PARAM_D/2);
		memmove(z_com , z_agg_com , sizeof(scalar)*PARAM_N*PARAM_D/2);
		memmove(z_hide , z_agg_hide , sizeof(scalar)*PARAM_N*PARAM_D/2);
		free(x_coefs);
		free(z_coefs);
		free(z_agg_com_coefs);
	}
	

	
	
	}
	free(T_coeffs_tmp);
	free(sch_comp_coeffs_tmp);
	free(cplx_T_coeffs_tmp);
	free(G_coefs);
	free(G_coefs_me);
	free(A_coeffs);
	free(T_coeffs);
	free(sch_comp_coeffs);
	free(cplx_T_coeffs);
	free(synd_coefs);

}