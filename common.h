#ifndef __COMMON_H__
#define __COMMON_H__
#include <inttypes.h>
#include <complex.h>
#include <math.h>

#define PARAM_Q 4054119029
#define PARAM_K 32
#define PARAM_N 16
#define PARAM_R 2 // number of irreducible factors of x^n + 1 in F_q[x]
#define PARAM_D 2 // rank of the module
#define PARAM_SIGMA 10 // Gaussian parameter (generation of the trapdoor)
#define PARAM_ALPHA (60.50) // Gaussian parameter (sampling perturbations)
#define PARAM_ZETA (112522.0) // Gaussian parameter (presampling)
#define PARAM_T 15 // Tailcut
#define PARAM_B 2
#define PARAM_Lambda 1
#define PARAM_NumTrapdoors 1024*4
#define PARAM_NumSubvectors 1
#define NTT_ALWAYS_REDUCE
#define SMALL_DEGREE (PARAM_N / PARAM_R) // degree of irreducible factors of x^n + 1 in F_q[x]
#define PARAM_M (PARAM_D * (PARAM_K + 2)) // parameter m ( = d(k+2) for the computational instantiation)
//#define PARAM_ALPHA (3.0 * PARAM_SIGMA) // Gaussian parameter (sampling perturbations)



// LOG_N = base 2 log of n, defines the depth of the complex factorisation tree
#if (PARAM_N == 256)
	#define LOG_N 8
#elif (PARAM_N == 512)
	#define LOG_N 9
#elif (PARAM_N == 8)
	#define LOG_N 3
#elif (PARAM_N == 128)
	#define LOG_N  7
#elif (PARAM_N == 64)
	#define LOG_N 6
#elif (PARAM_N == 32)
	#define LOG_N 5
#elif (PARAM_N == 16)
	#define LOG_N 4
#else
	#error "LOG_N is not defined for this value of PARAM_N"
#endif

// LOG_R = base 2 log of r, defines the depth of the factorisation tree
#if (PARAM_R == 64)
	#define LOG_R 6
#elif (PARAM_R == 32)
	#define LOG_R 5
#elif (PARAM_R == 16)
	#define LOG_R 4
#elif (PARAM_R == 8)
	#define LOG_R 3
#elif (PARAM_R == 2)
	#define LOG_R 1
#else
	#error "LOG_R is not defined for this value of PARAM_R"
#endif

#if ((1 << LOG_N) != PARAM_N)
	#error "LOG_N is not base 2 log of PARAM_N"
#endif

#if ((1 << LOG_R) != PARAM_R)
	#error "LOG_R is not base 2 log of PARAM_R"
#endif

/*#if ((1 << PARAM_K) < PARAM_Q) || ((1 << (PARAM_K - 1)) >= PARAM_Q)
	#error "PARAM_K is not base 2 log of PARAM_Q"
#endif*/

#if (PARAM_K < 2)
	#error "PARAM_K < 2, not enough space in cplx_p_coeffs (sample_perturb)"
	#endif

//#define GET_BIT(x, i) (((x) >> (i)) & 1)


#define GET_BIT(x, i) get_bit_b(x, i)
#define Q_BIT(i) get_bit_b(PARAM_Q, i)

//#define DOUBLE_ZERO (((double_scalar) PARAM_Q) << (8*sizeof(double_scalar) - PARAM_K - 1)) // a double_scalar that is 0 mod q, such that there will be no underflow nor overflow during computations
#define DOUBLE_ZERO (((double_scalar) PARAM_Q) * PARAM_Q)  // a double_scalar that is 0 mod q, such that there will be no underflow nor overflow during computations
#define SIGNED_DOUBLE_ZERO (((double_scalar) PARAM_Q) << (8*sizeof(double_scalar) - PARAM_K - 2))

typedef uint64_t scalar;
typedef unsigned __int128 double_scalar;
typedef int64_t signed_scalar;
typedef __int128 signed_double_scalar;
typedef _Complex __float128  cplx;
//typedef double complex cplx;
typedef cplx *cplx_poly;
typedef cplx *cplx_poly_matrix;

//#if (PARAM_Q > 0x7FFFFFFF)
//	#error "Modulus PARAM_Q does not fit in 31 bits"
//#endif

/*
#if (LOG_N - LOG_R + 2 * PARAM_K) >= 63
	#error "(n/r)*q^2 >= 2^63, there might be overflow in polynomial multiplication"
#endif
*/
// Polynomial in standard (coefficients) representation
typedef scalar *poly;

// Non-normalized polynomial with double_scalar coefficients
typedef double_scalar *double_poly;

// Polynomial with signed coefficients (used in Karatsuba)
typedef signed_scalar *signed_poly;

// Polynomial with signed_double_scalar coefficients (used in Karatsuba)
typedef signed_double_scalar *signed_double_poly;

// Matrix of polynomials (single pointer on scalars, everything is accessed through macros)
typedef scalar *poly_matrix;

// Matrix of signed polynomials, used when sampling since Gaussian values are signed
typedef signed_scalar *signed_poly_matrix;

#define poly_matrix_element(M, nb_col, i, j) (&M[(PARAM_N)*(((i)*(nb_col)) + (j))])

#define crt_poly_component(f, deg, i) (&(f)[(i)*(deg)])

#define triangular_poly_matrix_element(M, i, j) (&M[(PARAM_N) * ((i)*((i)+1)/2 + (j))])

extern scalar cyclotomic_factorisation_array[2*PARAM_R - 1];
extern scalar *cyclotomic_factorisation_tree[LOG_R + 1];
extern scalar bezout_coefficients_array[2*PARAM_R - 1];
extern scalar *bezout_coefficients_tree[LOG_R + 1];
extern cplx cplx_roots_of_unity[2*PARAM_N - 1];

//#define cplx_root_of_unity_of_order(k, i) (cplx_roots_of_unity[(2 * (PARAM_N) / (k) * (i))])

#define cplx_cyclotomic_factorisation_tree(i, j) (cplx_roots_of_unity[(1 << (i)) + (j) - 1])

/*
	Stuff for Gaussian sampling
*/
typedef long double RR_t;
typedef double real;
#define LDRMX ((RR_t) RAND_MAX)
#define LOG_2 0.6931471805599453094172321214581765680755001343602552541206800094933936219696947156058633269964186875420014810205706857336855202357581305570326707516L
#define SIGMA_1 0.84932180028801904272150283410288961971514109378435394286159953238339383120795466719298223538163406787061691601172910413284884326532697308797136114023L //sqrt(1/(2*log(2)))

#endif
