#ifndef poly_h
#define poly_h

/*  
    (m,t) = (12,64), (13,96), (13,119), (13,128)
*/
#define m 12
#define t 64

#define MAX_COEF_POLY_DEGREE 2 * m + 1
#define MAX_POLY_DEGREE 2 * t + 1

#define TRUE 0
#define FALSE 1

#define IN
#define OUT

#include <stdio.h>
#include <memory.h>
#include <math.h>
#include "time.h"

//==========================  ????? ????  =================================
// typedef struct _COEF_POLY_ {  
//     int coef_max_degree;  
//     int coef[MAX_COEF_POLY_DEGREE+1]; 
// } COEF_POLY;

typedef struct poly {     
   int max_degree;         
   int coef[MAX_POLY_DEGREE+1];        
}POLY;

typedef struct _ctx_ {     
    POLY mod_gx;
    int mod_coef; 
} CTX;
//=========================================================================

void POLY_init(POLY* fx, int maxdeg);
void POLY_set(POLY* fx, int a[], int poly_deg);

void ctx_init(CTX* ctx);
void ctx_set(CTX* ctx, int ft, int gx[t+1], int deg_gx);
int POLY_is_zero(POLY* fx);
void set_POLY_degree(POLY* dst);

//int COEF_POLY_equal(COEF_POLY* src1, COEF_POLY* src2);
int POLY_equal(POLY* src1, POLY* src2);

//int COEF_is_one(COEF_POLY* dst);
//void COEF_POLY_copy(COEF_POLY* dst, COEF_POLY* src);
void POLY_copy(OUT POLY* dst, IN POLY* src);
//-- 
void COEF_POLY_print(int ft);
void POLY_print(POLY* fx);            //X
//-- 
void int2vec(OUT int* vec, OUT int* vec_size, IN int zz);
void vec2int(OUT int* zz, IN int* vec, IN int vec_size);
//== 
void coef_modft_table(OUT int ft_table[], IN CTX* ctx);
void coef_squ(OUT int* asqu,IN int ft_table[], IN int a[], IN int a_size, IN CTX* ctx);
void gen_Ttable(OUT int* Ttable, OUT int* InvTtable, IN int ft_table[], IN CTX* ctx);
//--
void gen_Xitable(OUT POLY Xtable[], IN CTX* ctx, IN int ft_table[]);

//== 
void COEF_POLY_mod_ft(OUT IN int* dst, IN CTX* ctx, IN int fttable[]);
void POLY_mod_gx(OUT IN POLY* dst, IN CTX* ctx, IN POLY Xtable[], IN int ft_table[]); 
void X_sqrt(OUT POLY* x_sqrt, IN POLY Xtable[], IN POLY* src,IN int fttable[],IN int* Ttable, IN CTX* ctx);
//-- 
void COEF_POLY_mul_zzx(OUT int* dst, IN int src, IN int ft_table[], IN CTX* ctx);
void COEF_POLY_mul(OUT int* dst, IN int src1, IN int src2, IN int ft_table[], IN CTX* ctx) ;
void POLY_mul(OUT POLY* dst, IN POLY* src1, IN POLY* src2, IN CTX* ctx, IN int fttable[], IN POLY Xtable[]); 
//-- 
void MULscalar_zzx(IN OUT POLY* dst, IN int src, IN int ft_table[], IN CTX* ctx);
void MULscalar(OUT POLY* dst, IN POLY* src1, IN int src2, IN int ft_table[], IN CTX* ctx);
//--

void POLY_add_zzx(OUT POLY* dst, IN POLY* src);
void POLY_add(OUT POLY* dst, IN POLY* src1, IN POLY* src2);
//-- ...
//6개 -> 구현. -> 결과. 

#endif  