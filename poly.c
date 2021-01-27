#include "poly.h"


//======================================================================
//
//1 + X + X^3 + X^5 + X^6 + X^7 + X^12  (1,1,0,1,0,1,1,1,0,0,0,0,1)
//m=12  
//f(X) = 1 + X + X^3 + X^4 + X^13       (1,1,0,1,1,0,0,0,0,0,0,0,0,1)
//m=13
//x^96+x^6+x^5+x^3+x^2+x+1   
//t=96
//
//

//=====================================================================

void COEF_POLY_init(COEF_POLY* ft, int maxdeg)
{
    ft->coef_max_degree = maxdeg;
    memset(ft->coef,0,(MAX_COEF_POLY_DEGREE+1)*sizeof(int));
}

void POLY_init(POLY* fx, int maxdeg)
{ 
    fx->max_degree = maxdeg;
    for (int i = 0; i < MAX_POLY_DEGREE + 1; i++)
        COEF_POLY_init(&fx->coef[i],0);
}

void COEF_POLY_set(COEF_POLY* ft, int* a, int deg)
{
    ft->coef_max_degree = deg;
    memcpy(ft->coef, a, (deg + 1)*sizeof(int));
}

void POLY_set(POLY* fx, int a[][MAX_COEF_POLY_DEGREE+1], int poly_deg, int coef_deg)
{  
    fx->max_degree = poly_deg;
    for (int i = 0; i < poly_deg+1; i++)
        COEF_POLY_set(&fx->coef[i], a[i], coef_deg);
}

void ctx_init(CTX* ctx)
{
    POLY_init(&ctx->mod_gx,0);
    COEF_POLY_init(&ctx->mod_coef,0);
}

void ctx_set(CTX* ctx, int ft[m+1], int gx[t+1], int deg_ft, int deg_gx)
{ 
    COEF_POLY_set(&ctx->mod_coef, ft, deg_ft);
    
    ctx->mod_gx.max_degree = deg_gx;
    for (int i = 0; i < MAX_POLY_DEGREE; i++)
        if(gx[i] == 1)
        {
            ctx->mod_gx.coef[i].coef[0] = 1;
            ctx->mod_gx.coef[i].coef_max_degree = 0;
        }
}

//===================================================================================

void COEF_POLY_copy(COEF_POLY* dst, COEF_POLY* src)
{  
    dst->coef_max_degree = src->coef_max_degree;
    memcpy(dst->coef, src->coef, (src->coef_max_degree+1)*sizeof(int));
}

void POLY_copy(OUT POLY* dst, IN POLY* src)
{
    dst->max_degree = src->max_degree;
    for (int i = 0; i < dst->max_degree+1; i++)
        COEF_POLY_copy(&dst->coef[i], &src->coef[i]);
}

void COEF_POLY_print(COEF_POLY* ft){
    int first = 1;                           // 0이 아닌 수가 처음인지 아닌지 판별하기 위한 변수
    int max_deg = ft->coef_max_degree;       // max_deg 설정
    while(max_deg >= 0){                     
        if(ft->coef[max_deg] != 0){          // 해당 수가 0이 아니라면 실행
            if(first == 0){                  // 0이 아닌 수가 처음이 아니라면 '+'기호를 먼저 print함
                printf(" + ");
            }
            
            if(max_deg == 0){                // 0차의 수는 상수이므로 상수만 출력
                printf("%d ", ft->coef[max_deg]);
            }
            else if(max_deg == 1){           // 1차의 수는 () * t로 출력
                printf("%d * t ", ft->coef[max_deg]);
            }
            else{                            // n차의 수는 () * t^n로 출력 (n > 1)
                printf("%d * t^%d ", ft->coef[max_deg], max_deg);
            }
            first = 0;                       // 한 번이라도 if문이 실행되면 first를 0으로 변경
        }
        else{
            if(ft->coef_max_degree == 0)     // ft가 0이라면 그대로 0만 출력
                printf("%d ", ft->coef[max_deg]);
        }
        max_deg -= 1;                      
    }
}

int COEF_is_zero(COEF_POLY* ft) // ft가 0인지 아닌지 판별하는 함수, 0이라면 '0' 출력, 아니라면 '1'출력
{             
    int zero = 0;                            // ft가 0이라면 마지막까지 변하지 않음
    int max_deg = ft->coef_max_degree;
    while(max_deg >= 0){                     
        if(ft->coef[max_deg] != 0){          // 0이 아닌경우 zero에 1을 추가하고 바로 break
            zero += 1;
            break;
        }
        max_deg -= 1;
    }

    if(zero == 0)       return TRUE;
    else                return FALSE;
}

int POLY_is_zero(POLY* fx)
{
    if(fx->max_degree == 0 && COEF_is_zero(&fx->coef[0]) == TRUE)   return TRUE;
    else    return FALSE;
}

void POLY_print(POLY* fx)
{  
    if(POLY_is_zero(fx) == TRUE)
    {
        printf("0\n");
        return;
    }
    int first = 1;                                           // 0이 아닌 수가 처음인지 아닌지 판별하기 위한 변수
    int max_deg = fx->max_degree;
    while(max_deg >= 0)
    {
        //printf("\n %d \n", COEF_is_zero(&fx->coef[max_deg]));
        if(COEF_is_zero(&fx->coef[max_deg]) != 0)
        {           // 해당 수가 0이 아니라면 실행
            if(first == 0){                                  // 0이 아닌 수가 처음이 아니라면 '+'기호를 먼저 print함
                printf(" + ");
            }

            if(max_deg == 0){                                // 0차의 수는 상수이므로 상수만 출력
                printf("(");
                COEF_POLY_print(&fx->coef[max_deg]);
                printf(")");
            }
            else if(max_deg == 1){                            // 1차의 수는 () * X로 출력
                printf("(");
                COEF_POLY_print(&fx->coef[max_deg]);
                printf(") * X ");
            }
            else{                                             // n차의 수는 () * X^n로 출력 (n > 1)
                printf("(");
                COEF_POLY_print(&fx->coef[max_deg]);
                printf(") * X^%d ", max_deg);
            }
            first = 0;                                        // 한 번이라도 if문이 실행되면 first를 0으로 변경
        }
        else{
            if(fx->max_degree == 0){                          // ft가 0이라면 그대로 0만 출력
                COEF_POLY_print(&fx->coef[max_deg]);
            }
        }
        max_deg -= 1;
    }
}

void int2vec(OUT int* vec, OUT int* vec_size, IN int zz){ 

}

void vec2int(OUT int* zz, IN int* vec, IN int vec_size){

}

void coef_modft_table(OUT int* ft_table, IN CTX* ctx){  

}

void coef_squ(OUT int* asqu, IN COEF_POLY* a, IN CTX* ctx){  

}

void gen_Ttable(OUT int* Ttable, OUT int* InvTtable, IN CTX* ctx ){  

}

void ModExpX_i(OUT POLY* xi, IN POLY* x, IN int i, IN CTX* ctx){ 

}

void gen_Xitable(OUT int** Xtable, IN CTX* ctx){ 

}

void POLY_mod_gx(OUT POLY* dst, IN POLY* src, IN CTX* ctx){  
    
}

void X_sqrt(OUT POLY* x_sqrt, IN POLY* x, IN CTX* ctx){  //x^i

}

void COEF_POLY_mul(OUT COEF_POLY* ht,IN COEF_POLY* ft, IN COEF_POLY* gt, IN CTX* ctx) {

}

void POLY_MUL(OUT POLY* dst, IN POLY* src1, IN POLY* src2, IN CTX* ctx){

}

void MULscalar(OUT POLY* dst, IN POLY* src, IN int a, IN CTX ctx){

}

void COEF_POLY_add(COEF_POLY* ht,COEF_POLY* ft, COEF_POLY* gt, IN CTX ctx) {

}

void POLY_ADD(OUT POLY* dst, IN POLY* src1, IN POLY* src2, IN CTX ctx){

}










/*

void COEF_POLY_add(COEF_POLY* ht,COEF_POLY* ft, COEF_POLY* gt)   // ???????? ?????? xor?? ???
{
    POLY_init(ht);
    int mi ;
    mi=(ft->coef_max_degree < gt->coef_max_degree)? ft->coef_max_degree:gt->coef_max_degree;
    int i;
    for(i=0 ; i < mi;i++)
    {
        (*ht)->coef[i]=ft->coef[i]^gt->coef[i];    
    }
    if (mi==ft->coef_max_degree)
    {
        for(i=mi;i<gt->coef_max_degree;i++)
            (*ht)->coef[i]=gt->coef[i]; 
    }
    else
    {
        for(i=mi;i<ft->coef_max_degree;i++)
            (*ht)->coef[i]=ft->coef[i]; 
    }
}


void COEF_POLY_mul(COEF_POLY** ht,COEF_POLY* ft, COEF_POLY* gt) {

}

void set_POLY_zero(IN POLY** fx){   // fx ?? 0?? ????????? ??????. 


}


void POLY_set(POLY** fx, int a[]){      // ????? ?????? ????? ??????.???????. 


}


void POLY_delete(POLY** fx) {    // ????? ????

}


void ADDpoly(OUT POLY hx, IN POLY fx, IN POLY gx) {  // fx + gx = hx

}


void MULpoly(OUT POLY hx, IN POLY fx, IN POLY gx){      //?????? ????  fx * gx = hx

}


void gen_Ttable(OUT int* c, IN OUT int* a,IN int b, IN int mod_coef){    // T type table   a???????θ?? b???????
    //????? ????? 8192?? ????. 


    for(int i =0; i<pow(2,m);i++)
    {
        a[i]=pow(a[i],b);
        a[i]=a[i]%mod_coef;        
    }

// 2mt-1 ?????? ?????? ????? ?????? 
    for(int i=0;i<pow(2,m);i++)
    {
        if()

    }
}


void gen_Rtable(POLY fx, POLY gx) {    // R type table   a???????θ?? b???????

}


void squX(OUT POLY hx, IN POLY fx, IN POLY gx, IN int n) {  // fx^n mod gx = hx

}


void MULscalar(POLY* fx,int a){    // ?????? ?????? ????.

    for(int i =0; i<t+1;i++)
    {
        fx->coef[i]=(fx->coef[i])*(a)%(fx->mod_coef);
    }
}


void modulo(POLY fx, POLY gx){   //fx (mod gx)

}

*/