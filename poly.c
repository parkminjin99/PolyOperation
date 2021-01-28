#include "poly.h"
#define MAX(x,y)( (x)>(y)?(x):(y) ) 

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
                printf("+ ");
            }

            if(max_deg == 0){                                // 0차의 수는 상수이므로 상수만 출력
                //printf("(");
                COEF_POLY_print(&fx->coef[max_deg]);
                //printf(")");
            }
            else if(max_deg == 1){                            // 1차의 수는 () * X로 출력
                if(fx->coef[max_deg].coef_max_degree == 0 && fx->coef[max_deg].coef[0] == 1)
                    printf("X ");
                else
                {
                    printf("(");
                    COEF_POLY_print(&fx->coef[max_deg]);
                    printf(") * X ");
                }
            }
            else{                                             // n차의 수는 () * X^n로 출력 (n > 1)
                if(fx->coef[max_deg].coef_max_degree == 0 && fx->coef[max_deg].coef[0] == 1)
                    printf("X^%d ", max_deg);
                else
                {
                    printf("(");
                    COEF_POLY_print(&fx->coef[max_deg]);
                    printf(") * X^%d ", max_deg);   
                }
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
    printf("\n");
}

void int2vec(OUT int* vec, OUT int* vec_size, IN int zz) { // 정수를 벡터로
    int count = 0;
    
    for(int i=0;i<m+1;i++)
        vec[i]=0;
    while (1) {
        vec[count] = zz % 2;
        zz = zz / 2;
        if (zz == 0)
            break;
        count++;
    }
    *vec_size = count+1;
}

void vec2int(OUT int* zz, IN int* vec, IN int vec_size) { // 벡터를 정수로
    *zz = 0;
    for (int i = 0; i < vec_size; i++)
    {
        *zz += vec[i] * pow(2, i);
    }
}

void coef_modft_table(OUT int* ft_table, IN CTX* ctx){  
    
    if(m==12)  // (1,1,0,1,0,1,1,1,0,0,0,0,1)   
    {
        int tmp;
        int vec[m+1]={0,0,0,0,0,0,0,0,0,0,0,0,1};

        for(int i=0;i<m;i++)  // t^12 ~ t^22
        {            
            if(vec[12]==1)
            {
                vec[12]=0;
                vec[0]^=1;
                vec[1]^=1;
                vec[3]^=1;
                vec[5]^=1;
                vec[6]^=1;
                vec[7]^=1;
            }
            vec2int(&ft_table[i],vec,12);
            tmp=vec[m];
            for(int j=m-1;j>=0;j--)
            {
                vec[j+1]=vec[j];
            }
            vec[0]=tmp;
        }
    }
    else  // m=13 (1,1,0,1,1,0,0,0,0,0,0,0,0,1)
    {
        int tmp;
        int vec[m+1]={0,0,0,0,0,0,0,0,0,0,0,0,0,1};

        for(int i=0;i<m;i++)  // t^13 ~ t^24
        {            
            if(vec[13]==1)
            {
                vec[13]=0;
                vec[0]^=1;
                vec[1]^=1;
                vec[3]^=1;
                vec[4]^=1;
            }
            vec2int(&ft_table[i],vec,14);
            printf("%d ",ft_table[i]);
            for(int j=0;j<m+1;j++)
                printf("%d",vec[j]);
            printf("\n");  
            for(int j=m-1;j>=0;j--)
            {
                vec[j+1]=vec[j];
            }
            vec[0]=0;
        
  
        }
    }
}

void coef_squ(OUT int* asqu, IN COEF_POLY* a, IN CTX* ctx){  

}

void gen_Ttable(OUT int* Ttable, OUT int* InvTtable, IN CTX* ctx ){  

}

void ModExpX_i(OUT POLY* xi, IN POLY* x, IN int i, IN CTX* ctx){ 

}

void COEF_POLY_add_zzx(OUT COEF_POLY* dst, IN COEF_POLY* src, IN CTX* ctx)
{
    for(int i = 0; i <= src->coef_max_degree; i++)
    {
        dst->coef[i] ^= src->coef[i];
        if(dst->coef[i] != 0)   dst->coef_max_degree = i;
    }
}

void POLY_add_zzx(OUT POLY* dst, IN POLY* src, IN CTX* ctx) // t차 이상의 다항식은 안들어온다는 가정하에 덧셈 
{
    for (int i = 0; i < t; i++)
    {
        COEF_POLY_add_zzx(&dst->coef[i], &src->coef[i], ctx);
        if(COEF_is_zero(&dst->coef[i]) != TRUE)     dst->max_degree = i;  
    }
}

void gen_Xitable(OUT POLY Xtable[], IN CTX* ctx)
{ 
    POLY_copy(&Xtable[0],&ctx->mod_gx);
    COEF_POLY_init(&Xtable[0].coef[t], 0);
    for(int i = 1; i <= t; i++) // Xtable[i] POLY 생성 
    {
        POLY_init(&Xtable[i], 0);
        for (int j = 1; j <= t-1; j++)
        {
            COEF_POLY_copy(&Xtable[i].coef[j], &Xtable[i-1].coef[j-1]);
            if(COEF_is_zero(&Xtable[i].coef[j]) != TRUE)
                Xtable[i].max_degree = j;
        }
        if(COEF_is_zero(&Xtable[i-1].coef[t-1]) != TRUE)
            POLY_add_zzx(&Xtable[i], &Xtable[0], ctx);
    }
}

void POLY_mod_gx(OUT POLY* dst, IN POLY* src, IN CTX* ctx)
{  
    int max_deg = dst->max_degree;
    //for(int i = t; i <  )
}

void X_sqrt(OUT POLY* x_sqrt, IN POLY* x, IN CTX* ctx){  //x^i

}

void COEF_POLY_mul(OUT COEF_POLY* ht,IN COEF_POLY* ft, IN COEF_POLY* gt, IN CTX* ctx) {

}

void POLY_MUL(OUT POLY* dst, IN POLY* src1, IN POLY* src2, IN CTX* ctx){

}

void MULscalar(OUT POLY* dst, IN POLY* src, IN int a, IN CTX ctx){

}


void POLY_ADD(OUT POLY* dst, IN POLY* src1, IN POLY* src2, IN CTX ctx){

}


/*




void COEF_POLY_add(COEF_POLY* ht,COEF_POLY* ft, COEF_POLY* gt)   // 계수간의 덧셈은 xor로 처리
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

void set_POLY_zero(IN POLY** fx){   // fx 를 0인 다항식으로 만들기. 


}


void POLY_set(POLY** fx, int a[]){      // 구조체 계수를 원하는 값으로.설정하기. 


}


void POLY_delete(POLY** fx) {    // 구조체 삭제

}


void ADDpoly(OUT POLY hx, IN POLY fx, IN POLY gx) {  // fx + gx = hx

}


void MULpoly(OUT POLY hx, IN POLY fx, IN POLY gx){      //다항식간 곱셈  fx * gx = hx

}


void gen_Ttable(OUT int* c, IN OUT int* a,IN int b, IN int mod_coef){    // T type table   a개원소모두를 b제곱하는
    //테이블 크기는 8192로 고정. 


    for(int i =0; i<pow(2,m);i++)
    {
        a[i]=pow(a[i],b);
        a[i]=a[i]%mod_coef;        
    }

// 2mt-1 제곱이 빠를까 서치가 빠를까 
    for(int i=0;i<pow(2,m);i++)
    {
        if()

    }
}


void gen_Rtable(POLY fx, POLY gx) {    // R type table   a개원소모두를 b제곱하는

}


void squX(OUT POLY hx, IN POLY fx, IN POLY gx, IN int n) {  // fx^n mod gx = hx

}


void MULscalar(POLY* fx,int a){    // 다항식과 상수의 곱셈.

    for(int i =0; i<t+1;i++)
    {
        fx->coef[i]=(fx->coef[i])*(a)%(fx->mod_coef);
    }
}


void modulo(POLY fx, POLY gx){   //fx (mod gx)

}
*/
