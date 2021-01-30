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

void set_COEF_POLY_degree(COEF_POLY* dst)
{
    for(int i = m-1 ;i >= 0; i--)
    {
        if(dst->coef[i]!=0)
        {
            dst->coef_max_degree = i;
            break;
        }
    }
}

void set_POLY_degree(POLY* dst)
{
    for(int i = t-1; i >= 0; i--)
    {    
        if(COEF_is_zero(&dst->coef[i]) != TRUE)    
        {
            dst->max_degree = i;
            break;
        }
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

void COEF_POLY_print(COEF_POLY* ft)
{
    int first = 1;                           // 0이 아닌 수가 처음인지 아닌지 판별하기 위한 변수
    int max_deg = ft->coef_max_degree;       // max_deg 설정
    while(max_deg >= 0){                     
        if(ft->coef[max_deg] != 0){          // 해당 수가 0이 아니라면 실행
            if(first == 0){                  // 0이 아닌 수가 처음이 아니라면 '+'기호를 먼저 print함
                printf(" + ");
            }
            
            if(max_deg == 0){                // 0차의 수는 상수이므로 상수만 출력
                printf("%d", ft->coef[max_deg]);
            }
            else if(max_deg == 1)
            {           // 1차의 수는 () * t로 출력
                if(ft->coef[max_deg] == 1)      printf("t");
                else                            printf("%d * t", ft->coef[max_deg]);
            }
            else{                            // n차의 수는 () * t^n로 출력 (n > 1)
                if(ft->coef[max_deg])           printf("t^%d", max_deg);
                else                            printf("%d * t^%d", ft->coef[max_deg], max_deg);
            }
            first = 0;                       // 한 번이라도 if문이 실행되면 first를 0으로 변경
        }
        else{
            if(ft->coef_max_degree == 0)     // ft가 0이라면 그대로 0만 출력
                printf("%d", ft->coef[max_deg]);
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
                //printf("(");
                COEF_POLY_print(&fx->coef[max_deg]);
                //printf(")");
            }
            else if(max_deg == 1){                            // 1차의 수는 () * X로 출력
                if(fx->coef[max_deg].coef_max_degree == 0 && fx->coef[max_deg].coef[0] == 1)
                    printf("X");
                else
                {
                    printf("(");
                    COEF_POLY_print(&fx->coef[max_deg]);
                    printf(") * X");
                }
            }
            else{                                             // n차의 수는 () * X^n로 출력 (n > 1)
                if(fx->coef[max_deg].coef_max_degree == 0 && fx->coef[max_deg].coef[0] == 1)
                    printf("X^%d", max_deg);
                else
                {
                    printf("(");
                    COEF_POLY_print(&fx->coef[max_deg]);
                    printf(") * X^%d", max_deg);   
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

void coef_modft_table(OUT COEF_POLY* ft_table, IN CTX* ctx){  
    
    if(m==12)  // (1,1,0,1,0,1,1,1,0,0,0,0,1)   
    {
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
            
            COEF_POLY_set(&ft_table[i],vec,12);
            //vec2int(&ft_table[i],vec,13);
            
            for(int j=m-1;j>=0;j--)
            {
                vec[j+1]=vec[j];
            }
            vec[0]=0;
        }
    }
    else  // m=13 (1,1,0,1,1,0,0,0,0,0,0,0,0,1)
    {
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
            
            //vec2int(&ft_table[i],vec,14);
            COEF_POLY_set(&ft_table[i],vec,13);

            //printf("%d ",ft_table[i]);
            //for(int j=0;j<m+1;j++)
            //    printf("%d",vec[j]);
            //printf("\n");  
            for(int j=m-1;j>=0;j--)
            {
                vec[j+1]=vec[j];
            }
            vec[0]=0;
        }
    }
}

void coef_squ(OUT int* asqu,IN COEF_POLY* ft_table, IN COEF_POLY* a, IN CTX* ctx)
{  
    //int vec[MAX_COEF_POLY_DEGREE]={0,}, vec_size=0;
    COEF_POLY tmp;
    COEF_POLY_init(&tmp,0);
    //printf("\na = "); COEF_POLY_print(a); printf("\n");
    
    for(int i=a->coef_max_degree;i>=0;i--)
        tmp.coef[2*i]=a->coef[i];

    tmp.coef_max_degree=2*a->coef_max_degree;

    //printf("a^2 = "); COEF_POLY_print(&tmp); printf("  ");
    //printf("maxdeg= %d\n",tmp.coef_max_degree);

    for(int i=tmp.coef_max_degree;i>=m;i=i-2)
    {
        if(tmp.coef[i]==1)
        {
            tmp.coef[i]=0;

            //int2vec(vec,&vec_size,ft_table[i-m]);
            //printf(" %d ",ft_table[i-m]);
            for(int j=0;j<=ft_table[i-m].coef_max_degree;j++)
            {
                tmp.coef[j]^=ft_table[i-m].coef[j];
            }
        }
        //printf("a^2 mod g = "); COEF_POLY_print(&tmp); printf("\n");
        
    }
    //COEF_POLY_print(&tmp); printf("\n");
    vec2int(asqu, tmp.coef, m+1);
    //printf("int (a^2 mod g) = %d\n",*asqu);

}

void gen_Ttable(OUT int* Ttable, OUT int* InvTtable, IN COEF_POLY* ft_table, IN CTX* ctx ){  // f2m의 모든 원소.
    int vec[MAX_COEF_POLY_DEGREE]={0,}, vec_size=0, dst=0;
    COEF_POLY tmp;
    for(int i=0;i<pow(2,m);i++)
    {
        int2vec(vec,&vec_size,i);
        COEF_POLY_set(&tmp, vec, vec_size-1);
        coef_squ(&dst,ft_table,&tmp,ctx);
        Ttable[i]=dst;
        InvTtable[dst]=i;
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
            POLY_add_zzx(&Xtable[i], &Xtable[0]);
    }
}

void COEF_POLY_mod_ft(OUT IN COEF_POLY* dst, IN CTX* ctx, IN COEF_POLY fttable[])
{
    if(dst->coef_max_degree < m)    return;
    int max_deg = dst->coef_max_degree;
    for(int i = m; i<= max_deg; i++)
    {
        if(dst->coef[i] == 1)
        {
            dst->coef[i] = 0;
            COEF_POLY_add_zzx(dst, &fttable[i-m]);
        }
    }
    set_COEF_POLY_degree(dst);
}

void POLY_mod_gx(OUT IN POLY* dst, IN CTX* ctx, IN POLY Xtable[])
{  
    if(dst->max_degree < t)    return;
    int max_deg = dst->max_degree;
    for(int i = t; i <= max_deg; i++) 
    {
        for(int j = 0; j <= Xtable[i-t].max_degree; j++) // mul scalar 대체 부분 나중에 구현되면 적용 예정
            if(Xtable[i-t].coef[j].coef[0] != 0)
                COEF_POLY_add_zzx(&dst->coef[j], &dst->coef[i]);
        COEF_POLY_init(&dst->coef[i], 0);
        //POLY_print(dst);
    }
    set_POLY_degree(dst); // dst poly 차수 정하는 부분
}

void X_sqrt(OUT POLY* x_sqrt, IN POLY* x, IN CTX* ctx)
{  

}

void COEF_POLY_mul_zzx(OUT COEF_POLY* dst, IN COEF_POLY* src, IN COEF_POLY* ft_table, IN CTX* ctx)
{
    COEF_POLY temp;
    COEF_POLY_init(&temp, 0);
    COEF_POLY_copy(&temp, dst);
    COEF_POLY_init(dst, 0);
    COEF_POLY_mul(dst,&temp,src, ft_table, ctx);
}

void COEF_POLY_mul(OUT COEF_POLY* dst,IN COEF_POLY* src1, IN COEF_POLY* src2, IN COEF_POLY* ft_table, IN CTX* ctx) 
{
    int vec[MAX_COEF_POLY_DEGREE]={0,};
    for(int i = 0; i <= src1->coef_max_degree; i++)
        for(int j = 0; j <= src2->coef_max_degree; j++)
            dst->coef[i+j] ^= src1->coef[i] * src2->coef[j];

    for(int i = MAX_COEF_POLY_DEGREE; i >= 0; i--) 
    {
        if(dst->coef[i] != 0)    
        {
            dst->coef_max_degree = i;
            break;
        }
    }
    COEF_POLY_mod_ft(dst, ctx, ft_table);
}

void POLY_mul(OUT POLY* dst, IN POLY* src1, IN POLY* src2, IN CTX* ctx, IN COEF_POLY fttable[], IN POLY Xtable[])
{
    for(int i = 0;i <= src1->max_degree; i++)
        for(int j = 0; j <= src2->max_degree; j++)
            COEF_POLY_mul(&dst->coef[i+j], &src1->coef[i], &src2->coef[j], fttable, ctx);
    for(int i = MAX_POLY_DEGREE; i >= 0; i--)
    {
        if(COEF_is_zero(&dst->coef[i]) != TRUE)    
        {
            dst->max_degree = i;
            break;
        }
    }
    //POLY_print(dst);
    POLY_mod_gx(dst, ctx, Xtable);
}   

void MULscalar_zzx(IN OUT POLY* dst, IN COEF_POLY* src, IN COEF_POLY ft_table[], IN CTX* ctx)
{
    for(int i = 0; i <= dst->max_degree; i++)
        COEF_POLY_mul_zzx(&dst->coef[i], src, ft_table, ctx);
    set_POLY_degree(dst);
}

void MULscalar(OUT POLY* dst, IN POLY* src1, IN COEF_POLY* src2, IN COEF_POLY ft_table[], IN CTX* ctx)
{
    dst->max_degree = src1->max_degree;
    for(int i = 0; i <= src1->max_degree; i++)
        COEF_POLY_mul(&dst->coef[i], &src1->coef[i], src2, ft_table, ctx);
    set_POLY_degree(dst); // dst의 차수가 달라질 수 있으므로 
}

void COEF_POLY_add_zzx(OUT COEF_POLY* dst, IN COEF_POLY* src)
{
    int max_deg = MAX(src->coef_max_degree, dst->coef_max_degree);
    for(int i = 0; i <= max_deg; i++)
    {
        dst->coef[i] ^= src->coef[i];
        if(dst->coef[i] != 0)   dst->coef_max_degree = i;
    }
}

void POLY_add_zzx(OUT POLY* dst, IN POLY* src) // t차 이상의 다항식은 안들어온다는 가정하에 덧셈 
{
    for (int i = 0; i < t; i++)
    {
        COEF_POLY_add_zzx(&dst->coef[i], &src->coef[i]);
        if(COEF_is_zero(&dst->coef[i]) != TRUE)     dst->max_degree = i;  
    }
}

void COEF_POLY_add(OUT COEF_POLY* dst, IN COEF_POLY* src1, IN COEF_POLY* src2)
{
    int max_deg = MAX(src1->coef_max_degree, src2->coef_max_degree);
    for(int i = 0; i <= max_deg; i++)
    {
        dst->coef[i] = src1->coef[i] ^ src2->coef[i];
        if(dst->coef[i] != 0)   dst->coef_max_degree = i;
    }
}

void POLY_add(OUT POLY* dst, IN POLY* src1, IN POLY* src2)
{   
    int max_deg = MAX(src1->max_degree, src2->max_degree);
    for (int i = 0; i <= max_deg; i++)
    {
        COEF_POLY_add(&dst->coef[i], &src1->coef[i], &src2->coef[i]);
        if(COEF_is_zero(&dst->coef[i]) != TRUE)     dst->max_degree = i;  
    }
}
