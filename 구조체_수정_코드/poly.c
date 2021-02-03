#include "poly.h"
#define MAX(x,y)( (x)>(y)?(x):(y) ) 
#define get_j_th_bit(x,j) ((x>>j)&0x1)

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

void POLY_init(POLY* fx, int maxdeg)
{ 
    fx->max_degree = maxdeg;
    memset(fx->coef, 0, (MAX_POLY_DEGREE+1)*sizeof(int));
}
void POLY_set(POLY* fx, int a[], int poly_deg)
{  
    fx->max_degree = 0;
    for (int i = 0; i < poly_deg+1; i++)
    {   
        fx->coef[i] = a[i];
        if(fx->coef[i] != 0)    fx->max_degree = i;
    }
}

void ctx_init(CTX* ctx)
{
    POLY_init(&ctx->mod_gx,0);
    ctx->mod_coef = 0;
}

void ctx_set(CTX* ctx, int ft, int gx[t+1], int deg_gx)
{ 
    ctx->mod_coef = ft;
    
    ctx->mod_gx.max_degree = deg_gx;
    for (int i = 0; i < t+1; i++)
    {
        if(gx[i] != 0)
            ctx->mod_gx.coef[i] = gx[i];
    }
}

void set_POLY_degree(POLY* dst)
{
    dst->max_degree=0;
    for(int i = t-1; i >= 0; i--)
    {    
        if(dst->coef[i] != 0)    
        {
            dst->max_degree = i;
            break;
        }
    }
}

//===================================================================================

int POLY_equal(POLY* src1, POLY* src2)
{
    if(src1->max_degree != src2->max_degree)
        return FALSE;
    for(int i = 0; i <= src1->max_degree; i++)
    {
        if(src1->coef[i] != src2->coef[i])
        {
            printf("i = %d\n", i);
            return FALSE;
        }
    }
    return TRUE;
}

void POLY_copy(OUT POLY* dst, IN POLY* src)
{
    dst->max_degree = src->max_degree;
    for (int i = 0; i < dst->max_degree+1; i++)
        dst->coef[i] = src->coef[i];
}

void COEF_POLY_print(int ft)
{
    int first = 1;                           // 0이 아닌 수가 처음인지 아닌지 판별하기 위한 변수
    int max_deg = MAX_COEF_POLY_DEGREE;       // max_deg 설정
    int print_2_m[MAX_COEF_POLY_DEGREE+1] = {1,};
    for (int i = 1; i < MAX_COEF_POLY_DEGREE+1; i++)
        print_2_m[i] = print_2_m[i-1] << 1;
    
    while(max_deg >= 0)
    {              
        int temp = ft&print_2_m[max_deg];
        if(temp == print_2_m[max_deg])
        {          // 해당 수가 0이 아니라면 실행
            if(first == 0){                  // 0이 아닌 수가 처음이 아니라면 '+'기호를 먼저 print함
                printf(" + ");
            }
            
            if(max_deg == 0)
            {                // 0차의 수는 상수이므로 상수만 출력
                if(temp == print_2_m[max_deg])                   printf("1");
                else                                             printf("0");
            }
            else if(max_deg == 1)
            {           // 1차의 수는 () * t로 출력
                if(temp == print_2_m[max_deg])                   printf("t");
            }
            else
            {                            // n차의 수는 () * t^n로 출력 (n > 1)
                if(temp == print_2_m[max_deg])                   printf("t^%d", max_deg);
            }
            first = 0;                       // 한 번이라도 if문이 실행되면 first를 0으로 변경
        }
        else
        {
            if(ft == 0)     // ft가 0이라면 그대로 0만 출력
                printf("0");
        }
        max_deg -= 1;                      
    }
}

int POLY_is_zero(POLY* fx)
{
    if(fx->max_degree == 0 && fx->coef[0] == 0)   return TRUE;
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
        if(fx->coef[max_deg] != 0)
        {           // 해당 수가 0이 아니라면 실행
            if(first == 0){                                  // 0이 아닌 수가 처음이 아니라면 '+'기호를 먼저 print함
                printf(" + ");
            }

            if(max_deg == 0)
            {                                // 0차의 수는 상수이므로 상수만 출력
                //printf("(");
                COEF_POLY_print(fx->coef[max_deg]);
                //printf(")");
            }
            else if(max_deg == 1)
            {                            // 1차의 수는 () * X로 출력
                if(fx->coef[max_deg] == 1)
                    printf("X");
                else
                {
                    printf("(");
                    COEF_POLY_print(fx->coef[max_deg]);
                    printf(") * X");
                }
            }
            else{                                             // n차의 수는 () * X^n로 출력 (n > 1)
                if(fx->coef[max_deg] == 1)
                    printf("X^%d", max_deg);
                else
                {
                    printf("(");
                    COEF_POLY_print(fx->coef[max_deg]);
                    printf(") * X^%d", max_deg);   
                }
            }
            first = 0;                                        // 한 번이라도 if문이 실행되면 first를 0으로 변경
        }
        else
        {
            if(fx->max_degree == 0)                          // ft가 0이라면 그대로 0만 출력
                COEF_POLY_print(fx->coef[max_deg]);
        }
        max_deg -= 1;
    }
    printf("\n");
}

void int2vec(OUT int* vec, OUT int* vec_size, IN int zz) 
{ // 정수를 벡터로
    int count = 0;
    for(int i=0;i< *vec_size;i++)
        vec[i]=0;
    while (1) 
    {
        vec[count] = zz % 2;
        zz = zz / 2;
        if (zz == 0)
            break;
        count++;
    }
    *vec_size = count+1;
}

void vec2int(OUT int* zz, IN int* vec, IN int vec_size) 
{ // 벡터를 정수로
    *zz = 0;
    for (int i = 0; i < vec_size; i++)
        *zz += vec[i] * pow(2, i);
}

void coef_modft_table(OUT int ft_table[], IN CTX* ctx)
{  
    if(m==12)  // (1,1,0,1,0,1,1,1,0,0,0,0,1)   
    {
        int vec[m+1]={0,0,0,0,0,0,0,0,0,0,0,0,1};
        for(int i = 0; i <= m; i++)  // t^12 ~ t^22
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
            vec2int(&ft_table[i], vec, 12); //COEF_POLY_set(&ft_table[i],vec,12);
            for(int j=m-1;j>=0;j--)
                vec[j+1]=vec[j];
            vec[0]=0;
        }
    }
    else  // m=13 (1,1,0,1,1,0,0,0,0,0,0,0,0,1)
    {
        int vec[m+1]={0,0,0,0,0,0,0,0,0,0,0,0,0,1};
        for(int i = 0; i <= m; i++)  // t^13 ~ t^24
        {            
            if(vec[13]==1)
            {
                vec[13]=0;
                vec[0]^=1;
                vec[1]^=1;
                vec[3]^=1;
                vec[4]^=1;
            }
            vec2int(&ft_table[i], vec, 13); //COEF_POLY_set(&ft_table[i],vec,13);
            for(int j=m-1;j>=0;j--)
                vec[j+1]=vec[j];
            vec[0]=0;
        }
    }
}

void coef_squ(OUT int* asqu,IN int ft_table[], IN int a[], IN int a_size, IN CTX* ctx)
{  
    int vec[MAX_COEF_POLY_DEGREE]={0,}, vec_temp[MAX_COEF_POLY_DEGREE] = {0,}, vec_temp_size = 0;
    for(int i = a_size-1 ; i>=0; i--)
        vec[2*i] = a[i];
    for(int i = 2*a_size; i>=m; i=i-2)
    {
        if(vec[i]==1)
        {
            vec[i]=0;
            int2vec(vec_temp, &vec_temp_size, ft_table[i-m]);
            for(int j=0; j < vec_temp_size;j++)
                vec[j]^=vec_temp[j];
        }
    }
    vec2int(asqu, vec, m+1);
}

void gen_Ttable(OUT int* Ttable, OUT int* InvTtable, IN int ft_table[], IN CTX* ctx)
{  // f2m의 모든 원소.
    int vec[MAX_COEF_POLY_DEGREE]={0,}, vec_size=0, dst=0;
    for(int i=0;i<pow(2,m);i++)
    {
        int2vec(vec, &vec_size, i);
        coef_squ(&dst,ft_table,vec, vec_size, ctx);
        Ttable[i]=dst;
        InvTtable[dst]=i;
    }
}

void gen_Xitable(OUT POLY Xtable[], IN CTX* ctx, IN int ft_table[])
{ 
    POLY_copy(&Xtable[0],&ctx->mod_gx);
    Xtable[0].coef[t] = 0;      //COEF_POLY_init(&Xtable[0].coef[t], 0);
    set_POLY_degree(&Xtable[0]);
    POLY temp;  
    for(int i = 1; i <= t; i++) // Xtable[i] POLY 생성 
    {
        POLY_init(&Xtable[i], 0);
        for (int j = 1; j <= t-1; j++)
        {
            Xtable[i].coef[j] = Xtable[i-1].coef[j-1];         //COEF_POLY_copy(&Xtable[i].coef[j], &Xtable[i-1].coef[j-1]);
            if(Xtable[i].coef[j] != 0)
                Xtable[i].max_degree = j;
        }
        if(Xtable[i-1].coef[t-1] != 0)
        {
            POLY_init(&temp, 0);
            MULscalar(&temp, &Xtable[0], Xtable[i-1].coef[t-1], ft_table, ctx);
            POLY_add_zzx(&Xtable[i], &temp);
        }
    }
}


void COEF_POLY_mod_ft(OUT IN int* dst, IN CTX* ctx, IN int fttable[])
{
    for(int i = m; i <= MAX_COEF_POLY_DEGREE; i++)
    {
        if(get_j_th_bit(*dst, i) == 1)
        {
            *dst = *dst & (~(1<<i));
            *dst ^= fttable[i-m];   //COEF_POLY_add_zzx(dst, &fttable[i-m]);
        }
    }
}

void POLY_mod_gx(OUT IN POLY* dst, IN CTX* ctx, IN POLY Xtable[], IN int ft_table[])
{  
    if(dst->max_degree < t)    return;
    int max_deg = dst->max_degree;
    POLY temp;
    for(int i = t; i <= max_deg; i++) 
    {
        if(dst->coef[i] != 0)
        {
            POLY_init(&temp, 0);
            MULscalar(&temp, &Xtable[i-t], dst->coef[i], ft_table, ctx);
            POLY_add_zzx(dst, &temp);
            dst->coef[i] = 0;
        }
        // for(int j = 0; j <= Xtable[i-t].max_degree; j++) // mul scalar 대체 부분 나중에 구현되면 적용 예정
        // {
        //     if(Xtable[i-t].coef[j].coef[0] != 0)
        //         COEF_POLY_add_zzx(&dst->coef[j], &dst->coef[i]);
        // }
        // COEF_POLY_init(&dst->coef[i], 0);
    }
    set_POLY_degree(dst); // dst poly 차수 정하는 부분
}

void X_sqrt(OUT POLY* x_sqrt, IN POLY Xtable[], IN POLY* src, IN int fttable[], IN int* Ttable, IN CTX* ctx)
{  //x^i
    POLY tmp;
    POLY_init(x_sqrt,0);
    POLY_copy(x_sqrt,src);
    for(int count=0;count<m*t-1;count++)
    {
        POLY_init(&tmp,0);
        for(int i=x_sqrt->max_degree;i>=0;i--)
            tmp.coef[2*i] = Ttable[x_sqrt->coef[i]];
        tmp.max_degree = 2 * x_sqrt->max_degree;
        POLY_mod_gx(&tmp,ctx,Xtable, fttable);
        POLY_copy(x_sqrt,&tmp);
    }
}


void COEF_POLY_mul_zzx(OUT int* dst, IN int src, IN int ft_table[], IN CTX* ctx)
{
    int temp = *dst;
    *dst = 0;
    COEF_POLY_mul(dst, temp, src, ft_table, ctx);
}

void COEF_POLY_mul(OUT int* dst, IN int src1, IN int src2, IN int ft_table[], IN CTX* ctx) 
{
    for(int i = 0; i <= MAX_COEF_POLY_DEGREE; i++)
    {
        if(get_j_th_bit(src1,i) == 0)   continue;
        else    *dst ^= src2<<i;
    }
    COEF_POLY_mod_ft(dst, ctx, ft_table);
}

void POLY_mul(OUT POLY* dst, IN POLY* src1, IN POLY* src2, IN CTX* ctx, IN int fttable[], IN POLY Xtable[])
{
    for(int i = 0;i <= src1->max_degree; i++)
        for(int j = 0; j <= src2->max_degree; j++)
            COEF_POLY_mul(&dst->coef[i+j], src1->coef[i], src2->coef[j], fttable, ctx);
    dst->max_degree = 0;
    for(int i = MAX_POLY_DEGREE; i >= 0; i--)
    {
        if(dst->coef[i] != 0)    
        {
            dst->max_degree = i;
            break;
        }
    }
    POLY_mod_gx(dst, ctx, Xtable, fttable);
}   

void MULscalar_zzx(IN OUT POLY* dst, IN int src, IN int ft_table[], IN CTX* ctx)
{
    for(int i = 0; i <= dst->max_degree; i++)
        COEF_POLY_mul_zzx(&dst->coef[i], src, ft_table, ctx);
    set_POLY_degree(dst);
}

void MULscalar(OUT POLY* dst, IN POLY* src1, IN int src2, IN int ft_table[], IN CTX* ctx)
{
    if(src2 == 1)
    {
        POLY_copy(dst, src1);
        return;
    }
    dst->max_degree = src1->max_degree;
    for(int i = 0; i <= src1->max_degree; i++)
        COEF_POLY_mul(&dst->coef[i], src1->coef[i], src2, ft_table, ctx);
    set_POLY_degree(dst); // dst의 차수가 달라질 수 있으므로 
}


void POLY_add_zzx(OUT POLY* dst, IN POLY* src) // t차 이상의 다항식은 안들어온다는 가정하에 덧셈 
{
    for (int i = 0; i < t; i++)
    {
        dst->coef[i] ^= src->coef[i];   //COEF_POLY_add_zzx(&dst->coef[i], &src->coef[i]);
        if(dst->coef[i] != 0)     dst->max_degree = i;  
    }
}

void POLY_add(OUT POLY* dst, IN POLY* src1, IN POLY* src2)
{   
    int max_deg = MAX(src1->max_degree, src2->max_degree);
    dst->max_degree = 0;
    for (int i = 0; i <= max_deg; i++)
    {
        dst->coef[i] = src1->coef[i] ^ src2->coef[i]; // COEF_POLY_add(&dst->coef[i], &src1->coef[i], &src2->coef[i]);
        if(dst->coef[i] != 0)     dst->max_degree = i;  
    }
}
