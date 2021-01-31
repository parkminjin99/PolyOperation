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
    for (int i = 0; i < poly_deg+1; i++)
    {
        for(int j = coef_deg; j >= 0; j--)
        {
            if(a[i][j] != 0)
            {
                COEF_POLY_set(&fx->coef[i], a[i], j);
                break;
            }
        }
    }
    fx->max_degree = 0;
    for (int i = MAX_POLY_DEGREE; i >= 0; i--)
    {
        if(COEF_is_zero(&fx->coef[i]) != TRUE)
        {
            fx->max_degree = i;
            break;
        }
    }
    
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
    for (int i = 0; i < t+1; i++)
    {
        if(gx[i] == 1)
        {
            ctx->mod_gx.coef[i].coef[0] = 1;
            ctx->mod_gx.coef[i].coef_max_degree = 0;
        }
    }
}

void ctx_set_m_12(CTX* ctx, int ft[m+1], int gx[][m+1],int deg_ft, int deg_gx)
{
    COEF_POLY_set(&ctx->mod_coef, ft, deg_ft);
    ctx->mod_gx.max_degree = deg_gx;
    for (int i = 0; i < t+1; i++)
    {
        for (int j = m; j >= 0; j--)
        {
            if(gx[i][j] != 0)
            {
                COEF_POLY_set(&ctx->mod_gx.coef[i], gx[i], j);
                break;
            }
        }
    }
}

void set_COEF_POLY_degree(COEF_POLY* dst)
{
    dst->coef_max_degree = 0;
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
    dst->max_degree=0;
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

int COEF_POLY_equal(COEF_POLY* src1, COEF_POLY* src2)
{
    if(src1->coef_max_degree != src2->coef_max_degree)
    {
        printf("%d %d\n", src1->coef_max_degree, src2->coef_max_degree);
        return FALSE;
    }
    for (int i = 0; i <= src1->coef_max_degree; i++)
    {
        if(src1->coef[i] != src2->coef[i])  
            return FALSE;
    }
    return TRUE;
    
}

int POLY_equal(POLY* src1, POLY* src2)
{
    if(src1->max_degree != src2->max_degree)
    {
        return FALSE;
    }
    for(int i = 0; i <= src1->max_degree; i++)
    {
        if(COEF_POLY_equal(&src1->coef[i], &src2->coef[i]) == FALSE)
        {
            printf("i = %d\n", i);
            return FALSE;
        }
    }
    return TRUE;
}

int COEF_is_one(COEF_POLY* dst)
{
    if(dst->coef_max_degree == 0 && dst->coef[0] == 1)   return TRUE;
    else    return FALSE;
}

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
    int first = 1;                           // 0�� �ƴ� ���� ó������ �ƴ��� �Ǻ��ϱ� ���� ����
    int max_deg = ft->coef_max_degree;       // max_deg ����
    while(max_deg >= 0){                     
        if(ft->coef[max_deg] != 0){          // �ش� ���� 0�� �ƴ϶�� ����
            if(first == 0){                  // 0�� �ƴ� ���� ó���� �ƴ϶�� '+'��ȣ�� ���� print��
                printf(" + ");
            }
            
            if(max_deg == 0){                // 0���� ���� ����̹Ƿ� ����� ���
                printf("%d", ft->coef[max_deg]);
            }
            else if(max_deg == 1)
            {           // 1���� ���� () * t�� ���
                if(ft->coef[max_deg] == 1)      printf("t");
                else                            printf("%d * t", ft->coef[max_deg]);
            }
            else{                            // n���� ���� () * t^n�� ��� (n > 1)
                if(ft->coef[max_deg])           printf("t^%d", max_deg);
                else                            printf("%d * t^%d", ft->coef[max_deg], max_deg);
            }
            first = 0;                       // �� ���̶� if���� ����Ǹ� first�� 0���� ����
        }
        else{
            if(ft->coef_max_degree == 0)     // ft�� 0�̶�� �״�� 0�� ���
                printf("%d", ft->coef[max_deg]);
        }
        max_deg -= 1;                      
    }
}

int COEF_is_zero(COEF_POLY* ft) // ft�� 0���� �ƴ��� �Ǻ��ϴ� �Լ�, 0�̶�� '0' ���, �ƴ϶�� '1'���
{             
    int zero = 0;                            // ft�� 0�̶�� ���������� ������ ����
    int max_deg = ft->coef_max_degree;
    while(max_deg >= 0){                     
        if(ft->coef[max_deg] != 0){          // 0�� �ƴѰ�� zero�� 1�� �߰��ϰ� �ٷ� break
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
    int first = 1;                                           // 0�� �ƴ� ���� ó������ �ƴ��� �Ǻ��ϱ� ���� ����
    int max_deg = fx->max_degree;
    while(max_deg >= 0)
    {
        //printf("\n %d \n", COEF_is_zero(&fx->coef[max_deg]));
        if(COEF_is_zero(&fx->coef[max_deg]) != 0)
        {           // �ش� ���� 0�� �ƴ϶�� ����
            if(first == 0){                                  // 0�� �ƴ� ���� ó���� �ƴ϶�� '+'��ȣ�� ���� print��
                printf(" + ");
            }

            if(max_deg == 0){                                // 0���� ���� ����̹Ƿ� ����� ���
                //printf("(");
                COEF_POLY_print(&fx->coef[max_deg]);
                //printf(")");
            }
            else if(max_deg == 1){                            // 1���� ���� () * X�� ���
                if(fx->coef[max_deg].coef_max_degree == 0 && fx->coef[max_deg].coef[0] == 1)
                    printf("X");
                else
                {
                    printf("(");
                    COEF_POLY_print(&fx->coef[max_deg]);
                    printf(") * X");
                }
            }
            else{                                             // n���� ���� () * X^n�� ��� (n > 1)
                if(fx->coef[max_deg].coef_max_degree == 0 && fx->coef[max_deg].coef[0] == 1)
                    printf("X^%d", max_deg);
                else
                {
                    printf("(");
                    COEF_POLY_print(&fx->coef[max_deg]);
                    printf(") * X^%d", max_deg);   
                }
            }
            first = 0;                                        // �� ���̶� if���� ����Ǹ� first�� 0���� ����
        }
        else
        {
            if(fx->max_degree == 0)                          // ft�� 0�̶�� �״�� 0�� ���
                COEF_POLY_print(&fx->coef[max_deg]);
        }
        max_deg -= 1;
    }
    printf("\n");
}

void int2vec(OUT int* vec, OUT int* vec_size, IN int zz) 
{ // ������ ���ͷ�
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
{ // ���͸� ������
    *zz = 0;
    for (int i = 0; i < vec_size; i++)
        *zz += vec[i] * pow(2, i);
}

void coef_modft_table(OUT COEF_POLY* ft_table, IN CTX* ctx)
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
            COEF_POLY_set(&ft_table[i],vec,12);
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
            COEF_POLY_set(&ft_table[i],vec,13);
            for(int j=m-1;j>=0;j--)
                vec[j+1]=vec[j];
            vec[0]=0;
        }
    }
}

void coef_squ(OUT int* asqu,IN COEF_POLY* ft_table, IN COEF_POLY* a, IN CTX* ctx)
{  
    COEF_POLY tmp;
    COEF_POLY_init(&tmp,0);
    
    for(int i=a->coef_max_degree;i>=0;i--)
        tmp.coef[2*i]=a->coef[i];

    tmp.coef_max_degree=2*a->coef_max_degree;
    for(int i=tmp.coef_max_degree;i>=m;i=i-2)
    {
        if(tmp.coef[i]==1)
        {
            tmp.coef[i]=0;
            for(int j=0;j<=ft_table[i-m].coef_max_degree;j++)
                tmp.coef[j]^=ft_table[i-m].coef[j];
        }
    }
    vec2int(asqu, tmp.coef, m+1);
}

void gen_Ttable(OUT int* Ttable, OUT int* InvTtable, IN COEF_POLY* ft_table, IN CTX* ctx )
{  // f2m�� ��� ����.
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

void gen_Xitable(OUT POLY Xtable[], IN CTX* ctx, IN COEF_POLY ft_table[])
{ 
    POLY_copy(&Xtable[0],&ctx->mod_gx);
    COEF_POLY_init(&Xtable[0].coef[t], 0);
    set_POLY_degree(&Xtable[0]);
    POLY temp;  
    for(int i = 1; i <= t; i++) // Xtable[i] POLY ���� 
    {
        POLY_init(&Xtable[i], 0);
        for (int j = 1; j <= t-1; j++)
        {
            COEF_POLY_copy(&Xtable[i].coef[j], &Xtable[i-1].coef[j-1]);
            if(COEF_is_zero(&Xtable[i].coef[j]) != TRUE)
                Xtable[i].max_degree = j;
        }
        if(COEF_is_zero(&Xtable[i-1].coef[t-1]) != TRUE)
        {
            POLY_init(&temp, 0);
            MULscalar(&temp, &Xtable[0], &Xtable[i-1].coef[t-1], ft_table, ctx);
            POLY_add_zzx(&Xtable[i], &temp);
        }
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

void POLY_mod_gx(OUT IN POLY* dst, IN CTX* ctx, IN POLY Xtable[], IN COEF_POLY ft_table[])
{  
    if(dst->max_degree < t)    return;
    int max_deg = dst->max_degree;
    POLY temp;
    for(int i = t; i <= max_deg; i++) 
    {
        if(COEF_is_zero(&dst->coef[i])==FALSE)
        {
            POLY_init(&temp, 0);
            MULscalar(&temp, &Xtable[i-t], &dst->coef[i], ft_table, ctx);
            POLY_add_zzx(dst, &temp);
            COEF_POLY_init(&dst->coef[i], 0);
        }
        // for(int j = 0; j <= Xtable[i-t].max_degree; j++) // mul scalar ��ü �κ� ���߿� �����Ǹ� ���� ����
        // {
        //     if(Xtable[i-t].coef[j].coef[0] != 0)
        //         COEF_POLY_add_zzx(&dst->coef[j], &dst->coef[i]);
        // }
        // COEF_POLY_init(&dst->coef[i], 0);
    }
    set_POLY_degree(dst); // dst poly ���� ���ϴ� �κ�
}

void X_sqrt(OUT POLY* x_sqrt, IN POLY Xtable[], IN POLY* src,IN COEF_POLY* fttable, IN int* Ttable, IN CTX* ctx)
{  //x^i
    POLY tmp,aaa;
    POLY x;
    int k,n;
    int vec[MAX_COEF_POLY_DEGREE]={0,}, vec_size=0;

    POLY_init(&x,0);
    POLY_copy(&x,src);
    //POLY_print(&x);
    for(int count=0;count<m*t-1;count++)
    {
        POLY_init(&tmp,0);
        for(int i=x.max_degree;i>=0;i--)
        {
            for(int j=0;j<x.coef[i].coef_max_degree+1;j++)
                vec[j]=x.coef[i].coef[j];
            vec2int(&k, vec,x.coef[i].coef_max_degree+1);
            n = Ttable[k];
            int2vec(vec,&vec_size,n);
            COEF_POLY_set(&tmp.coef[2*i],vec,vec_size-1);
        }
        // tmp.max_degree = 0;
        // for(int i = MAX_POLY_DEGREE; i >= 0; i--)
        // {
        //     if(COEF_is_zero(&tmp.coef[i]) != TRUE)
        //     {
        //         tmp.max_degree = i;
        //         break;
        //     }
        // }
        //printf("tmp = ");   POLY_print(&tmp);
        tmp.max_degree = 2 * x.max_degree;
        POLY_mod_gx(&tmp,ctx,Xtable, fttable);
        // tmp.max_degree = 2 * x.max_degree;
        // for(int i=tmp.max_degree;i>=t;i=i-2)
        // {
        //     POLY_init(&aaa,0);
        //     if(COEF_is_zero(&tmp.coef[i])==FALSE)
        //     {
        //         MULscalar(&aaa,&Xtable[i-t],&tmp.coef[i],fttable,ctx);
        //         //printf("xt="); POLY_print(&Xtable[i-t]); 
        //         //printf("  %d  ",aaa.max_degree);
        //         for(int j=0;j<tmp.coef[i].coef_max_degree+1;j++)
        //         {
        //             tmp.coef[i].coef[j]=0; // �ϴ� i 0����
        //         }
        //         //POLY_print(&aaa); printf("\n");
        //         POLY_add_zzx(&tmp,&aaa); 
        //     }
        // //printf("a^2 mod g = "); COEF_POLY_print(&tmp); printf("\n");
        // }
       //POLY_print(&tmp); printf("\n");
        POLY_copy(&x,&tmp);
    }
    POLY_copy(x_sqrt,&x);

    //printf("int (a^2 mod g) = %d\n",*asqu);
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
    for(int i = 0; i <= src1->coef_max_degree; i++)
        for(int j = 0; j <= src2->coef_max_degree; j++)
            dst->coef[i+j] ^= src1->coef[i] * src2->coef[j];
    dst->coef_max_degree = 0;
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
    dst->max_degree = 0;
    for(int i = MAX_POLY_DEGREE; i >= 0; i--)
    {
        if(COEF_is_zero(&dst->coef[i]) != TRUE)    
        {
            dst->max_degree = i;
            break;
        }
    }
    POLY_mod_gx(dst, ctx, Xtable, fttable);
}   

void MULscalar_zzx(IN OUT POLY* dst, IN COEF_POLY* src, IN COEF_POLY ft_table[], IN CTX* ctx)
{
    for(int i = 0; i <= dst->max_degree; i++)
        COEF_POLY_mul_zzx(&dst->coef[i], src, ft_table, ctx);
    set_POLY_degree(dst);
}

void MULscalar(OUT POLY* dst, IN POLY* src1, IN COEF_POLY* src2, IN COEF_POLY ft_table[], IN CTX* ctx)
{
    if(COEF_is_one(src2) == TRUE)
    {
        POLY_copy(dst, src1);
        return;
    }
    dst->max_degree = src1->max_degree;
    for(int i = 0; i <= src1->max_degree; i++)
        COEF_POLY_mul(&dst->coef[i], &src1->coef[i], src2, ft_table, ctx);
    set_POLY_degree(dst); // dst�� ������ �޶��� �� �����Ƿ� 
}

void COEF_POLY_add_zzx(OUT COEF_POLY* dst, IN COEF_POLY* src)
{
    int max_deg = MAX(src->coef_max_degree, dst->coef_max_degree);
    dst->coef_max_degree = 0;
    for(int i = 0; i <= max_deg; i++)
    {
        dst->coef[i] ^= src->coef[i];
        if(dst->coef[i] != 0)   dst->coef_max_degree = i;
    }
}

void POLY_add_zzx(OUT POLY* dst, IN POLY* src) // t�� �̻��� ���׽��� �ȵ��´ٴ� �����Ͽ� ���� 
{
    for (int i = 0; i < t; i++)
    {
        COEF_POLY_add_zzx(&dst->coef[i], &src->coef[i]);
        if(COEF_is_zero(&dst->coef[i]) != TRUE)     dst->max_degree = i;  
    }
}

void COEF_POLY_add(OUT COEF_POLY* dst, IN COEF_POLY* src1, IN COEF_POLY* src2)
{
    if(COEF_is_zero(src1) == TRUE && COEF_is_zero(src2) == TRUE)
    {
        COEF_POLY_init(dst, 0);
        return;
    }
    if(COEF_is_zero(src1) == TRUE)
    {
        COEF_POLY_copy(dst, src2);
        return;
    }
    if(COEF_is_zero(src2) == TRUE)
    {
        COEF_POLY_copy(dst, src1);
        return;
    }
    int max_deg = MAX(src1->coef_max_degree, src2->coef_max_degree);
    dst->coef_max_degree = 0;
    for(int i = 0; i <= max_deg; i++)
    {
        dst->coef[i] = src1->coef[i] ^ src2->coef[i];
        if(dst->coef[i] != 0)   dst->coef_max_degree = i;
    }
}

void POLY_add(OUT POLY* dst, IN POLY* src1, IN POLY* src2)
{   
    int max_deg = MAX(src1->max_degree, src2->max_degree);
    dst->max_degree = 0;
    for (int i = 0; i <= max_deg; i++)
    {
        COEF_POLY_add(&dst->coef[i], &src1->coef[i], &src2->coef[i]);
        if(COEF_is_zero(&dst->coef[i]) != TRUE)     dst->max_degree = i;  
    }
}