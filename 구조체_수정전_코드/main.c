#include "poly.h"
void test_coef_poly()
{
    COEF_POLY A;
    int a[10] = {1,1,1,1,1,1,1,1,1,1};
    COEF_POLY B;
    int b[7] = {0,0,0,0,0,1,0};    
    COEF_POLY C;
    int c[1] = {1};    
    COEF_POLY D;
    int d[1] = {0};
    
    COEF_POLY_init(&A, 0); 
    COEF_POLY_set(&A, a, 9);
    COEF_POLY_init(&B, 0); 
    COEF_POLY_set(&B, b, 6);
    COEF_POLY_init(&C, 0); 
    COEF_POLY_set(&C, c, 0);
    COEF_POLY_init(&D, 0); 
    COEF_POLY_set(&D, d, 0);
 
    printf("\n=========COEF_POLY===========\n");
    printf("\nA == ");
    COEF_POLY_print(&A);
    //zero = COEF_is_zero(&A);
    //printf("\nA == %d", zero);
    
    printf("\nB == ");
    COEF_POLY_print(&B);
    //zero = COEF_is_zero(&B);
    //printf("\nB == %d", zero);
    
    printf("\nC == ");
    COEF_POLY_print(&C);
    //zero = COEF_is_zero(&C);
    //printf("\nC == %d", zero);
    
    printf("\nD == ");
    COEF_POLY_print(&D);
    //zero = COEF_is_zero(&D);
    //printf("\nD == %d", zero);
}

void test_poly()
{
    printf("\n========= POLY set ===========\n");
    POLY AA, AA_copy;
    int aa[3+1][MAX_COEF_POLY_DEGREE+1] = {{1,1,1,1,1,1,}, {1,0,0,0,0,1,}, {1,0,1,0,1,1,}, {0,1,0,1,0,1,}};
    POLY_init(&AA, 0);
    POLY_set(&AA, aa, 3, 5);
    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < MAX_COEF_POLY_DEGREE+1; j++)
            printf("%d ", AA.coef[i].coef[j]);
        printf("\n");
    }
    POLY_print(&AA);

    printf("\n========= POLY copy ===========\n");
    POLY_init(&AA_copy, 0);
    POLY_copy(&AA_copy, &AA);
    POLY_print(&AA_copy);
}

void set_CTX(CTX* ctx)
{
    printf("\n========= Set CTX ===========\n");
#if m == 12
    int ft[m+1]={1,1,0,1,0,1,1,1,0,0,0,0,1};
#elif m == 13
    int ft[m+1]={1,1,0,1,1,0,0,0,0,0,0,0,0,1};
#endif
    ctx_init(ctx);
#if t == 64 // 64託ぞぞぞぞぞぞぞぞ 鯵旭精暗:)
    // x^64 +
    // (t^11 + t^9 + t^8 + t^4 + t^3 + t)*x^3 + 
    // (t^10 + t^9 + t^8 + t^7 + t^6 + t^4 + t^3 + 1)*x^2 + 
    // (t^10 + t^9 + t^8 + t^7 + t^3 + t^2 + 1)*x + 
    // t^5 + t^4 + t^3 + t^2
    int gx[t+1][m+1] = {0,};
    gx[0][5] = gx[0][4] = gx[0][3] = gx[0][2] = 1;
    gx[1][10] = gx[1][9] = gx[1][8] = gx[1][7] = gx[1][3] = gx[1][2] = gx[1][0] = 1;
    gx[2][10] = gx[2][9] = gx[2][8] = gx[2][7] = gx[2][6] = gx[2][4] = gx[2][3] = gx[2][0] = 1;
    gx[3][11] = gx[3][9] = gx[3][8] = gx[3][4] = gx[3][3] = gx[3][1] = 1;
    gx[64][0] = 1;
    ctx_set_m_12(ctx, ft, gx, m, t);
#elif t == 96
    int gx[t+1] = {0,};
    gx[0]=1; gx[1]=1; gx[2]=1; gx[3]=1; gx[5]=1; gx[6]=1; gx[96]=1;
    ctx_set(ctx, ft, gx, m, t);
#elif t == 119
    int gx[t+1] = {0,};
    gx[0]=1; gx[8]=1; gx[119]=1;
    ctx_set(ctx, ft, gx, m, t);
#elif t == 128
    int gx[t+1] = {0,};
    gx[0]=1; gx[1]=1; gx[2]=1; gx[7]=1; gx[128]=1;
    ctx_set(ctx, ft, gx, m, t);
#endif

    printf("PRINT mod_gx\n");
    POLY_print(&ctx->mod_gx);
    printf("PRINT mod_coef\n");
    COEF_POLY_print(&ctx->mod_coef);
    printf("\n");
}

void test_gen_Xtable(POLY* Xtable, CTX* ctx, COEF_POLY fttable[])
{
    printf("\n========= GEN Xitable ===========\n");
    gen_Xitable(Xtable, ctx, fttable);
    for(int i = 0; i <= t; i++)
    {
        printf("[X^%d]\t\t", i+t);
        POLY_print(&Xtable[i]);
    }
}

void test_gen_fttable(COEF_POLY* fttable, CTX* ctx)
{
    printf("\n========= GEN fttable ===========\n");
    for(int i=0;i<=m;i++)
    {
        COEF_POLY_init(&fttable[i],0);
    }   
    coef_modft_table(fttable, ctx);
    for(int i = 0; i <= m; i++)
    {
        //printf("[t^%d]\t\t\n", i+m);    COEF_POLY_print(&fttable[i]);   printf("\n");
    }
}

void test_gen_Ttable(int Ttable[],int InvTtable[], COEF_POLY* fttable, CTX* ctx)
{
    printf("\n========= GEN Ttable ===========\n");
    
    gen_Ttable(Ttable, InvTtable, fttable, ctx ); // f2m税 乞窮 据社.
    printf("i   ttable   invttable\n");
    for(int i = 0; i < pow(2,m); i++)
    {
        //printf("[%d]\t\t%d\t\t%d\n", i,Ttable[i],InvTtable[i]);  
    }

}

void test_poly_add()
{
    printf("\n========= POLY_ADD ==========\n");
    POLY AA, BB, CC;
    int aa[3+1][MAX_COEF_POLY_DEGREE+1] = {{1,1,1,1,1,1,}, {1,0,0,0,0,1,}, {1,0,1,0,1,1,}, {0,1,0,1,0,1,}};
    int bb[5+1][MAX_COEF_POLY_DEGREE+1] = {{1,0,1,0,1,0,1,}, {1,1,0,0,1,1,}, {1,0,0,0,1,0,}, {1,0,0,1,0,1,}, {1,0,1,1,1,1,}, {1,0,1,1,0,1,}};
    POLY_init(&AA, 0);      POLY_init(&BB, 0);      POLY_init(&CC, 0);
    POLY_set(&AA, aa, 3, 5);    POLY_set(&BB, bb, 5, 6);

    POLY_add(&CC, &AA, &BB);   
    printf("A = ");         POLY_print(&AA);      
    printf("B = ");         POLY_print(&BB);      
    printf("A + B = ");     POLY_print(&CC);    
}

void test_coef_poly_mul(COEF_POLY fttable[], CTX* ctx)
{
    COEF_POLY A;
    int a[10] = {0,0,0,0,0,0,0,1,1,1};
    COEF_POLY B;
    int b[7] = {0,0,1,0,0,1,0};    
    COEF_POLY C;

    COEF_POLY_init(&A, 0); 
    COEF_POLY_set(&A, a, 9);
    COEF_POLY_init(&B, 0); 
    COEF_POLY_set(&B, b, 5);
    COEF_POLY_init(&C, 0); 

    printf("\n========= COEF_POLY MUL ===========\n");
    printf("A = ");     COEF_POLY_print(&A); printf("\n");
    printf("B = ");     COEF_POLY_print(&B); printf("\n");
    COEF_POLY_mul(&C, &A, &B, fttable, ctx) ;
    printf("A*B = ");     COEF_POLY_print(&C); printf("\n");
}

void test_poly_mul(POLY Xtable[], COEF_POLY fttable[], CTX* ctx)
{
    POLY A, B, C;
    int aa[2+1][MAX_COEF_POLY_DEGREE+1] = {{0,0,0,0,0,0,}, {0,0,0,0,0,1,}, {0,0,0,0,0,0,1,1,}};
    int bb[t+1][MAX_COEF_POLY_DEGREE+1] = {{0,0,0,0,0,0,0,}, {0,0,0,0,0,0,}, {0,0,0,0,0,0,}, {1,0,0,0,0,0,0,1,}};
    bb[t][4] = 1, bb[t][7] = 1;
    POLY_init(&A, 0);   POLY_init(&B, 0);   POLY_init(&C, 0);  
    POLY_set(&A, aa, 2, 7);    POLY_set(&B, bb, t, 7);

    printf("\n========= POLY MUL ===========\n");
    printf("A = ");     POLY_print(&A);
    printf("B = ");     POLY_print(&B);
    POLY_mul(&C, &A, &B, ctx, fttable, Xtable);
    printf("A*B = ");   POLY_print(&C);
}

void test_scalar_mul(COEF_POLY fttable[], CTX* ctx)
{
    COEF_POLY x;
    int a[10] = {0,0,0,0,0,0,0,1,1,1};
    COEF_POLY_init(&x, 0); 
    COEF_POLY_set(&x, a, 9);
    POLY A, B;
    int aa[3+1][MAX_COEF_POLY_DEGREE+1] = {{1,1,1,1,1,1,}, {1,0,0,0,0,1,}, {1,0,1,0,1,1,}, {0,1,0,1,0,1,}};
    POLY_init(&A, 0);       POLY_init(&B, 0);
    POLY_set(&A, aa, 3, 5);

    printf("\n========= scalar MUL ===========\n");
    printf("x = ");     COEF_POLY_print(&x);    printf("\n");
    printf("A = ");     POLY_print(&A);
    MULscalar(&B, &A, &x, fttable, ctx);
    printf("x*A = ");     POLY_print(&B);

    MULscalar_zzx(&A, &x, fttable, ctx);
    printf("x*A = ");   POLY_print(&A);
}

void test_xsqrt(COEF_POLY fttable[], int Ttable[], POLY Xtable[], CTX* ctx)
{
    printf("\n========= POLY Xsqrt ===========\n");
    POLY AA, AA_sqrt,AA_square;
    int aa[3+1][MAX_COEF_POLY_DEGREE+1] = {{1,0,0,0,0,0,}, {0,0,0,0,0,0,}, {0,0,0,0,0,0,}, {1,0,0,0,0,0,}};
    POLY_init(&AA, 0);     POLY_init(&AA_sqrt, 0);  POLY_init(&AA_square,0);
    POLY_set(&AA, aa, 3, 5);
    //printf("AA = ");        POLY_print(&AA);

    X_sqrt(&AA_sqrt, Xtable, &AA, fttable, Ttable,ctx);
    //printf("AA_sqrt = ");    POLY_print(&AA_sqrt);

    POLY_mul(&AA_square,&AA_sqrt,&AA_sqrt,ctx,fttable,Xtable);
    //printf("AA = ");        POLY_print(&AA_square);
    if(POLY_equal(&AA, &AA_square)== TRUE)
        printf("TRUE\n");
    else
        printf("FALSE\n");
}

int main()
{
    test_coef_poly();
    test_poly();

    CTX ctx;
    set_CTX(&ctx);
    COEF_POLY fttable[m+1];
    test_gen_fttable(fttable, &ctx);

    int Ttable[8192]={0,};
    int InvTtable[8192]={0,};
    test_gen_Ttable(Ttable, InvTtable, fttable, &ctx);
    
    POLY Xtable[t+1];
    test_gen_Xtable(Xtable, &ctx, fttable);

    test_poly_add();
    test_coef_poly_mul(fttable, &ctx);
    test_poly_mul(Xtable, fttable, &ctx);
    test_scalar_mul(fttable, &ctx);
    test_xsqrt(fttable, Ttable, Xtable, &ctx);

    //==============奄沙 魁 =====================
    //==Ri
    clock_t start, end; 
    double result_ac, result_af;

    POLY Ri[t/2];   //   Ri砺戚鷺 持失 
    POLY S,S_sqrt;
    for(int i=0;i<t/2;i++)
    {
        POLY_init(&S,0); POLY_init(&S_sqrt,0);
        S.coef[2*i+1].coef[0]=1;
        S.max_degree=2*i+1;
        S.coef[2*i+1].coef_max_degree=0;
        X_sqrt(&S_sqrt,Xtable,&S,fttable, Ttable,&ctx);
        POLY_copy(&Ri[i],&S_sqrt);
    }

    //==== algorithm 1==
    printf("\n========= Algorithm1 ===========\n");
    POLY Qx ,Qx_al1,Qx_al1_ver;
    int qq[t-1+1][MAX_COEF_POLY_DEGREE+1] = {{0,0,0,0,0,0,0,}, {0,0,0,0,0,0,}, {0,0,0,0,0,0,}, {1,0,0,0,0,0,0,1,}};
    qq[t-3][4] = 1;
    qq[t-10][4] = 1;
    qq[t-1][4] = 1, qq[t-1][7] = 1;
    for(int i = 0; i < t; i++)
    {
        for(int j = 0; j < 2 ; j++)
        {
            qq[i][j]=1;
        }
    }
    POLY_init(&Qx, 0); 
    POLY_init(&Qx_al1,0); 
    POLY_init(&Qx_al1_ver,0);
    POLY_set(&Qx, qq, t-1, 7);

    start = clock();
    X_sqrt(&Qx_al1, Xtable, &Qx, fttable, Ttable, &ctx);
    end = clock();
    result_ac = (double)(end - start)/(double)CLOCKS_PER_SEC;
    //printf("%f\n",result_ac);
    //printf("Qx="); POLY_print(&Qx); 
    //printf("Qx_result="); POLY_print(&Qx_al1); 
    POLY_mul(&Qx_al1_ver,&Qx_al1,&Qx_al1,&ctx,fttable,Xtable);
    //printf("Qx_ver="); POLY_print(&Qx_al1_ver);
    if(POLY_equal(&Qx, &Qx_al1_ver)== TRUE)
        printf("TRUE\n");
    else
        printf("FALSE\n");

    printf("\n========= Algorithm2 ===========\n");
    POLY Qx_al2,Qx_al2_ver;
    POLY Qx_odd;
    int z,a;
    POLY_init(&Qx_al2,0); 
    POLY_init(&Qx_al2_ver,0);
    COEF_POLY tmp;

    start = clock();
    for(int i=0;i<=(t-1)/2;i++)
    {
        vec2int(&a,Qx.coef[2*i].coef,Qx.coef[2*i].coef_max_degree+1);
        z = InvTtable[a];
        int2vec(Qx_al2.coef[i].coef,&Qx_al2.coef[i].coef_max_degree,z);
    }
    for(int i=0;i<=(t/2)-1;i++)
    {
        COEF_POLY_init(&tmp,0);
        POLY_init(&Qx_odd,0);
        vec2int(&a,Qx.coef[2*i+1].coef,Qx.coef[2*i+1].coef_max_degree+1);
        z = InvTtable[a];
        int2vec(tmp.coef,&tmp.coef_max_degree,z);
        MULscalar(&Qx_odd,&Ri[i],&tmp,fttable,&ctx);
        POLY_add_zzx(&Qx_al2,&Qx_odd);
    }
    end = clock();
    result_af = (double)(end - start)/(double)CLOCKS_PER_SEC;   
    //printf("%f\n",result_af);

    //printf("Qx_result="); POLY_print(&Qx_al2); 
    POLY_mul(&Qx_al2_ver,&Qx_al2,&Qx_al2,&ctx,fttable,Xtable);
   // printf("Qx_ver="); POLY_print(&Qx_al2_ver); 
    if(POLY_equal(&Qx, &Qx_al2_ver)== TRUE)
        printf("TRUE\n");
    else
        printf("FALSE\n");
    
    printf("[al1 | al2] %f %f\n", result_ac, result_af);
    return 0;
}