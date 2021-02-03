#include "poly.h"

void test_poly()
{
    printf("\n========= POLY set ===========\n");
    POLY AA, AA_copy;
    int aa[3+1] = {0b111111, 0b100001, 0b110101, 0b101010};
    POLY_init(&AA, 0);
    POLY_set(&AA, aa, 3);
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
    int ft = 0b1000011101011; 
#elif m == 13
    int ft = 0b10000000011011;
#endif

    ctx_init(ctx);
#if t == 64 // 64託ぞぞぞぞぞぞぞぞ 鯵旭精暗:)
    // x^64 +
    // (t^11 + t^9 + t^8 + t^4 + t^3 + t)*x^3 + 
    // (t^10 + t^9 + t^8 + t^7 + t^6 + t^4 + t^3 + 1)*x^2 + 
    // (t^10 + t^9 + t^8 + t^7 + t^3 + t^2 + 1)*x + 
    // t^5 + t^4 + t^3 + t^2
    int gx[t+1] = {0,};
    gx[0] = 0b111100;                                   //gx[0][5] = gx[0][4] = gx[0][3] = gx[0][2] = 1;
    gx[1] = 0b11110001101;                            //gx[1][10] = gx[1][9] = gx[1][8] = gx[1][7] = gx[1][3] = gx[1][2] = gx[1][0] = 1;
    gx[2] = 0b11111011001;                          //gx[2][10] = gx[2][9] = gx[2][8] = gx[2][7] = gx[2][6] = gx[2][4] = gx[2][3] = gx[2][0] = 1;
    gx[3] = 0b101100011010;                         //gx[3][11] = gx[3][9] = gx[3][8] = gx[3][4] = gx[3][3] = gx[3][1] = 1;
    gx[64] = 1;
    ctx_set(ctx, ft, gx, t);
#elif t == 96
    int gx[t+1] = {0,};
    gx[0]=1; gx[1]=1; gx[2]=1; gx[3]=1; gx[5]=1; gx[6]=1; gx[96]=1;
    ctx_set(ctx, ft, gx, t);
#elif t == 119
    int gx[t+1] = {0,};
    gx[0]=1; gx[8]=1; gx[119]=1;
    ctx_set(ctx, ft, gx, t);
#elif t == 128
    int gx[t+1] = {0,};
    gx[0]=1; gx[1]=1; gx[2]=1; gx[7]=1; gx[128]=1;
    ctx_set(ctx, ft, gx, t);
#endif

    printf("PRINT mod_gx\n");
    POLY_print(&ctx->mod_gx);
    printf("PRINT mod_coef\n");
    COEF_POLY_print(ctx->mod_coef);
    printf("\n");
}

void test_gen_Xtable(POLY* Xtable, CTX* ctx, int fttable[])
{
    printf("\n========= GEN Xitable ===========\n");
    gen_Xitable(Xtable, ctx, fttable);
    for(int i = 0; i <= t; i++)
    {
        printf("[X^%d]\t\t", i+t);
        POLY_print(&Xtable[i]);
    }
}

void test_gen_fttable(int fttable[], CTX* ctx)
{
    printf("\n========= GEN fttable ===========\n");
    for(int i=0;i<=m;i++)
        fttable[i] = 0;
    coef_modft_table(fttable, ctx);
    for(int i = 0; i <= m; i++)
    {
        printf("[t^%d]\t\t", i+m);    COEF_POLY_print(fttable[i]);   printf("\n");
    }
}

void test_gen_Ttable(int Ttable[],int InvTtable[], int fttable[], CTX* ctx)
{
    printf("\n========= GEN Ttable ===========\n");
    
    gen_Ttable(Ttable, InvTtable, fttable, ctx); // f2m税 乞窮 据社.
    printf("i   ttable   invttable\n");
    for(int i = 0; i < pow(2,m); i++)
    {
        printf("[%d]\t\t%d\t\t%d\n", i,Ttable[i],InvTtable[i]);  
    }

}

void test_poly_add()
{
    printf("\n========= POLY_ADD ==========\n");
    POLY AA, BB, CC;
    int aa[3+1] = {0b111111, 0b100001, 0b110101, 0b101010};
    int bb[5+1] = {0b1010101, 0b110011, 0b10001, 0b101001, 0b111101, 0b101101};
    POLY_init(&AA, 0);      POLY_init(&BB, 0);      POLY_init(&CC, 0);
    POLY_set(&AA, aa, 3);    POLY_set(&BB, bb, 5);

    POLY_add(&CC, &AA, &BB);   
    printf("A = ");         POLY_print(&AA);      
    printf("B = ");         POLY_print(&BB);      
    printf("A + B = ");     POLY_print(&CC);    
}

void test_coef_poly_mul(int fttable[], CTX* ctx)
{
    int A, B, C;
    A = 0b1110000000;
    B = 0b100100;

    printf("\n========= COEF_POLY MUL ===========\n");
    printf("A = ");     COEF_POLY_print(A); printf("\n");
    printf("B = ");     COEF_POLY_print(B); printf("\n");
    COEF_POLY_mul(&C, A, B, fttable, ctx) ;
    printf("A*B = ");     COEF_POLY_print(C); printf("\n");
}

void test_poly_mul(POLY Xtable[], int fttable[], CTX* ctx)
{
    POLY A, B, C;
    int aa[2+1] = {0b0, 0b100000, 0b11000000};
    int bb[t-1+1] = {0b0, 0b0, 0b0, 0b10000001,};
    bb[t-1] = 0b1000001;
    POLY_init(&A, 0);   POLY_init(&B, 0);   POLY_init(&C, 0);  
    POLY_set(&A, aa, 2);    POLY_set(&B, bb, t);

    printf("\n========= POLY MUL ===========\n");
    printf("A = ");     POLY_print(&A);
    printf("B = ");     POLY_print(&B);
    POLY_mul(&C, &A, &B, ctx, fttable, Xtable);
    printf("A*B = ");   POLY_print(&C);
}

void test_scalar_mul(int fttable[], CTX* ctx)
{
    int x = 0b1110000000;
    POLY A, B;
    int aa[3+1] = {0b111111, 0b100001, 0b110101, 0b101010};
    POLY_init(&A, 0);       POLY_init(&B, 0);
    POLY_set(&A, aa, 3);

    printf("\n========= scalar MUL ===========\n");
    printf("x = ");     COEF_POLY_print(x);    printf("\n");
    printf("A = ");     POLY_print(&A);
    MULscalar(&B, &A, x, fttable, ctx);
    printf("x*A = ");     POLY_print(&B);

    MULscalar_zzx(&A, x, fttable, ctx);
    printf("x*A = ");   POLY_print(&A);
}

void test_xsqrt(int fttable[], int Ttable[], POLY Xtable[], CTX* ctx)
{
    printf("\n========= POLY Xsqrt ===========\n");
    POLY AA, AA_sqrt,AA_square;
    int aa[5+1] = {0b1010101, 0b110011, 0b10001, 0b101001, 0b111101, 0b101101};
    POLY_init(&AA, 0);     POLY_init(&AA_sqrt, 0);  POLY_init(&AA_square,0);
    POLY_set(&AA, aa, 5);
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
    test_poly();

    CTX ctx;
    set_CTX(&ctx);

    int fttable[m+1];
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
    POLY S;
    //printf("R[i]============\n");
    for(int i=0;i<t/2;i++)
    {
        POLY_init(&S,0);
        //POLY_init(&S_sqrt,0);
        S.coef[2*i+1] = 1;
        S.max_degree = 2*i+1;
        X_sqrt(&Ri[i],Xtable,&S,fttable, Ttable,&ctx);
        //POLY_copy(&Ri[i],&S_sqrt);
        //POLY_print(&Ri[i]);
    }

    //==== algorithm 1==
    printf("\n========= Algorithm1 ===========\n");
    POLY Qx ,Qx_al1,Qx_al1_ver;
    int qq[t-1+1] = {0b0, 0b0, 0b0, 0b10000001,};
    // for(int i = 0; i < t; i++)
    // {
    //     for(int j = 0; j < 4; j++)
    //     {
    //         qq[i] ^= 1 << j;
    //     }
    // }
    POLY_init(&Qx, 0); 
    POLY_init(&Qx_al1,0); 
    POLY_init(&Qx_al1_ver,0);
    POLY_set(&Qx, qq, t-1);

    start = clock();
    X_sqrt(&Qx_al1, Xtable, &Qx, fttable, Ttable, &ctx);
    end = clock();
    result_ac = (double)(end - start)/(double)CLOCKS_PER_SEC;
    //printf("%f\n",result_ac);
    printf("Qx="); POLY_print(&Qx); 
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
    POLY_init(&Qx_al2,0); 
    POLY_init(&Qx_al2_ver,0);

    start = clock();
    for(int i=0;i<=(t-1)/2;i++)
    {
        Qx_al2.coef[i] = InvTtable[Qx.coef[2*i]];
    }

    for(int i=0;i<=(t/2)-1;i++)
    {
        POLY_init(&Qx_odd,0);
        //tmp = Qx.coef[2*i+1];
        //printf("\nri="); POLY_print(&Ri[i]);
        //COEF_POLY_print(Qx.coef[2*i+1]);
        MULscalar(&Qx_odd,&Ri[i],InvTtable[Qx.coef[2*i+1]],fttable,&ctx);
        //printf("Qxodd="); POLY_print(&Qx_odd);
        POLY_add_zzx(&Qx_al2,&Qx_odd);
    }
    //printf("Qx_Al2="); POLY_print(&Qx_al2);

    end = clock();
    result_af = (double)(end - start)/(double)CLOCKS_PER_SEC;   
    //printf("%f\n",result_af);

    //printf("Qx_result="); POLY_print(&Qx_al2); 
    POLY_mul(&Qx_al2_ver,&Qx_al2,&Qx_al2,&ctx,fttable,Xtable);
    printf("Qx_ver="); POLY_print(&Qx_al2_ver); 
    if(POLY_equal(&Qx_al1, &Qx_al2)== TRUE)
        printf("TRUE\n");
    else
        printf("FALSE\n");

    if(POLY_equal(&Qx, &Qx_al2_ver)== TRUE)
        printf("TRUE\n");
    else
        printf("FALSE\n");
    
    printf("[al1 | al2] %f %f\n", result_ac, result_af);

    return 0;
}