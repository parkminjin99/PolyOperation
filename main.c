#include "poly.h"

void test_coef_poly()
{
    int zero;
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

int main()
{
    test_coef_poly();
    test_poly();

    printf("\n========= Set CTX ===========\n");
    int fttable[12];
    CTX ctx;
    int ft[m+1]={1,1,0,1,1,0,0,0,0,0,0,0,0,1};
    int gx[t+1] = {0,}; //x^96+x^6+x^5+x^3+x^2+x+1   
    gx[0]=1; gx[1]=1; gx[2]=1; gx[3]=1; gx[5]=1; gx[6]=1; gx[96]=1;
    ctx_init(&ctx);
    ctx_set(&ctx, ft, gx, m, t);
    printf("print mod_gx\n");
    POLY_print(&ctx.mod_gx);
    printf("print mod_coef\n");
    COEF_POLY_print(&ctx.mod_coef);
    printf("\n");

    printf("\n========= GEN Xitable ===========\n");
    POLY Xtable[t+1];
    gen_Xitable(Xtable, &ctx);
    

    return 0;
}
