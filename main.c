#include "poly.h"

int main()
{
    int zero;
    printf("Start\n");
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
    
    printf("\n========= POLY set ===========\n");
    POLY AA, AA_copy; // 3차 poly (계수는 4차)
    int aa[3+1][MAX_COEF_POLY_DEGREE+1] = {{1,1,1,1,1,1,}, {1,0,0,0,0,1,}, {1,0,1,0,1,1,}, {0,1,0,1,0,1,}};
    POLY_init(&AA, 0); // init 함수 사용 시에는 항상 최고차수를 0으로 설정하기 
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
    printf("\nEnd\n");
    return 0;
}
