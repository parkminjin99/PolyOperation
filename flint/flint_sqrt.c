 /**********************
 * ����ؾ��ϴ� �Ķ����
 * m  12  13  13    13
 * t  64  96  128   119
 ***********************/

#include <time.h>
#include <stdlib.h>
#include <stdio.h>

#include "fq_nmod.h"
#include "ulong_extras.h"
#include "long_extras.h"
#include "fq_nmod_poly.h"
#include "nmod_poly.h"
#include "fmpz.h"

#define UPBOUND 10000  // 2^m -1 ���� ���̺� ����
#define FALSE 0    // �� ���� ���ϴ� �������� false�� ��� 0�� �����.

int main()
{
    double cloal=0, cloex=0;

    int cou=1;

    for(int count=0;count<cou;count++){

    fmpz_t q, twom, twom_1, twomt_1, two;
    slong  m, z, t;
    ulong  u1, u2, u3;
    clock_t c0, c1;
    fq_nmod_ctx_t ctx;
    fq_nmod_t qp[UPBOUND], qi[UPBOUND], coef, table[UPBOUND],table2[UPBOUND];
    fq_nmod_poly_t Qx, g, Sq, x, r[UPBOUND], result, qod, qev, q2t;
    nmod_poly_t np1, np2;
    //flint_rand_t state;
    //flint_randinit(state);
    FLINT_TEST_INIT(state);    // ���� ������ ���� �õ� ����.

    //============���׽� ����ü ����========================================
    fq_nmod_poly_init(Qx, ctx);  // Q(x)
    fq_nmod_poly_init(g, ctx);   // g(x)
    fq_nmod_poly_init(Sq, ctx);  // Sq(x)
    fq_nmod_poly_init(x, ctx);   // x(x) = x
    //============����ڰ� �����ϴ� ����===============================
    fmpz_init_set_ui(q, 2);       // Fqm���� q
    m = 13;                       // Fqm���� m
    t = 96;                       // Fqm(X)�� ���� ����. g(x)�� t, Q(x)�� t-1�� ���׽�.
    //============================================================
    fq_nmod_ctx_init_conway(ctx, q, m, "t");                   // ctx�� Fqm����.
    //printf("\n======ctx info======\n");    fq_nmod_ctx_print(ctx);   printf("====================\n"); // ����� ctx���� Ȯ��
    fq_nmod_poly_gen(x, ctx);                                  // 1�����׽� x(x)=x����.
    // ���׽� Qx,g�� �����ϰ� ����.
    //fq_nmod_poly_randtest(Qx, state, t, ctx);                  // t���� ���� ���� Qx. �� t-1���� ��.
    //fq_nmod_poly_print_pretty(Qx,"X",ctx);

    //fq_nmod_poly_randtest_irreducible(g, state, t+1, ctx);     // t+1���� ���� ���� g. �� t���� ��. �� g�� �����׽����� ����.
    

    //=============================Gx ����=======================
    fmpz_t fmz,fm;
    fmpz_init(fmz);
    fmpz_init(fm);

    fmpz_one(fmz);
    fq_nmod_poly_set_coeff_fmpz(g, 96, fmz, ctx);
    fq_nmod_poly_set_coeff_fmpz(g, 0, fmz, ctx);
    fq_nmod_poly_set_coeff_fmpz(g, 1, fmz, ctx);
    fq_nmod_poly_set_coeff_fmpz(g, 2, fmz, ctx);
    fq_nmod_poly_set_coeff_fmpz(g, 3, fmz, ctx);
    fq_nmod_poly_set_coeff_fmpz(g, 5, fmz, ctx);
    fq_nmod_poly_set_coeff_fmpz(g, 6, fmz, ctx); //Gx ���� �Ϸ�

    fmpz_init(fmz);
    ulong c=2;
    fmpz_set_ui(fmz, c);
    //fmpz_init_set_ui(fmz,c);
    fq_nmod_poly_set_coeff_fmpz(Qx, 95, fmz, ctx);
    c=44;
    fmpz_set_ui(fm, c);
    fmpz_print(fm); printf("\n");
    fq_nmod_poly_set_coeff_fmpz(Qx, 90, fm, ctx);   //Qx������ �̻�

    printf("gx = "); fq_nmod_poly_print_pretty(g,"X",ctx);     printf("\n");
    printf("Qx = "); fq_nmod_poly_print_pretty(Qx,"X",ctx);    printf("\n");

    
    //printf("%d",fq_nmod_poly_is_irreducible(g,ctx));
    //==================== pre-computation ===========================================================
    //=========================== Ri(X) ������� LINE 2 ==================================
    fmpz_init(two);               // 2 ����
    fmpz_init(twomt_1);           // 2^(mt-1) ����
    fmpz_set_ui(two, 2);
    fmpz_pow_ui(twomt_1, two, (m*t-1));
    fq_nmod_poly_powmod_fmpz_binexp(Sq, x, twomt_1, g, ctx);    // Sq(x) = x^(2^(mt-1)) mod g(x)
    //printf("Sq(x) = ");     fq_nmod_poly_print_pretty(Sq, "X", ctx);    printf("\n");   // Sq Ȯ��
    for(ulong i = 0; i <= (t/2 - 1); i++)
    {
        fq_nmod_poly_init(r[i], ctx);                           // r[i] ����ü ����.
        fq_nmod_poly_powmod_ui_binexp(r[i], x, i, g, ctx);      // r[i] = x^i mod g(x)
        fq_nmod_poly_mulmod(r[i], r[i], Sq, g, ctx);            // r[i] = r[i]*Sq mod g(x)
        //printf("ri(x) = ");     fq_nmod_poly_print_pretty(r[i],"X",ctx);    printf("\n");   //ri Ȯ��
    }
    //======== ����������̺�(Fqm�� ��� ���ҿ� ���� 2^(m-1)������ �̸� �س��� ���̺�) LINE 3 =========
    fmpz_init(twom_1);                                      // twom_1 = 2^(m-1). (fmpz)
    fmpz_init(twom);                                        // twom = 2^m        (fmpz)
    fmpz_set_ui(twom_1, pow(2,m-1));
    fmpz_set_ui(twom, pow(2,m));
    u3 = fmpz_get_ui(twom);                                 // u3 = 2^m (ulong)
    //printf("2^(m-1) = ");  fmpz_print(twom_1);  printf("\n");

    fq_nmod_init(coef,ctx);
    for(ulong i = 0; i < u3; i++)
    {
        nmod_poly_init(np1, m);      // ������ m�� �̸�.
        fq_nmod_init(table[i], ctx); // fqm���� ��� ���ҿ� ���� 2^twom_1 �� table�� ����.
        u1 = i;
        for(ulong j = 0; j < m; j++) // �ش� For���� �־��� j��� ������ ����ü 2^m���� ���ҷ� �ٲٴ� ������ ��ģ��.
        {
            u2 = u1 % 2;
            u1 = u1 / 2;
            nmod_poly_set_coeff_ui(np1, j, u2); // 2�� �����鼭 fqm�� ���ҷ� ǥ��.
        }
        // ========���̺� ���� Ȯ��=========
        // printf("%lu  ", i);                                // i��° ���̺�
        // nmod_poly_print_pretty(np1, "x");  printf("  ");   // nmod ���� ���
        fq_nmod_set_nmod_poly(coef, np1, ctx);                // nmod�� fqnmod ���·�
        fq_nmod_pow(table[i], coef, twom_1 ,ctx);             // table = coef^twom_1
        //fq_nmod_print_pretty(coef, ctx); printf(" ");      // fqnmod ���� ���
        //fq_nmod_print_pretty(table[i], ctx); printf("\n"); // table ���
        nmod_poly_clear(np1);
    }
    //======================= q^2 �������
    //printf("2^(m-1) = ");  fmpz_print(twom_1);  printf("\n");

    for(ulong i = 0; i < u3; i++)
    {
        nmod_poly_init(np1, m);      // ������ m�� �̸�.
        fq_nmod_init(table2[i], ctx); // fqm���� ��� ���ҿ� ���� 2^twom_1 �� table�� ����.
        u1 = i;
        for(ulong j = 0; j < m; j++) // �ش� For���� �־��� j��� ������ ����ü 2^m���� ���ҷ� �ٲٴ� ������ ��ģ��.
        {
            u2 = u1 % 2;
            u1 = u1 / 2;
            nmod_poly_set_coeff_ui(np1, j, u2); // 2�� �����鼭 fqm�� ���ҷ� ǥ��.
        }
        // ========���̺� ���� Ȯ��=========
        // printf("%lu  ", i);                                // i��° ���̺�
        // nmod_poly_print_pretty(np1, "x");  printf("  ");   // nmod ���� ���
        fq_nmod_set_nmod_poly(coef, np1, ctx);                // nmod�� fqnmod ���·�
        fq_nmod_pow(table2[i], coef, two ,ctx);             // table = coef^twom_1
        // fq_nmod_print_pretty(coef, ctx); printf(" ");      // fqnmod ���� ���
        //fq_nmod_print_pretty(table2[i], ctx); printf("\n"); // table ���
        nmod_poly_clear(np1);
    }
    //===================================================================
    //x^t�Ѵ� ���� mod g(x)����ϱ�.
    fmpz_t extwo[UPBOUND];
    fq_nmod_poly_t xpow[UPBOUND];

    for(int i=0;i<t;i++)
    {
        fmpz_init(extwo[i]);
        fq_nmod_poly_init(xpow[i],ctx);
        fmpz_mul_ui(extwo[i],two,i);  // x^t ~ x^2t-2 ������ ������ �� �� ���
        fq_nmod_poly_powmod_fmpz_binexp(xpow[i],x,extwo[i],g,ctx); // �����̿��� x^t mod gx ���
        //fq_nmod_poly_print_pretty(xpow[i],"x",ctx);
        //printf("\n");
    }

    //================================================================================================

    //printf("x    = ");     fq_nmod_poly_print_pretty(x, "X", ctx);     printf("\n");  // xȮ��
    //printf("g(x) = ");     fq_nmod_poly_print_pretty(g, "X", ctx);     printf("\n");  // gȮ��
    //printf("Q(x) = ");     fq_nmod_poly_print_pretty(Qx, "X", ctx);    printf("\n");  // QȮ��

    nmod_poly_init(np2, t);
    for(slong i=0; i < t; i++)
    {
        fq_nmod_init(qi[i], ctx);                          // fqm���� qi[i]
        fq_nmod_init(qp[i], ctx);
    }
    fq_nmod_poly_init(result, ctx);
    fq_nmod_poly_init(qev, ctx); fq_nmod_poly_init(qod, ctx);

    int pow2[20]={0,};
    for(int j=0;j<m;j++)
    {
        pow2[j]=pow(2, j);
    }

    //============= LINE1 ==========================
    c0 = clock();                                            // time start
    // LINE 1�� LINE 2�� �����̳� LINE 2�� ����������� ���������Ƿ� ������ �ʿ����.
    // printf("Sq(x) = ");     fq_nmod_poly_print_pretty(Sq, "X", ctx);    printf("\n");  //Sq Ȯ��

    //============= LINE2 ==========================
    // LINE 2�� ��������� ���� ����.

    //============= LINE3 ==========================
    for(slong i=0; i < t; i++)
    {
        //nmod_poly_init(np2, t);                            // t-1�� ���� -> �ʱ�ȭ
        nmod_poly_zero(np2);                               // t-1�� ���� -> �ʱ�ȭ
        //fq_nmod_init(qi[i], ctx);                          // fqm���� qi[i]
        fq_nmod_poly_get_coeff(qi[i], Qx, i, ctx);         // �Է¹��� Qx�� ����� ����. ����� fqm���� ����
        fq_nmod_get_nmod_poly(np2, qi[i], ctx);            // �ش� ���Ҹ� x�̿� ������ ���·� ��ȯ
        //printf("coef = "); nmod_poly_print_pretty(np2, "x");

        z = 0;                                             // ���� Fqm���Ҹ� ���ڷ� �ٲپ� ����.
        for(slong j = 0; j < m; j++)                       // �ش� for���� Fqm���Ҹ� ���ڷ� �ٲپ� ����.
        {
            z += nmod_poly_get_coeff_ui(np2, j) * pow2[j];                   // �������·� ��ȯ
        }
        //==============��� Ȯ�� ===============================
        //printf(" z=%ld ",z);                                                   // ���� x������ ��������
        //printf("table= "); fq_nmod_print_pretty(table[z], ctx); printf("\n");  // �ش����̺� ����Ȱ�
        //fq_nmod_init(qp[i], ctx);
        fq_nmod_set(qp[i], table[z], ctx);                                     // ��������� ���� qp�� ����
        //printf("qp= "); fq_nmod_print_pretty(qp[i], ctx);

        //nmod_poly_print_pretty(np2, "x"); printf("\n");
        //nmod_poly_clear(np2);
    }

    //============= LINE4, 5, 6 ==========================
    //fq_nmod_poly_init(result, ctx);   // ���� ��� ���� ����
    fq_nmod_poly_zero(result, ctx);   // result�� 0���� �ʱ�ȭ

    for(int i = 0; i < t; i++)
    {
        //printf("qp= "); fq_nmod_print_pretty(qp[i], ctx); printf("  ");            // ��Ȯ��
        if(i % 2 == 0)  // index�� ¦���� ��� LINE 4�� ���� ���� 
        {
            //fq_nmod_poly_init(qev, ctx);
            fq_nmod_poly_zero(qev,ctx);                                          // qp�� index�� ¦���� ��� �߰��� ���忡 ���
            fq_nmod_poly_set_coeff(qev, i/2, qp[i], ctx);                            // qev�� qp[i]x^(i/2)�� ����.
            fq_nmod_poly_add(result, result, qev, ctx);                              // �ش� ���� ��������� ����
            //printf("q'= "); fq_nmod_print_pretty(qp[i], ctx);                        // ���Ȯ��
            //printf("qev= "); fq_nmod_poly_print_pretty(qev, "X", ctx); printf("  "); // ���Ȯ��
            //fq_nmod_poly_clear(qev, ctx);                                            // ���� �ʱ�ȭ
        }
        else // index�� Ȧ���� ��� LINE 5�� ���� ����
        {
            //fq_nmod_poly_init(qod, ctx);                                                  // qp�� index�� Ȧ���� ��� �߰��� ���忡 ���
            fq_nmod_poly_zero(qod,ctx);
            fq_nmod_poly_scalar_mul_fq_nmod(qod, r[i/2], qp[i], ctx);                     // qod�� r[i/2]�� qp[i]�� �� ���׽����� ����
            fq_nmod_poly_add(result, result, qod, ctx);                                   // �ش� ���� ��������� ����
            //printf("q'= "); fq_nmod_print_pretty(qp[i], ctx);                             // ���Ȯ��
            //printf("r[i]= "); fq_nmod_poly_print_pretty(r[i/2], "X", ctx); printf("  ");  // ���Ȯ��
            //printf("qr[i]="); fq_nmod_poly_print_pretty(qod, "X", ctx);    printf("  ");  // ���Ȯ��
            //fq_nmod_poly_clear(qod, ctx);
        }
        //printf("return(x) = "); fq_nmod_poly_print_pretty(result, "X", ctx);  printf("\n");   // ri
    }
    c1 = clock();                           // time end
    cloal += (double) (c1 - c0) / CLOCKS_PER_SEC;

    nmod_poly_clear(np2);
    fq_nmod_poly_clear(qev, ctx);     // ���� �ʱ�ȭ
    fq_nmod_poly_clear(qod, ctx);


    //printf("\n\n");
    //printf("g(x)      = ");    fq_nmod_poly_print_pretty(g, "X", ctx);      printf("\n\n");  // gȮ��
    //printf("Q(x)      = ");    fq_nmod_poly_print_pretty(Qx, "X", ctx);     printf("\n\n");  // QȮ��
    //printf("return(x) = ");    fq_nmod_poly_print_pretty(result, "X", ctx); printf("\n\n");  // sqrt(Qx) ���

    fq_nmod_poly_t calre;                               // ��갪�� ���� ���� ������ ���ϱ� ����.  ��������� result�� ������ Q(x)�� �������� Ȯ��

    fq_nmod_poly_init(calre, ctx);
    fq_nmod_poly_mulmod(calre, result, result, g, ctx); //result�� ������ mod g(x)�� ���� ����Ͽ� calre�� ����.

    //printf("cal       = ");    fq_nmod_poly_print_pretty(calre, "X", ctx);  printf("\n");  //ri

    int a,b;

    a = fq_nmod_poly_equal(Qx,calre,ctx);               // Qx == calre �� ��� 1���, �ƴϸ� 0���

    //printf("\n\n");
    //printf("m = %ld, t = %ld \n", m, t);
    //printf("%d",a);
    if(a==FALSE)
    {
        printf("\n FALSE \n");
        break;
    }
    else
        printf("\n TRUE Q(x) == return^2 11\n");

    //printf("time algori : %fs\n", cloal);     // �˰��� ���ð�

    //================== ������ ���� ================================
    // fqm[X]/g(X)�� ���� f(X)�� ���� f(X)^(2^(mt-1)) mod g(X)�� sqrt(f(X))�� ���� �̿�.

    //fq_nmod_poly_t expone;
    //fq_nmod_poly_init(expone,ctx);
    //fq_nmod_poly_powmod_fmpz_binexp(expone, Qx, twomt_1, g, ctx);
    fq_nmod_poly_init(q2t, ctx);
    //fq_nmod_poly_zero(q2t, ctx);                                          // qp�� index�� ¦���� ��� �߰��� ���忡 ���

    nmod_poly_init(np2, t);                            // t-1�� ���� -> �ʱ�ȭ
    c0 = clock();                       // time start
    for(int x=0;x<m*t-1;x++){
        fq_nmod_poly_zero(q2t,ctx);

        for(slong i=0; i < t; i++)
        {
            fq_nmod_poly_get_coeff(qi[i], Qx, i, ctx);         // �Է¹��� Qx�� ����� ����. ����� fqm���� ����
            fq_nmod_get_nmod_poly(np2, qi[i], ctx);            // �ش� ���Ҹ� x�̿� ������ ���·� ��ȯ
            z = 0;                                             // ���� Fqm���Ҹ� ���ڷ� �ٲپ� ����.
            for(slong j = 0; j < m; j++)                       // �ش� for���� Fqm���Ҹ� ���ڷ� �ٲپ� ����.
            {
                z += nmod_poly_get_coeff_ui(np2, j) * pow2[j];                   // �������·� ��ȯ
            }



            fq_nmod_set(qp[i], table2[z], ctx);                                     // ��������� ���� qp�� ����

            if(2*i<t)
                fq_nmod_poly_set_coeff(q2t, 2*i, qp[i], ctx);
            else
            {
                fq_nmod_poly_scalar_addmul_fq_nmod(q2t,xpow[i],qp[i],ctx);
                //fq_nmod_poly_set_coeff(q2t, 2*i, qp[i], ctx);
            }

        }
        //fq_nmod_poly_print_pretty(q2t,"x",ctx);
        fq_nmod_poly_set(Qx,q2t,ctx);
        //fq_nmod_poly_rem(Qx,q2t,g,ctx);
        //fq_nmod_poly_print_pretty(Qx,"x",ctx);
    }

    //printf("Sq(x) = ");     fq_nmod_poly_print_pretty(h,"X",ctx);    printf("\n");  //Sq Ȯ��

    c1 = clock();                       // time end
    nmod_poly_clear(np2);

    cloex += (double) (c1 - c0) / CLOCKS_PER_SEC;

    b = fq_nmod_poly_equal(result,Qx,ctx);          // ���� �˰���� ���� ����� �������� Ȯ��.

    if(b==FALSE) //��� Ȯ�ο�
    {
        printf("\n FALSE \n");
        break;
    }
    else
        printf(" TRUE sqrt(Q(x)) == return \n");

    //printf("time expone : %fs\n", cloex); // ������ ���ð�
    //printf("\n\n");

    //====================== �޸� ���� ============================
    fmpz_clear(q);  fmpz_clear(twom);  fmpz_clear(twom_1);  fmpz_clear(twomt_1);  fmpz_clear(two);

    for (int i = 0; i < t; i++)
    {
        fq_nmod_clear(qp[i],ctx);
        fq_nmod_clear(qi[i],ctx);
    }
    for (int i = 0; i < u3; i++)
    {
        fq_nmod_clear(table[i],ctx);         fq_nmod_clear(table2[i],ctx);
    }
    for (int i = 0; i < t/2 ; i++)
        fq_nmod_poly_clear(r[i],ctx);

    for(int i=0;i<t;i++)
    {
        fmpz_clear(extwo[i]);
        fq_nmod_poly_clear(xpow[i],ctx);
    }

    fq_nmod_clear(coef,ctx);

    fq_nmod_poly_clear(Qx,ctx);
    fq_nmod_poly_clear(g,ctx);
    fq_nmod_poly_clear(Sq,ctx);
    fq_nmod_poly_clear(x,ctx);
    fq_nmod_poly_clear(result,ctx);
    fq_nmod_poly_clear(calre,ctx);
    //fq_nmod_poly_clear(expone,ctx);
    fq_nmod_poly_clear(q2t,ctx);
    fq_nmod_ctx_clear(ctx);
    printf("%d\n",count);

    }

    printf("\n");
    printf("time algori : %fs\n", cloal/cou);     // �˰��� ���ð�
    printf("time expone : %fs\n", cloex/cou);     // ������ ���ð�


    return EXIT_SUCCESS;
}