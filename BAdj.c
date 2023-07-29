#include <gmp.h>
#include <stdio.h>
#include <time.h>

void AtomicBlockJacobianDoubling(mpz_t jx1, mpz_t jy1, mpz_t jz1, mpz_t p);

int main()
{
    //Parámetros de la curva P-192
    const char primo[] = "6277101735386680763835789423207666416083908700390324961279";
    const char EC_B[] = "64210519 e59c80e7 0fa7e9ab 72243049 feb8deec c146b9b1";
    const char EC_X1[] = "188da80e b03090f6 7cbf20eb 43a18800 f4ff0afd 82ff1012";
    const char EC_Y1[] = "07192b95 ffc8da78 631011ed 6b24cdd5 73f977a1 1e794811";
    const char SEED[] = "3045ae6f c8422f64 ed579528 d38120ea e12196d5";

    mpz_t p, a, b, dos, tres, seed, x1, y1, jx1, jy1, jz1;
    clock_t time1, time2;

    mpz_init(p);
    mpz_init(a);
    mpz_init(b);
    mpz_init(dos);
    mpz_init(tres);
    mpz_init(seed);
    mpz_init(x1);
    mpz_init(y1);
    mpz_init(jx1);
    mpz_init(jy1);
    mpz_init(jz1);
    

    //Set de parámetros de la curva: p, a, b;
    mpz_set_str(p, primo, 10);
    mpz_set_ui(dos, 2);
    mpz_set_ui(tres, 3);
    mpz_sub(a, p, tres);
    mpz_set_str(b, EC_B, 16);
    mpz_set_str(seed, SEED, 16);

    //Set del punto base G: (x1,y1)
    mpz_set_str(x1, EC_X1, 16);
    mpz_set_str(y1, EC_Y1, 16);

    mpz_set_str(jx1, EC_X1, 16);
    mpz_set_str(jy1, EC_Y1, 16);
    mpz_set_str(jz1, "1", 16);

    time1 = clock();
    AtomicBlockJacobianDoubling(jx1, jy1, jz1, p);
    time2 = clock();
    printf("Tiempo de ejecución: %f\n", (double)(time2 - time1) / CLOCKS_PER_SEC);



    return 0;
}

void AtomicBlockJacobianDoubling(mpz_t jx1, mpz_t jy1, mpz_t jz1, mpz_t p){
    //mpz R1 a R8
    mpz_t R1, R2, R3, R4, R5, R6, R7, R8, jx3, jy3, jz3;
    mpz_init(R1);
    mpz_init(R2);
    mpz_init(R3);
    mpz_init(R4);
    mpz_init(R5);
    mpz_init(R6);
    mpz_init(R7);
    mpz_init(R8);

    //R1 = jx1, R2 = jy1, R3 = jz1
    mpz_set(R1, jx1);
    mpz_set(R2, jy1);
    mpz_set(R3, jz1);

    //Block 1
    mpz_pow_ui(R4, R3, 2);  //R4 = R3^2
    mpz_neg(R5, R4);        //R5 = -R4
    mpz_add(R6, R1, R4);    //R6 = R1+R4
    mpz_add(R4, R1, R5);    //R4 = R1+R5
    mpz_mul(R5, R6, R4);    //R5 = R6*R4
    mpz_add(R4, R5, R5);    //R4 = R5+R5

    //Block 2
    mpz_pow_ui(R6, R2, 2);  //R6 = R2^2
    mpz_neg(R7, R1);        //R7 = -R1
    mpz_add(R1, R7, R7);    //R1 = R7+R7
    mpz_add(R7, R6, R6);    //R7 = R6+R6
    mpz_mul(R6, R1, R7);    //R6 = R1*R7
    mpz_add(R1, R5, R4);    //R1 = R5+R4

    //Block 3
    mpz_pow_ui(R4, R1, 2);  //R4 = R1^2
    mpz_neg(R5, R1);        //R5 = -R1
    mpz_add(R8, R6, R6);    //R8 = R6+R6
    mpz_add(R1, R4, R8);    //R1 = R4+R8

    mpz_mod(R1, R1, p);     //R1 mod p

    mpz_mul(R4, R2, R3);    //R4 = R2*R3
    mpz_add(R3, R4, R4);    //R3 = R4+R4

    mpz_mod(R3, R3, p);     //R3 mod p

    //Block 4
    mpz_pow_ui(R8, R7, 2);  //R6 = R2^2
    mpz_neg(R2, R8);        //R7 = -R1
    mpz_add(R8, R1, R6);    //R1 = R7+R7
    mpz_add(R4, R2, R2);    //R7 = R6+R6
    mpz_mul(R6, R5, R8);    //R6 = R1*R7
    mpz_add(R2, R6, R4);    //R1 = R5+R4
    
    mpz_mod(R2, R2, p);     //R1 mod p

    //jx3 = R1, jy3 = R2, jz3 = R3

    mpz_set(jx3, R1);
    mpz_set(jy3, R2);
    mpz_set(jz3, R3);

    printf("Resultado bloque atomico 1: Doblado Jacobiano\n");

    gmp_printf("G+2P: (%#Zd, %#Zd, %#Zd) \n", jx3, jy3, jz3);

    mpz_clear(R1);
    mpz_clear(R2);
    mpz_clear(R3);
    mpz_clear(R4);
    mpz_clear(R5);
    mpz_clear(R6);
    mpz_clear(R7);
    mpz_clear(R8);
    mpz_clear(jx3);
    mpz_clear(jy3);
    mpz_clear(jz3);
}