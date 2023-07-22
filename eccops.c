//
//Autores: Agustín Martínez - Sergio Meirone
//Descripción: Implementación de las operaciones sobre curvas elípticas
//Curso: Criptografía I
//


#include <stdio.h>
#include <time.h>
#include <gmp.h>

void echo(mpz_t *a, mpz_t *b);

int main()
{
    //Parámetros de la curva P-192
    const char primo[] = "6277101735386680763835789423207666416083908700390324961279";
    const char EC_B[] = "64210519 e59c80e7 0fa7e9ab 72243049 feb8deec c146b9b1";
    const char EC_X1[] = "188da80e b03090f6 7cbf20eb 43a18800 f4ff0afd 82ff1012";
    const char EC_Y1[] = "07192b95 ffc8da78 631011ed 6b24cdd5 73f977a1 1e794811";
    const char SEED[] = "3045ae6f c8422f64 ed579528 d38120ea e12196d5";

    //Variables y estado aleatorio

    mpz_t p, a, b, dos, tres, seed, x1, y1, x2, y2, lambda1, inverso;
    mpz_t x3, y3;
    mpz_t jx1, jy1, jz1, alpha, beta, jx2, jy2, jz2, aux, jx3, jy3, jz3;

    gmp_randstate_t state;

    //Inicialización
    mpz_init(p);
    gmp_randinit_mt(state);

    mpz_init(a);
    mpz_init(b);
    mpz_init(tres);
    mpz_init(x1);
    mpz_init(y1);
    mpz_init(lambda1);
    mpz_init(seed);
    mpz_init(inverso);
    mpz_init(x2);
    mpz_init(y2);
    mpz_init(dos);

    mpz_init(x3);
    mpz_init(y3);

    mpz_init(jx1);
    mpz_init(jy1);
    mpz_init(jz1);
    mpz_init(alpha);
    mpz_init(beta);
    mpz_init(jx2);
    mpz_init(jy2);
    mpz_init(jz2);
    mpz_init(aux);
    mpz_init(jx3);
    mpz_init(jy3);
    mpz_init(jz3);


    //Seed del random state
    //gmp_randseed_ui(state, time(NULL));


    //Variables de tiempo para la medición
    clock_t time1, time2;

    //time1 = clock();

    // Set de parámetros de la curva: p, a, b;
    mpz_set_str(p, primo, 10);
    mpz_set_ui(dos, 2);
    mpz_set_ui(tres, 3);
    mpz_sub(a, p, tres);
    mpz_set_str(b, EC_B, 16);
    mpz_set_str(seed, SEED, 16);

    //printf("This is a test\n");
    //echo(&dos, &tres);

    //Seed del random state
    gmp_randseed(state, seed);

    //Set del punto base G: (x1,y1)
    mpz_set_str(x1, EC_X1, 16);
    mpz_set_str(y1, EC_Y1, 16);

    //Cálculo de lambda1 para operación doblado de G
    mpz_powm(lambda1, x1, dos, p);
    mpz_mul(lambda1, lambda1, tres);
    mpz_add(lambda1, lambda1, a);
    mpz_mul(inverso, y1, dos);
    mpz_invert(inverso, inverso, p);
    mpz_mul(lambda1, lambda1, inverso);
    mpz_mod(lambda1, lambda1, p);


    //Cálculo de x2
    mpz_mul(x2, lambda1, lambda1);
    mpz_sub(x2, x2, x1);
    mpz_sub(x2, x2, x1);
    mpz_mod(x2, x2, p);

    //Cálculo de y2
    mpz_sub(y2, x1, x2);
    mpz_mul(y2, y2, lambda1);
    mpz_sub(y2, y2, y1);
    mpz_mod(y2, y2, p);

    printf("-----------Parámetros de la EC-----------\n");

    gmp_printf("primo: %#Zd\n", p);
    gmp_printf("Elliptic Curve: y^2 = x^3 + %Zd *x + %Zd \n", a, b);
    gmp_printf("Punto base G: (%#Zd, %#Zd) \n", x1, y1);

    printf("-----------Doblado con fórmulas explícitas en coordenadas afín-----------\n");
    
    gmp_printf("2G = (%#Zd, %#Zd)\n", x2, y2);

    //Cálculo de G + 2G con fórmula explícita
    mpz_sub(lambda1, y1, y2);
    mpz_sub(inverso, x1, x2);
    mpz_invert(inverso, inverso, p);
    mpz_mul(lambda1, lambda1, inverso);

    //Cálculo de x3
    mpz_mul(x3, lambda1, lambda1);
    mpz_sub(x3, x3, x1);
    mpz_sub(x3, x3, x2);
    mpz_mod(x3, x3, p);

    //Cálculo de y2
    mpz_sub(y3, x1, x3);
    mpz_mul(y3, y3, lambda1);
    mpz_sub(y3, y3, y1);
    mpz_mod(y3, y3, p);

    printf("-----------Suma con fórmulas explícitas en coordenadas afín-----------\n");

    gmp_printf("G + 2G = (%#Zd, %#Zd)\n", x3, y3);


    //Set de coordenadas jacobianas punto J
    mpz_set(jx1, x1);
    mpz_set(jy1, y1);
    mpz_set_ui(jz1, 1);

    //Cálculo de 2J
    //Cálculo de alpha
    mpz_powm(alpha, jx1, dos, p);
    mpz_mul_ui(alpha, alpha, 3);
    mpz_add(alpha, alpha, a);

    //Cálculo de beta
    mpz_powm(beta, jy1, dos, p);
    mpz_mul(beta, beta, jx1);
    mpz_mul_ui(beta, beta, 4);

    //Cálculo de jx2
    mpz_powm(jx2, alpha, dos, p);
    mpz_submul(jx2, beta, dos);
    mpz_mod(jx2, jx2, p);

    //Cálculo de jy2
    mpz_sub(jy2, beta, jx2);
    mpz_mul(jy2, jy2, alpha);
    mpz_pow_ui(aux, jy1, 4);
    mpz_mul_ui(aux, aux, 8);
    mpz_sub(jy2, jy2, aux);
    mpz_mod(jy2, jy2, p);

    //Cálculo de jz2
    mpz_add(jz2, jy1, jy1);
    mpz_mod(jz2, jz2, p);

    printf("\n-----------Doblado con fórmulas explícitas en coordenadas jacobianas-----------\n");

    gmp_printf("2J = (%#Zd, %#Zd, %#Zd)\n", jx2, jy2, jz2);

    //Conversión a coordenadas afines
/*    mpz_invert(inverso, jz2, p);
    mpz_mul(inverso, inverso, inverso);
    mpz_mul(jx2, jx2, inverso);
    mpz_mod(jx2, jx2, p);
    mpz_mul(inverso, inverso, inverso);
    mpz_mul(inverso, inverso, jz2);
    mpz_mul(jy2, jy2, inverso);
    mpz_mod(jy2, jy2, p); */

    printf("\n-----------Doblado con Jacobianas en coordenadas afín-----------\n");

    gmp_printf("2J afin= (%#Zd, %#Zd)\n", jx2, jy2);


    //Cálculo de G + 2J con suma mixta

    //Cálculo de alpha
    mpz_powm(alpha, jz2, tres, p);
    mpz_mul(alpha, alpha, y1);
    mpz_sub(alpha, alpha, jy2);


    //Cálculo de beta
    mpz_powm(beta, jz2, dos, p);
    mpz_mul(beta, beta, x1);
    mpz_sub(beta, beta, jx2);

    //Cálculo de jx3
    mpz_powm(jx3, alpha, dos, p);
    mpz_powm(aux, beta, tres, p);
    mpz_submul(jx3, jx3, aux);
    mpz_powm(aux, beta, dos, p);
    mpz_mul(aux, aux, jx2);
    mpz_submul(jx3, aux, dos);
    mpz_mod(jx3, jx3, p);

    //Cálculo de jy3
    mpz_powm(aux, beta, dos, p);
    mpz_mul(jy3, jx2, aux);
    mpz_sub(jy3, jy3, jx3);
    mpz_mul(jy3, jy3, alpha);
    mpz_mul(aux, aux, beta);
    mpz_mul(aux, aux, jy2);
    mpz_sub(jy3, jy3, aux);
    mpz_mod(jy3, jy3, p);

    //Cálculo de jz3
    mpz_mul(jz3, jz2, beta);
    mpz_mod(jz3, jz3, p);

    printf("\n-----------Suma con fórmulas mixtas en coordenadas jacobianas-----------\n");

    gmp_printf("G + 2J = (%#Zd, %#Zd, %#Zd)\n", jx3, jy3, jz3);

     //Conversión a coordenadas afines
    mpz_invert(inverso, jz3, p);
    mpz_mul(inverso, inverso, inverso);
    mpz_mul(jx3, jx3, inverso);
    mpz_mod(jx3, jx3, p);
    mpz_mul(inverso, inverso, inverso);
    mpz_mul(inverso, inverso, jz2);
    mpz_mul(jy3, jy3, inverso);
    mpz_mod(jy3, jy3, p);

    printf("\n-----------Suma con fórmulas mixtas en coordenadas afín-----------\n");

    gmp_printf("G + 2J afin= (%#Zd, %#Zd)\n", jx3, jy3);



    //Liberar memoria
    mpz_clear(a);
    mpz_clear(p);
    mpz_clear(b);
    mpz_clear(x1);
    mpz_clear(y1);
    mpz_clear(x2);
    mpz_clear(y2);
    mpz_clear(seed);
    mpz_clear(lambda1);
    mpz_clear(inverso);
    mpz_clear(dos);
    mpz_clear(tres);


    gmp_randclear(state);

    return 0;
}

void echo(mpz_t *a, mpz_t *b)
{
    mpz_t c;
    mpz_init(c);
    mpz_add(c, *a, *b);
    gmp_printf("This is the result: %#Zd\n", c);
    mpz_clear(c);
}