//
//Autores: Agustín Martínez - Sergio Meirone
//Descripción: Implementación de las operaciones sobre curvas elípticas
//Curso: Criptografía I
//


#include <stdio.h>
#include <time.h>
#include <gmp.h>

void dobladoAfin(mpz_t x1, mpz_t y1, mpz_t p, mpz_t a, mpz_t x2, mpz_t y2);
void sumaAfin(mpz_t x1, mpz_t y1, mpz_t x2, mpz_t y2, mpz_t x3, mpz_t y3, mpz_t p);

void dobladoJacobiano(mpz_t jx1, mpz_t jy1, mpz_t jz1, mpz_t jx2,mpz_t jy2, mpz_t jz2, mpz_t a, mpz_t p);
void sumaMixta(mpz_t x1, mpz_t y1, mpz_t jx2, mpz_t jy2, mpz_t jz2, mpz_t jx3, mpz_t jy3, mpz_t jz3, mpz_t p);

void AtomicBlockJacobianDoubling(mpz_t jx1, mpz_t jy1, mpz_t jz1, mpz_t p);
void AtomicBlockMixedAddition(mpz_t jx1, mpz_t jy1, mpz_t jz1, mpz_t x1, mpz_t y1, mpz_t p, mpz_t inverso);



int main()
{
    //Parámetros de la curva P-192
    const char primo[] = "6277101735386680763835789423207666416083908700390324961279";
    const char EC_B[] = "64210519 e59c80e7 0fa7e9ab 72243049 feb8deec c146b9b1";
    const char EC_X1[] = "188da80e b03090f6 7cbf20eb 43a18800 f4ff0afd 82ff1012";
    const char EC_Y1[] = "07192b95 ffc8da78 631011ed 6b24cdd5 73f977a1 1e794811";

    //Variables
    mpz_t p, a, b, dos, tres, seed, x1, y1, x2, y2, lambda1, inverso;
    mpz_t x3, y3;
    mpz_t jx1, jy1, jz1, alpha, beta, jx2, jy2, jz2, aux, jx3, jy3, jz3;


    //Inicialización
    mpz_init(p);

    mpz_init(a);
    mpz_init(b);
    mpz_init(tres);
    mpz_init(x1);
    mpz_init(y1);
    mpz_init(lambda1);
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



    //Variables de tiempo para la medición
    clock_t time1, time2;

    // Set de parámetros de la curva: p, a, b;
    mpz_set_str(p, primo, 10);
    mpz_set_ui(dos, 2);
    mpz_set_ui(tres, 3);
    mpz_sub(a, p, tres);
    mpz_set_str(b, EC_B, 16);

    //Set del punto base G: (x1,y1)
    mpz_set_str(x1, EC_X1, 16);
    mpz_set_str(y1, EC_Y1, 16);

    printf("-----------Parámetros de la EC-----------\n");
    gmp_printf("primo: %#Zd\n", p);
    gmp_printf("Elliptic Curve: y^2 = x^3 + %Zd *x + %Zd \n", a, b);
    gmp_printf("Punto base G: (%#Zd, %#Zd) \n", x1, y1);

    time1 = clock();
    //Doblado de punto G en coordenadas afines
    dobladoAfin(x1, y1, p, a, x2, y2);
    time2 = clock();

    printf("-----------Doblado con fórmulas explícitas en coordenadas afín-----------\n");
    gmp_printf("2G = (%#Zd, %#Zd)\n", x2, y2);
    printf("Tiempo de ejecución: %f\n", (double)(time2 - time1) / CLOCKS_PER_SEC);
    
    time1 = clock();
    //Suma coordenadas afines G+2G
    sumaAfin(x1, y1, x2, y2, x3, y3, p);
    time2 = clock();

    printf("-----------Suma con fórmulas explícitas en coordenadas afín-----------\n");
    gmp_printf("G + 2G = (%#Zd, %#Zd)\n", x3, y3);
    printf("Tiempo de ejecución: %f\n", (double)(time2 - time1) / CLOCKS_PER_SEC);


    //Set de coordenadas jacobianas punto J
    mpz_set(jx1, x1);
    mpz_set(jy1, y1);
    mpz_set_ui(jz1, 1);

    time1 = clock();
    //Doblado 2G en coordenadas jacobianas
    dobladoJacobiano(jx1, jy1, jz1, jx2, jy2, jz2, a, p);
    time2 = clock();

    printf("\n-----------Doblado con fórmulas explícitas en coordenadas jacobianas-----------\n");
    gmp_printf("2J = (%#Zd, %#Zd, %#Zd)\n", jx2, jy2, jz2);
    printf("Tiempo de ejecución: %f\n", (double)(time2 - time1) / CLOCKS_PER_SEC);


    //Conversión a coordenadas afines
/*  mpz_invert(inverso, jz2, p);
    mpz_mul(inverso, inverso, inverso);
    mpz_mul(jx2, jx2, inverso);
    mpz_mod(jx2, jx2, p);
    mpz_mul(inverso, inverso, inverso);
    mpz_mul(inverso, inverso, jz2);
    mpz_mul(jy2, jy2, inverso);
    mpz_mod(jy2, jy2, p); */

    printf("\n-----------Doblado con Jacobianas en coordenadas afín-----------\n");
    gmp_printf("2J afin= (%#Zd, %#Zd)\n", jx2, jy2);

    time1 = clock();
    //Cálculo de G + 2J con suma mixta
    sumaMixta(x1, y1, jx2, jy2, jz2, jx3, jy3, jz3, p);
    time2 = clock();

    printf("\n-----------Suma con fórmulas mixtas en coordenadas jacobianas-----------\n");
    gmp_printf("G + 2J = (%#Zd, %#Zd, %#Zd)\n", jx3, jy3, jz3);
    printf("Tiempo de ejecución: %f\n", (double)(time2 - time1) / CLOCKS_PER_SEC);

     //Conversión a coordenadas afines
    mpz_invert(inverso, jz3, p);
    mpz_mul(inverso, inverso, inverso);
    mpz_mul(jx3, jx3, inverso);
    mpz_mod(jx3, jx3, p);
    mpz_mul(inverso, inverso, inverso);
    mpz_mul(inverso, inverso, jz3);
    mpz_mul(jy3, jy3, inverso);
    mpz_mod(jy3, jy3, p);

    printf("\n-----------Suma con fórmulas mixtas en coordenadas afín-----------\n");
    gmp_printf("G + 2J afin= (%#Zd, %#Zd)\n", jx3, jy3);


    //Doblado con bloques atómicos 2G
    time1 = clock();
    AtomicBlockJacobianDoubling(jx1, jy1, jz1, p);
    time2 = clock();
    printf("Tiempo de ejecución: %f\n", (double)(time2 - time1) / CLOCKS_PER_SEC);

    //Suma mixta con bloques atómicos G + 2G
    time1 = clock();
    AtomicBlockMixedAddition(jx2, jy2, jz2, x1, y1, p, inverso);
    time2 = clock();
    printf("Tiempo de ejecución: %f\n", (double)(time2 - time1) / CLOCKS_PER_SEC);



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

    return 0;
}

void dobladoAfin(mpz_t x1, mpz_t y1, mpz_t p, mpz_t a, mpz_t x2, mpz_t y2)
{
    //Declaración de variables
    mpz_t lambda1, inverso;

    //Inicialización de variables
    mpz_init(lambda1);
    mpz_init(inverso);

    //Cálculo de lambda1 para operación doblado de G
    mpz_powm_ui(lambda1, x1, 2, p);
    mpz_mul_ui(lambda1, lambda1, 3);
    mpz_add(lambda1, lambda1, a);
    mpz_mul_ui(inverso, y1, 2);
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

    //Liberar la memoria
    mpz_clear(lambda1);
    mpz_clear(inverso);
}

void sumaAfin(mpz_t x1, mpz_t y1, mpz_t x2, mpz_t y2, mpz_t x3, mpz_t y3, mpz_t p)
{
    //Declaración de variables
    mpz_t lambda1, inverso;

    //Inicialización de variables
    mpz_init(lambda1);
    mpz_init(inverso);

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

    //Liberar la memoria
    mpz_clear(lambda1);
    mpz_clear(inverso);
}

void dobladoJacobiano(mpz_t jx1, mpz_t jy1, mpz_t jz1, mpz_t jx2,mpz_t jy2, mpz_t jz2, mpz_t a, mpz_t p)
{
    //Declaración de variables
    mpz_t alpha, beta, aux;

    //Inicialización de variables
    mpz_init(alpha);
    mpz_init(beta);
    mpz_init(aux);

    //Cálculo de alpha
    mpz_powm_ui(alpha, jx1, 2, p);
    mpz_mul_ui(alpha, alpha, 3);
    mpz_add(alpha, alpha, a);

    //Cálculo de beta
    mpz_powm_ui(beta, jy1, 2, p);
    mpz_mul(beta, beta, jx1);
    mpz_mul_ui(beta, beta, 4);

    //Cálculo de jx2
    mpz_powm_ui(jx2, alpha, 2, p);
    mpz_submul_ui(jx2, beta, 2);
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

    //Liberar la memoria
    mpz_clear(alpha);
    mpz_clear(beta);
    mpz_clear(aux);

}

void sumaMixta(mpz_t x1, mpz_t y1, mpz_t jx2, mpz_t jy2, mpz_t jz2, mpz_t jx3, mpz_t jy3, mpz_t jz3, mpz_t p)
{
    //Declaración de variables
    mpz_t alpha, beta, aux;

    //Inicialización de variables
    mpz_init(alpha);
    mpz_init(beta);
    mpz_init(aux);

    //Cálculo de alpha
    mpz_powm_ui(alpha, jz2, 2, p);
    mpz_mul(alpha, alpha, y1);
    mpz_sub(alpha, alpha, jy2);

    //Cálculo de alpha
    mpz_powm_ui(alpha, jz2, 3, p);
    mpz_mul(alpha, alpha, y1);
    mpz_sub(alpha, alpha, jy2);


    //Cálculo de beta
    mpz_powm_ui(beta, jz2, 2, p);
    mpz_mul(beta, beta, x1);
    mpz_sub(beta, beta, jx2);

    //Cálculo de jx3
    mpz_powm_ui(jx3, alpha, 2, p); //Alpha al cuadrado
    mpz_powm_ui(aux, beta, 3, p); //Beta al cubo
    mpz_sub(jx3, jx3, aux); //Alpha al cuadrado - Beta al cubo
    mpz_powm_ui(aux, beta, 2, p); //Beta al cuadrado
    mpz_mul(aux, aux, jx2); //Beta al cuadrado * jx2
    mpz_submul_ui(jx3, aux, 2); //Alpha al cuadrado - Beta al cubo - 2*Beta al cuadrado * jx2
    mpz_mod(jx3, jx3, p); //Modulo p

    //Cálculo de jy3
    mpz_powm_ui(aux, beta, 2, p); //Beta al cuadrado
    mpz_mul(jy3, jx2, aux); //jx2 * Beta al cuadrado
    mpz_sub(jy3, jy3, jx3); //jx2 * Beta al cuadrado - jx3
    mpz_mul(jy3, jy3, alpha);
    mpz_mul(aux, aux, beta); 
    mpz_mul(aux, aux, jy2);
    mpz_sub(jy3, jy3, aux);
    mpz_mod(jy3, jy3, p);

    //Cálculo de jz3
    mpz_mul(jz3, jz2, beta);
    mpz_mod(jz3, jz3, p);

    //Liberar la memoria
    mpz_clear(alpha);
    mpz_clear(beta);
    mpz_clear(aux);
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
    mpz_init(jx3);
    mpz_init(jy3);
    mpz_init(jz3);


    //R1 = jx1, R2 = jy1, R3 = jz1
    mpz_set(R1, jx1);
    mpz_set(R2, jy1);
    mpz_set(R3, jz1);

    //Block 1
    mpz_pow_ui(R4, R3, 2);  //R4 = R3^2     S
    mpz_neg(R5, R4);        //R5 = -R4      N
    mpz_add(R6, R1, R4);    //R6 = R1+R4    A
    mpz_add(R4, R1, R5);    //R4 = R1+R5    A
    mpz_mul(R5, R6, R4);    //R5 = R6*R4    M
    mpz_add(R4, R5, R5);    //R4 = R5+R5    A

    //Block 2
    mpz_pow_ui(R6, R2, 2);  //R6 = R2^2     S
    mpz_neg(R7, R1);        //R7 = -R1      N
    mpz_add(R1, R7, R7);    //R1 = R7+R7    A
    mpz_add(R7, R6, R6);    //R7 = R6+R6    A
    mpz_mul(R6, R1, R7);    //R6 = R1*R7    M
    mpz_add(R1, R5, R4);    //R1 = R5+R4    A

    //Block 3
    mpz_pow_ui(R4, R1, 2);  //R4 = R1^2     S
    mpz_neg(R5, R1);        //R5 = -R1      N
    mpz_add(R8, R6, R6);    //R8 = R6+R6    A
    mpz_add(R1, R4, R8);    //R1 = R4+R8    A

    mpz_mod(R1, R1, p);     //R1 mod p

    mpz_mul(R4, R2, R3);    //R4 = R2*R3    M
    mpz_add(R3, R4, R4);    //R3 = R4+R4    A

    mpz_mod(R3, R3, p);     //R3 mod p

    //Block 4
    mpz_pow_ui(R8, R7, 2);  //R6 = R2^2     S
    mpz_neg(R2, R8);        //R7 = -R1      N
    mpz_add(R8, R1, R6);    //R1 = R7+R7    A
    mpz_add(R4, R2, R2);    //R7 = R6+R6    A
    mpz_mul(R6, R5, R8);    //R6 = R1*R7    M
    mpz_add(R2, R6, R4);    //R1 = R5+R4    A
    
    mpz_mod(R2, R2, p);     //R1 mod p

    //jx3 = R1, jy3 = R2, jz3 = R3

    mpz_set(jx3, R1);
    mpz_set(jy3, R2);
    mpz_set(jz3, R3);

    printf("Resultado bloque atomico 1: Doblado Jacobiano\n");

    gmp_printf("(%#Zd, %#Zd, %#Zd) \n", jx3, jy3, jz3);

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

void AtomicBlockMixedAddition(mpz_t jx1, mpz_t jy1, mpz_t jz1, mpz_t x1, mpz_t y1, mpz_t p, mpz_t inverso){

    //mpz R1 a R12
    mpz_t R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, jx3, jy3, jz3;

    mpz_init(R1);
    mpz_init(R2);
    mpz_init(R3);
    mpz_init(R4);
    mpz_init(R5);
    mpz_init(R6);
    mpz_init(R7);
    mpz_init(R8);
    mpz_init(R9);
    mpz_init(R10);
    mpz_init(R11);
    mpz_init(R12);
    mpz_init(jx3);
    mpz_init(jy3);
    mpz_init(jz3);

        //Asignación
    mpz_set(R1, jx1);
    mpz_set(R2, jy1);
    mpz_set(R3, jz1);
    mpz_set(R4, x1);
    mpz_set(R5, y1);

    //Bloque 1
    mpz_pow_ui(R6, R3, 2);  //S
    mpz_neg(R1, R1);        //N
    mpz_add(R2, R2, R2);    //A
    mpz_add(R8, R1, R1);    //A
    mpz_mul(R4, R6, R4);    //M
    mpz_add(R1, R4, R1);    //A

    //Bloque 2
    mpz_pow_ui(R4, R1, 2);  //S
    mpz_neg(R9, R2);        //N
    mpz_add(R10, R6, R4);   //A
    mpz_add(R11, R3, R1);   //A
    mpz_mul(R6, R6, R3);    //M
    mpz_add(R3, R8, R8);    //A

    //Bloque 3
    mpz_pow_ui(R8, R4, 2);  //S
    mpz_neg(R10, R10);      //N
    mpz_add(R12, R4, R1);   //A
    mpz_add(R1, R8, R4);    //A
    mpz_mul(R6, R6, R5);    //M
    mpz_add(R9, R6, R9);    //A

    //Bloque 4
    mpz_pow_ui(R8, R12, 2); //S
    mpz_neg(R1, R1);        //N
    mpz_add(R1, R8, R1);    //A
    mpz_add(R9, R9, R6);    //A
    mpz_mul(R6, R3, R4);    //M
    mpz_add(R3, R6, R6);    //A

    //Bloque 5
    mpz_pow_ui(R4, R9, 2); //S
    mpz_neg(R1, R1);        //N
    mpz_add(R8, R1, R1); //A
    mpz_add(R3, R4, R3); //A
    mpz_mul(R4, R2, R8); //M
    mpz_add(R1, R3, R8); //A

    mpz_mod(R1, R1, p);

    //Bloque 6
    mpz_pow_ui(R3, R11, 2); //S
    mpz_neg(R2, R9);    

    mpz_add(R7, R6, R1); //A
    mpz_add(R3, R3, R10); //A

    mpz_mod(R3, R3, p);

    mpz_mul(R8, R2, R7); //M
    mpz_add(R2, R8, R4); //A

    mpz_mod(R2, R2, p);

    //jx3 = R1, jy3 = R2, jz3 = R3
    mpz_set(jx3, R1);
    mpz_set(jy3, R2);
    mpz_set(jz3, R3);

    printf("Resultado bloque atomico 2: Suma mixta en coordenadas jacobianas\n");

    gmp_printf("G+2P: (%#Zd, %#Zd, %#Zd) \n", jx3, jy3, jz3);

    //Conversión a coordenadas afines
    mpz_invert(inverso, jz3, p);
    mpz_mul(inverso, inverso, inverso);
    mpz_mul(jx3, jx3, inverso);
    mpz_mod(jx3, jx3, p);
    mpz_mul(inverso, inverso, inverso);
    mpz_mul(inverso, inverso, jz3);
    mpz_mul(jy3, jy3, inverso);
    mpz_mod(jy3, jy3, p);

    printf("Resultado bloque atomico 2: Suma mixta en coordenadas afines\n");

    gmp_printf("G+2P: (%#Zd, %#Zd) \n", jx3, jy3);

    mpz_clear(R1);
    mpz_clear(R2);
    mpz_clear(R3);
    mpz_clear(R4);
    mpz_clear(R5);
    mpz_clear(R6);
    mpz_clear(R7);
    mpz_clear(R8);
    mpz_clear(R9);
    mpz_clear(R10);
    mpz_clear(R11);
    mpz_clear(R12);
    mpz_clear(jx3);
    mpz_clear(jy3);
    mpz_clear(jz3);

}