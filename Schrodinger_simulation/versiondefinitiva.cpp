#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "complex.h"
#define PI 3.14159265358979323846

using namespace std;

int main()
{

    int i, j, N;
    double nciclos, lambda, ko, s;
    double xo, sigma, yo, gamma;
    double norma;
    fcomplex a, b;
    vector<double> potencial, envolvente;
    vector<fcomplex> funcion, alfa, beta, be, chi, Ao;
    double modulo;

    ofstream salida1, salida2;

    salida1.open("C:/Users/sanch/Documents/SARA/PRUEBAS/2026/norma.dat");
    salida2.open("C:/Users/sanch/Documents/SARA/PRUEBAS/2026/datos.dat");

    N = 500;
    nciclos = 25.0;

    potencial.resize(N + 1);
    funcion.resize(N + 1);
    envolvente.resize(N + 1);
    chi.resize(N + 1);
    Ao.resize(N + 1);
    alfa.resize(N);
    // Cambio 1:
    // beta.resize(N);
    beta.resize(N + 1);
    // be.resize(N);
    be.resize(N + 1);

    // Caracterizar la función de onda, onda plana con amplitud gaussiana anchura sigma y centro en xo

    ko = (2.0 * PI * nciclos) / (N);
    s = (1.0) / (4.0 * ko * ko);
    xo = N / 4.0;
    sigma = N / 16.0;

    // Funcion de onda inicial

    norma = 0.0;
    funcion[0] = Complex(0.0, 0.0);
    funcion[N] = Complex(0.0, 0.0);

    for (j = 1; j < N; j++)
    {
        funcion[j] = Complex(cos(ko * j) * exp((-(j - xo) * (j - xo)) / (2.0 * sigma * sigma)), sin(ko * j) * exp((-(j - xo) * (j - xo)) / (2.0 * sigma * sigma)));
        norma += funcion[j].r * funcion[j].r + funcion[j].i * funcion[j].i;
    }

    // Normalizamos la función de onda inicial

    for (j = 1; j < N; j++)
    {
        funcion[j].r = (funcion[j].r) / (sqrt(norma));
        funcion[j].i = (funcion[j].i) / (sqrt(norma));
    }

    // Caracterizar el potencial, anchura gamma y centro en yo

    lambda = 1.0;
    yo = N / 2.0;
    gamma = N / 5.0;

    for (j = 0; j <= N; j++)
    {
        if ((j >= (yo - (gamma / 2.0))) && (j <= (yo + (gamma / 2.0))))
            potencial[j] = lambda * ko * ko;
        else
            potencial[j] = 0.0;
    }

    // Paso1: Calcular alfa que es constante en el tiempo

    for (j = 0; j <= N; j++)
    {
        Ao[j] = Complex(-2.0 - potencial[j], 2.0 / s);
    }

    a = Complex(-1.0, 0.0);

    for (j = (N - 1); j >= 0; j--)
    {
        if (j == (N - 1))
        {
            alfa[j] = Complex(0.0, 0.0);
        }
        else
        {
            alfa[j] = Cdiv(a, Cadd(Ao[j + 1], alfa[j + 1]));
        }
    }

    // Recurrencia

    for (i = 0; i <= 500; i++)
    {

        // Paso2: Calcular beta

        b = Complex(1.0, 0.0);
        for (j = 0; j <= N; j++)
        {
            be[j] = Complex((-4.0 * funcion[j].i) / s, (4.0 * funcion[j].r) / s);
        }

        beta[N - 1] = Complex(0.0, 0.0);
        // beta[N]=Complex(0.0,0.0);

        for (j = (N - 2); j >= 0; j--)
        {
            // beta[j] = Cmul(Cdiv(b,Cadd(Ao[j+1],alfa[j+1])),Csub(be[j+1],beta[j+1]));
            beta[j] = Cmul(Cmul(a, alfa[j]), Csub(be[j + 1], beta[j + 1]));
        }
        beta[0] = Complex(0.0, 0.0);

        // Paso3: Calcular chi

        chi[0] = Complex(0.0, 0.0);
        chi[N] = Complex(0.0, 0.0);

        for (j = 0; j <= (N - 2); j++)
        {
            chi[j + 1] = Cadd(Cmul(alfa[j], chi[j]), beta[j]);
        }

        // Paso4: Calcular phi y su norma

        norma = 0.0;
        funcion[0] = Complex(0.0, 0.0);
        funcion[N] = Complex(0.0, 0.0);
        envolvente[0] = funcion[0].r * funcion[0].r + funcion[0].i * funcion[0].i;
        envolvente[N] = funcion[N].r * funcion[N].r + funcion[N].i * funcion[N].i;

        salida2 << 0 << "," << envolvente[0] << "," << potencial[0] << endl;

        for (j = 1; j < N; j++)
        {
            funcion[j] = Csub(chi[j], funcion[j]);
            norma += funcion[j].r * funcion[j].r + funcion[j].i * funcion[j].i;
        }

        // Normalizamos la función de onda

        for (j = 1; j < N; j++)
        {

            envolvente[j] = funcion[j].r * funcion[j].r + funcion[j].i * funcion[j].i;
            salida2 << j << "," << envolvente[j] << "," << potencial[j] << endl;
        }

        salida2 << N << "," << envolvente[N] << "," << potencial[N] << endl;
        salida2 << endl;

        salida1 << norma << endl;
    }

    salida1.close();
    salida2.close();

    return 0;
}