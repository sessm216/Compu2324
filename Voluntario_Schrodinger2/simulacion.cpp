//Código que calcula el coeficiete de transmisión para un  valor del parámetro lambda que parametriza el potencial V(x) y diferentes particiones N

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib> 
#include <ctime>   
#include <fstream>
#define PI 3.14159265358979323846

#include "complex.h"
using namespace std;

int main()
{
    int i, j, experimento, mT = 0, m = 1000; // número de experimentos
    double nciclos, ko, s;
    double xo, sigma, yo, gamma;
    double ProbD, ProbIz, k, q;
    int nD = 100; // intervalo de tiempo para aplicar detectores
    fcomplex a, b;
    vector<double> potencial;
    vector<fcomplex> funcion, alfa, beta, be, chi, Ao;
    ofstream salida0;
    
    salida0.open("/home/sara/Escritorio/compu2324/Voluntarios/Voluntario_Schrodinger2/transmision_50.dat");
    

    srand(time(NULL)); // Inicializar la semilla para rand() con el tiempo actual

    nciclos=50.0;
    double lambda= 5.0; //Es necesario ir cambiando el valor de lambda
    
    for (int N : {250, 500, 750, 1000, 2000})
    {
        mT = 0; // Reiniciar el contador de detecciones para cada valor de N

        potencial.resize(N + 1);
        funcion.resize(N + 1);
        chi.resize(N + 1);
        Ao.resize(N + 1);
        alfa.resize(N);
        beta.resize(N);
        be.resize(N);

        // Caracterizar la función de onda, onda plana con amplitud gaussiana anchura sigma y centro en xo
        ko = (2.0 * PI * nciclos) / N;
        s = 1.0 / (4.0 * ko * ko);
        xo = N / 4.0;
        sigma = N / 16.0;

        // Caracterizar el potencial, anchura gamma y centro en yo
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

        // Realizar m experimentos
        for (experimento = 0; experimento < m; experimento++)
        {

            bool particulaDetectada = false;

            // Inicializar la función de onda
            q = 0.0;
            funcion[0] = Complex(0.0, 0.0);
            funcion[N] = Complex(0.0, 0.0);

            for (j = 1; j < N; j++)
            {
                funcion[j] = Complex(cos(ko * j) * exp((-(j - xo) * (j - xo)) / (2.0 * sigma * sigma)), sin(ko * j) * exp((-(j - xo) * (j - xo)) / (2.0 * sigma * sigma)));
                q += funcion[j].r * funcion[j].r + funcion[j].i * funcion[j].i;
            }
            
            q = sqrt(1.0 / q);

            // Normalizar
            for (j = 0; j <= N; j++)
            {
                funcion[j].r *= q;
                funcion[j].i *= q;
            }

             while (!particulaDetectada)
            {
            
              // Evolucionar nD pasos
              for (i = 0; i <= nD; i++)
              {
                // Paso2: Calcular beta
                b = Complex(1.0, 0.0);

                for (j = (N - 1); j >= 0; j--)
                {
                    if (j == (N - 1))
                    {
                        beta[j] = Complex(0.0, 0.0);
                    }
                    else
                    {
                        be[j] = Complex((-4.0 * funcion[j].i) / s, (4.0 * funcion[j].r) / s);
                        beta[j] = Cmul(Cdiv(b, Cadd(Ao[j + 1], alfa[j + 1])), Csub(be[j + 1], beta[j + 1]));
                    }
                }

                // Paso3: Calcular chi
                chi[0] = Complex(0.0, 0.0);
                chi[N] = Complex(0.0, 0.0);

                for (j = 0; j <= (N - 2); j++)
                {
                    chi[j + 1] = Cadd(Cmul(alfa[j], chi[j]), beta[j]);
                }

                // Paso4: Calcular phi 
                funcion[0] = Complex(0.0, 0.0);
                funcion[N] = Complex(0.0, 0.0);

                for (j = 1; j < N; j++)
                {
                    funcion[j] = Csub(chi[j], funcion[j]);
                }

              }

            // Calcular PD
            ProbD = 0.0;
            for (j = 4 * N / 5; j <= N; j++)
            {
                ProbD += funcion[j].r * funcion[j].r + funcion[j].i * funcion[j].i;
            }


            // Generar número aleatorio x
            double x = (double)rand() / RAND_MAX;

            if (x < ProbD)
            {
                mT++;
                particulaDetectada = true;
                break;
            }

            // Proyección si no se detecta
            for (j = 4 * N / 5; j <= N; j++)
            {
                funcion[j] = Complex(0.0, 0.0);
            }

            // Calcular norma
            k = 0.0;
            for (j = 0; j <= N; j++)
            {
                k += funcion[j].r * funcion[j].r + funcion[j].i * funcion[j].i;
            }

            k = sqrt(1.0 / k);

            // Normalizar
            for (j = 0; j <= N; j++)
            {
                funcion[j].r *= k;
                funcion[j].i *= k;
            }

            // Calcular PI
            ProbIz = 0.0;
            for (j = 0; j <= N / 5; j++)
            {
                ProbIz += funcion[j].r * funcion[j].r + funcion[j].i * funcion[j].i;
            }

            // Generar número aleatorio x
            x = (double)rand() / RAND_MAX;

            if (x < ProbIz)
            {
                particulaDetectada = true;
                break;
            }

            // Proyección si no se detecta
            for (j = 0; j <= N / 5; j++)
            {
                funcion[j] = Complex(0.0, 0.0);
            }

            // Calcular norma
            k = 0.0;
            for (j = 0; j <= N; j++)
            {
                k += funcion[j].r * funcion[j].r + funcion[j].i * funcion[j].i;
            }

            k = sqrt(1.0 / k);

            // Normalizar
            for (j = 0; j <= N; j++)
            {
                funcion[j].r *= k;
                funcion[j].i *= k;
            }

            
          }

            
        }

            salida0 << lambda << " " << N << " " << (double)mT / m << endl;
        }


    salida0.close();

    return 0;
}
