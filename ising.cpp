// Versión original con configuración inicial alineada
using namespace std;
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <chrono>
#include <random>

//

int main()
{
    int N;
    N = 20;
    double epp;
    double T, deltae;
    T = 4.5;
    int reticula[N][N];
    double x;
    double p;
    double xi;
    double paso;
    paso = N * N;
    int iteraciones;
    int contadorMC;
    ofstream salida1, salida2;
    double aux, magnetizacion;
    aux = 0;
    magnetizacion = 0;

    salida1.open("spins_v2.dat");
    salida2.open("magnetizacion_v2.dat");

    deltae = 0.0;
    int i, j, k;
    i = 1.0;
    j = 1.0;
    k = 0;
    contadorMC = 0;

    unsigned seed1 = chrono::system_clock::now().time_since_epoch().count(); // obtain a seed from the system clock
    mt19937_64 generator(seed1);                                             // initialize the random number engine with seed1
    uniform_real_distribution<double> r_distribution(0., 1.);                // initialize the distribution r_distribution

    /// Configuración inicial:
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            reticula[i][j] = 1;
        }
    }

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            if (j == N - 1)
            {
                salida1 << reticula[i][j] << " " << endl;
            }
            else
            {
                salida1 << reticula[i][j] << ",";
            }
        }
    }

    // for(T=0.001;T<=5.2;T+0.5)
    //{
    salida1 << "\n";

    for (contadorMC = 0; contadorMC < 10000; contadorMC++)
    {
        for (iteraciones = 0; iteraciones < paso; iteraciones++)
        {
            // Elegir punto inicial al azar:
            // Genero dos enteros aleatorios para seleccionar una celda al azar:
            unsigned seed2 = chrono::system_clock::now().time_since_epoch().count(); // obtain a seed from the system clock
            mt19937_64 generator2(seed2);                                            // initialize the random number engine with seed1
            uniform_int_distribution<int> j_distribution(0., N - 1);                 // initialize the distribution r_distribution

            i = j_distribution(generator2);

            j = j_distribution(generator2);

            /// Calcular exponencial y calcular DeltaE:
            if (i == 0)
            {
                if (j == 0)
                {
                    deltae = 2 * reticula[i][j] * (reticula[i + 1][j] + reticula[N - 1][j] + reticula[i][j + 1] + reticula[i][N - 1]);
                }
                if (j == N - 1)
                {
                    deltae = 2 * reticula[i][j] * (reticula[i + 1][j] + reticula[N - 1][j] + reticula[i][j - 1] + reticula[i][0]);
                }
            }
            else if (i == N - 1)
            {
                if (j == 0)
                {
                    deltae = 2 * reticula[i][j] * (reticula[0][j] + reticula[i - 1][j] + reticula[i][j + 1] + reticula[i][N - 1]);
                }
                if (j == N - 1)
                    ;
                {
                    deltae = 2 * reticula[i][j] * (reticula[0][j] + reticula[i - 1][j] + reticula[i][j - 1] + reticula[i][0]);
                }
            }
            else
            {
                deltae = 2 * reticula[i][j] * (reticula[i + 1][j] + reticula[i - 1][j] + reticula[i][j + 1] + reticula[i][j - 1]);
            }

            epp = exp((-1.0) * deltae / T);
            if (1 < epp)
            {
                p = 1.0;
            }
            else
            {
                p = epp;
            }

            xi = r_distribution(generator);
            if (xi < p)
            {
                reticula[i][j] = reticula[i][j] * (-1);
            }
        }

        for (i = 0; i < N; i++)
        {
            for (j = 0; j < N; j++)
            {
                if (j == N - 1)
                {
                    salida1 << reticula[i][j] << " " << endl;
                }
                else
                {
                    salida1 << reticula[i][j] << ",";
                }
            }
            aux = aux + reticula[i][j];
        }
        if (aux < 0)
            aux = aux * (-1);
        magnetizacion = magnetizacion + aux;
        aux = 0;
        salida1 << endl;
        k = k + 1;
    }
    //}
    magnetizacion = magnetizacion * 1.0 / (10000 * N * N);
    salida2 << T << "  " << magnetizacion << endl;
    magnetizacion = 0;
    aux = 0;

    salida1.close();
    salida2.close();

    return 0;
}