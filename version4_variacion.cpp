#include <iostream>
#include <fstream>
#include <cmath>
#define G 6.67E-11
#define Mt 5.9736E+24
#define Ml 0.07349E+24
#define RT 6.378160E+6
#define dTl 3.844E+8
#define w 2.6617E-6
#define PI 3.14159265

using namespace std;

void AlgoritmoRunge(double f[5], double h, double delta, double mu);

int main(void)
{

  int j;
  int h;
  h = 10.0;

  // Variables del algoritmo:

  double delta, mu;
  delta = (G * Mt) / (dTl * dTl * dTl);
  mu = (Ml) / (Mt);

  // Variables de la energía:
  double Ho, H, rL0, rL, H2;

  // POsiciones delcohete y la luna:
  double x0, y0, xLuna0, yLuna0; // COndiciones iniciales
  double x, y, xLuna, yLuna;

  // Variables de la trayectoria:
  double r0, phi0, v0, pr0, pphi0, theta0, t0;
  double f[5];
  ofstream salida1, salida2, salida3;

  salida1.open("/home/sara/Escritorio/compu2324/cohete/naveangulo_6.dat");
  salida2.open("/home/sara/Escritorio/compu2324/cohete/Hnormalizado_6.dat");
  salida3.open("/home/sara/Escritorio/compu2324/cohete/Hangulo_6.dat");

  // Tenemos 4 ecuaciones diferenciales, f1 es r, f2 es phi,f3 es pr,f4 es pphi,f5 es t

  // Paso1:Dar condiciones iniciales

  r0 = (RT / dTl);    // Radio de la Tierra
  v0 = (11200 / dTl); // Velocidad de escape de la Tierra 11.2 Km/s
  phi0 = 0.44;        // Polo Norte
  theta0 = 0.44;      // Angulo de lanzamiento 52.5
  t0 = 0.0;           // Instante inicial

  // Paso1.1:Calcular momentos iniciales, consideramos masa de la nave igual a la unidad además ya están reescalados

  pr0 = v0 * cos(theta0 - phi0);
  pphi0 = v0 * r0 * sin(theta0 - phi0);

  f[0] = r0;
  f[1] = phi0;
  f[2] = pr0;
  f[3] = pphi0;
  f[4] = t0;

  // Parámetros de las ecuaciones de movimiento

  // Paso4:Calculamos las posiciones iniciales del cohete y de la luna

  x0 = f[0] * cos(f[1]);
  y0 = f[0] * sin(f[1]);
  xLuna0 = 1.0; // Distancia dTl
  yLuna0 = 0.0;
  rL0 = sqrt(dTl * dTl + pow(f[0], 2) - 2 * f[0] * dTl * cos(f[1] - w * f[4])) / dTl;

  Ho = pow(f[2], 2) * dTl * dTl / (2) + pow(f[3], 2) * dTl * dTl / (2.0 * pow(f[0], 2)) - G * Mt / (f[0] * dTl) - G * Ml / (rL0)-w * f[3] * dTl * dTl;

  salida2 << Ho / Ho << endl;
  salida3 << Ho << endl;

  salida1 << x0 << "," << y0 << endl;
  salida1 << xLuna0 << "," << yLuna0 << endl;
  salida1 << endl;

  // Algoritmo Runge-Kutta

  for (j = 0; j < 1E+6; j++)
  {

    AlgoritmoRunge(f, h, delta, mu);

    // Paso4:Calculamos las posiciones del cohete y de la luna

    x = f[0] * cos(f[1]);
    y = f[0] * sin(f[1]);
    xLuna = cos(w * f[4]);
    yLuna = sin(w * f[4]);
    rL = sqrt(dTl * dTl + pow(f[0], 2) - 2 * f[0] * dTl * cos(f[1] - w * f[4])) / dTl;

    H = (pow(f[2], 2) * dTl * dTl / (2) + pow(f[3], 2) * dTl * dTl / (2.0 * pow(f[0], 2)) - G * Mt / (f[0] * dTl) - G * Ml / (rL)-w * f[3] * dTl * dTl) / Ho;
    H2 = (pow(f[2], 2) * dTl * dTl / (2) + pow(f[3], 2) * dTl * dTl / (2.0 * pow(f[0], 2)) - G * Mt / (f[0] * dTl) - G * Ml / (rL)-w * f[3] * dTl * dTl);

    if (j % 60 == 0)
    {
      salida1 << x << "," << y << endl;
      salida1 << xLuna << "," << yLuna << endl;
      salida1 << endl;
    }

    salida2 << H << endl;
    salida3 << H2 << endl;
  }

  salida1.close();
  salida2.close();
  salida3.close();

  return 0;
}

void AlgoritmoRunge(double f[5], double h, double delta, double mu)
{

  int i;
  double aux1, aux2, aux3, aux4;
  double k1[4], k2[4], k3[4], k4[4];

  // Paso2:Calcular vectores k

  aux1 = sqrt(1 + pow(f[0], 2) - 2 * f[0] * cos(f[1] - w * f[4]));

  k1[0] = h * f[2];
  k1[1] = h * f[3] / (pow(f[0], 2));
  k1[2] = h * ((pow(f[3], 2) / pow(f[0], 3)) - delta * (1 / (pow(f[0], 2)) + mu * (f[0] - cos(f[1] - w * f[4])) / pow(aux1, 3)));
  k1[3] = h * (-delta * mu * f[0] * sin(f[1] - w * f[4]) / (pow(aux1, 3)));

  aux2 = sqrt(1 + pow(f[0] + k1[0] / 2, 2) - 2 * (f[0] + k1[0] / 2) * cos(f[1] + k1[1] / 2 - w * (f[4] + h / 2)));

  k2[0] = h * (f[2] + k1[2] / 2.0);
  k2[1] = h * (f[3] + k1[3] / 2.0) / (pow(f[0] + k1[0] / 2.0, 2));
  k2[2] = h * (pow(f[3] + k1[3] / 2, 2) / (pow(f[0] + k1[0] / 2, 3)) - delta * (1 / (pow(f[0] + k1[0] / 2, 2)) + mu * (f[0] + k1[0] / 2 - cos(f[1] + k1[1] / 2 - w * (f[4] + h / 2))) / pow(aux2, 3)));
  k2[3] = h * (-delta * mu * (f[0] + k1[0] / 2) * sin(f[1] + k1[1] / 2 - w * (f[4] + h / 2)) / (pow(aux2, 3)));

  aux3 = sqrt(1 + pow(f[0] + k2[0] / 2, 2) - 2 * (f[0] + k2[0] / 2) * cos(f[1] + k2[1] / 2 - w * (f[4] + h / 2)));

  k3[0] = h * (f[2] + k2[2] / 2.0);
  k3[1] = h * (f[3] + k2[3] / 2.0) / (pow(f[0] + k2[0] / 2.0, 2));
  k3[2] = h * (pow(f[3] + k2[3] / 2, 2) / (pow(f[0] + k2[0] / 2, 3)) - delta * (1 / (pow(f[0] + k2[0] / 2, 2)) + mu * (f[0] + k2[0] / 2 - cos(f[1] + k2[1] / 2 - w * (f[4] + h / 2))) / pow(aux3, 3)));
  k3[3] = h * (-delta * mu * (f[0] + k2[0] / 2) * sin(f[1] + k2[1] / 2 - w * (f[4] + h / 2)) / (pow(aux3, 3)));

  aux4 = sqrt(1 + pow(f[0] + k3[0], 2) - 2 * (f[0] + k3[0]) * cos(f[1] + k3[1] - w * (f[4] + h)));

  k4[0] = h * (f[2] + k3[2]);
  k4[1] = h * (f[3] + k3[3]) / (pow(f[0] + k3[0], 2));
  k4[2] = h * (pow(f[3] + k3[3], 2) / (pow(f[0] + k3[0], 3)) - delta * (1 / (pow(f[0] + k3[0], 2)) + mu * (f[0] + k3[0] - cos(f[1] + k3[1] - w * (f[4] + h))) / pow(aux4, 3)));
  k4[3] = h * (-delta * mu * (f[0] + k3[0]) * sin(f[1] + k3[1] - w * (f[4] + h)) / (pow(aux4, 3)));

  // Paso3:Modificar vector f

  for (i = 0; i < 4; i++)
  {
    f[i] = f[i] + (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6.0;
  }
  f[4] = f[4] + h;

  return;
}