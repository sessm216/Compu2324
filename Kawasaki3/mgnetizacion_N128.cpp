using namespace std; 
# include <stdio.h>
# include <iostream>
# include <cmath>
# include <fstream>
# include <vector>
#include <chrono>
#include <random>

//#define N 32
//const int N=32 ; 

void Intercambio(int reticula[][32], double T); 
//void signo_s(int (&s)[N][N],  double T);


int main()
{
    int N = 32; 
    int i,j, k; 
    //int N; N=32; 
    int reticula[32][32]; 

    double energia, p, epp, xi, T, magnetizacion, x, direction;

    double msup, minf, aux1,aux2,aux3; 
    T=0.0; 
    aux1=0.0;
    aux2=0.0; 

    magnetizacion=0.0; 
    msup=0.0; 
    minf=0.0; 

    int paso, contadorMC, iteraciones; 

    paso=N*N; 

    bool calcularmagnetizacion = true, densidad=false, simulacion =false; 


    ofstream salida1, salida2, salida3;  
    if(calcularmagnetizacion=true)     salida2.open("/home/sara/Escritorio/compu2324/Voluntarios/Kawasaki3/magnetizacion_N32_MC1000_it1000.dat"); ; 
    if(densidad==true)     salida3.open("/home/sara/Escritorio/compu2324/Voluntarios/Voluntario_Kawasaki/Lunesmanana/densidad55_10000_1.dat");
    if(simulacion==true) salida1.open("/home/sara/Escritorio/compu2324/Voluntarios/Voluntario_Kawasaki2/pruebapedro_3.dat"); 


    int l,m,posmax;
    l=N;
    m=N; 
    posmax=N*N/2; 
    int i1,i2,j1,j2; 


   for(i=0; i<N;i++)
   {
        for(j=0;j<N;j++)
        {
            if(i<N/2) 
                {reticula[i][j]=-1.0; 
                if (calcularmagnetizacion==true) minf+=reticula[i][j];
            }
            else 
            {   reticula[i][j]=1; 
                if (calcularmagnetizacion==true)msup+=reticula[i][j]; 
            }
        }
   }
    if(calcularmagnetizacion==true) 
    {   
        msup=msup/(N*N);
        minf=minf/(N*N); 
        salida2 << T << "  " << minf <<  "   " << msup << "  " << magnetizacion << endl; }


    //Paso Montecarlo: dentro de cada uno hago N*N Kawasakis (no segura)
    
    for(T=0.5;T<=4.5; T=T+0.5)
    {
        aux1=0;
        aux2=0; 
        //msup=0;
        //minf=0; 
        //magnetizacion=0; 
        energia=0.0; 
        for(contadorMC=0; contadorMC<1000; contadorMC++)
        {
            /*for(i=0; i<N; i++)
            {
            for(j=0; j<N; j++)
                {
                    if(j==N-1) 
                    {
                        salida1 << reticula[i][j] << " " << endl; 
                    }
                    else
                    { 
                    salida1 << reticula[i][j] << ","; 
                    }
                }
            }*/
            salida1 << endl;
    
            for(iteraciones=0;iteraciones<paso;iteraciones++)
            {
                Intercambio(reticula,T);
            }
        
            if(calcularmagnetizacion==true)
            {
                for(i=0; i<N;i++)
                {
                    for(j=0;j<N;j++)
                    {
                        magnetizacion+=reticula[i][j];
                        if(i<N/2) 
                        {aux1+=reticula[i][j];
                        }
                        else 
                        {
                            aux2+=reticula[i][j]; 
                        }
                    }
                }
            }
        }
        if(calcularmagnetizacion==true)
        {
            minf=aux1*1.0/(1000*N*N);
            msup=aux2*1.0/(1000*N*N); 
            magnetizacion=magnetizacion*1.0/(1000*N*N); 

            salida2 << T << "  " << minf << "   " << msup << "  " << magnetizacion << endl;
        }

    }

    //salida1.close(); 
    salida2.close(); 

return 0; 
}


void Intercambio(int reticula[][32], double T)
{
    int N=32; 

      unsigned seed1 = chrono::system_clock::now().time_since_epoch().count(); // obtain a seed from the system clock
    mt19937_64 generator(seed1);  // initialize the random number engine with seed1
    uniform_real_distribution<double> r_distribution(0., 1.); //initialize the distribution r_distribution 

    unsigned seed4 = chrono::system_clock::now().time_since_epoch().count(); // obtain a seed from the system clock
                mt19937_64 generator4(seed4);  // initialize the random number engine with seed1
                uniform_int_distribution<int> s_distribution(0,3); //initialize the distribution r_distribution

    unsigned seed2 = chrono::system_clock::now().time_since_epoch().count(); // obtain a seed from the system clock
                mt19937_64 generator2(seed2);  // initialize the random number engine with seed1
                uniform_int_distribution<int> j_distribution(0., N-1); //initialize the distribution r_distribution
    
    unsigned seed3 = chrono::system_clock::now().time_since_epoch().count(); // obtain a seed from the system clock
                mt19937_64 generator3(seed3);  // initialize the random number engine with seed1
                uniform_int_distribution<int> i_distribution(2., N-3); //initialize the distribution r_distribution

    double epp = 0.0, p=0.0 , xi=0.0, energia=0.0;

    double suma1=0.0, suma2=0.0; 
    int direction, i1,j1,i2,j2, espin1, espin2;
    
    xi=r_distribution(generator);
    direction=s_distribution(generator4);
    j1=j_distribution(generator2); 
    i1=i_distribution(generator3); 
    i2=i1;
    j1=j2; 

    //Determinar el vecino: aplicar condiciones de contorno
    if(direction==0) i2=i1+1; //Abajo
    if(direction==1) i2=i1-1; 
    if (direction == 2) j2 = (j1 + 1) % N; // Derecha
    if (direction == 3) j2 = (j1 - 1 + N) % N; // Izquierda
    espin2=reticula[i2][j2];    
    espin1=reticula[i1][j1]; 
    
    if(espin1 != espin2)
    {

    suma1=reticula[(i1+1)%N][j1]+reticula[(i1-1+N)%N][j1]+reticula[i1][j1+1]+reticula[i1][(j1-1+N)%N]-espin2; 
    suma2=reticula[(i2+1)%N][j2]+reticula[(i2-1+N)%N][j2]+reticula[i2][j2+1]+reticula[i2][(j2-1+N)%N]-espin1; 

    energia=-0.5*(espin2*(suma1+espin1)+espin1*(suma2+espin2) -(espin1*(suma1+espin2)+espin2*(suma2+espin1))); 
    epp=exp((-1.0)*energia/T);
    if(1<epp) p=1.0; 
    else p=epp; 

    if(xi<p)
    {reticula[i1][j1]=espin2; 
    reticula[i2][j2]=espin1; 
    }
    }
    


}
