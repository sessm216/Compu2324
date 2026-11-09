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
    int pasoMontecarlo=1000; 
    int paso, contadorMC, iteraciones; 
    paso=N*N; 

    double energia, p, epp, xi, T, magnetizacion, x, direction;
    double energiamedia = 0.0, energiacuadrado = 0.0, calor=0.0, chisup=0.0, chiinf=0.0; 

    double msup, minf, aux1,aux2,aux3, msupcuadrado, minfcuadrado, aux4; 

    T=1.5; 
    aux1=0.0;
    aux2=0.0; 

    magnetizacion=0.0; 
    msup=0.0; 
    minf=0.0; 

    int l,m,posmax;
    l=N;
    m=N; 
    posmax=N*N/2; 
    int i1,i2,j1,j2; 

    bool calcularmagnetizacion, densidadcalcular, simulacion, calcularenergia, calcularenergiamedia, calcularcalorespecifico, calcularsusceptibilidad; 

    densidadcalcular=false;
    calcularenergiamedia=false;  
    calcularmagnetizacion=false;
    simulacion=false; 
    calcularenergia=false; 
    calcularcalorespecifico=false; 
    calcularsusceptibilidad=true; 

    ofstream salida1, salida2, salida3, salida4, salida5, salida6 ;  
    if(simulacion==true) salida1.open("/home/sara/Escritorio/compu2324/Voluntarios/Voluntario_Kawasaki2/pruebapedro_3.dat"); 
    if(calcularmagnetizacion==true)     salida2.open("/home/sara/Escritorio/compu2324/Voluntarios/Kawasaki3/magnetizacion_N32_MC1000_it1000.dat"); ; 
    if(densidadcalcular==true)     salida3.open("/home/sara/Escritorio/compu2324/Voluntarios/Kawasaki3/densidad_N32_T35.dat");
    if(calcularenergiamedia==true) salida4.open("/home/sara/Escritorio/compu2324/Voluntarios/Kawasaki3/energiamedia_N32_T05.dat"); 
    if(calcularcalorespecifico==true) salida5.open("/home/sara/Escritorio/compu2324/Voluntarios/Kawasaki3/calorespecifico_N32_T05_3_ultimaversion.dat"); 
    if(calcularsusceptibilidad==true) salida6.open("/home/sara/Escritorio/compu2324/Voluntarios/Kawasaki3/susc_N32_T15.dat"); 



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

   double densidad[N]; 

   if(densidadcalcular==true) 
   {
    for(i=0;j<N;i++) 
    {
        if(reticula[i][5]==1) densidad[i]=1.0*reticula[i][5]; 
        else densidad[i]=0; 
    }
   }

    /*if(calcularmagnetizacion==true) 
    {   
        msup=msup/(N*N);
        minf=minf/(N*N); 
        salida2 << T << "  " << minf <<  "   " << msup << "  " << magnetizacion << endl; }*/ 


    //Paso Montecarlo: dentro de cada uno hago N*N Kawasakis (no segura)
    
    //for(T=0.5;T<=4.5; T=T+0.5)
    {
        aux1=0;
        aux2=0; 
        msup=0;
        minf=0; 
        magnetizacion=0; 
        energia=0.0; 
        energiacuadrado=0.0; 
        calor=0.0; 
        for(contadorMC=0; contadorMC<pasoMontecarlo; contadorMC++)
        {
            if(simulacion==true)
            {
            for(i=0; i<N; i++)
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
            }}
            salida1 << endl;
            aux1=0.;
            aux2=0; 
            for(iteraciones=0;iteraciones<paso;iteraciones++)
            {
                Intercambio(reticula,T);

            }

            if(calcularenergia==true)
            {
                for(i=0;i<N;i++)
                    for(j=0;j<N;j++)
                    {
                        energia+=(-0.5)*reticula[i][j]*(reticula[(i+1)%N][j]+reticula[(i-1+N)%N][j]+reticula[i][(j-1+N)%N]+reticula[i][(j+1+N)%N]); 
                    }

            }
            if(calcularcalorespecifico==true)
                {
                energiacuadrado+=energia*energia; 
                }
            if(densidadcalcular==true)
            {
                for(i=0;i<N;i++)
                {
                if(reticula[i][5]==1) densidad[i]+=1.0*reticula[i][5];  
                }
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
            if(calcularsusceptibilidad==true)
            {
              for(i=0; i<N;i++)
                {
                    for(j=0;j<N;j++)
                    {
                        magnetizacion+=reticula[i][j];
                        if(i<N/2) 
                        {aux1+=reticula[i][j];
                        minfcuadrado+=aux1*aux1; 
                        }
                        else 
                        {
                        aux2+=reticula[i][j]; 
                        msupcuadrado+=aux2*aux2; 
                        }
                    }
                }
                if(aux1<0) aux1=aux1*(-1);
                if(aux2<0) aux2=aux2*(-1); 
                aux1=0;
                aux2=0;

            }
        }

        if(calcularsusceptibilidad==true)
        {

            minf=aux1*1.0/(paso*pasoMontecarlo);
            msup=aux2*1.0/(pasoMontecarlo*paso);
            minfcuadrado=minfcuadrado/(pasoMontecarlo*paso); 
            msupcuadrado=msupcuadrado/(pasoMontecarlo*paso); 
            chisup=(msupcuadrado-msup)/(paso*T*T);
            chiinf=(minfcuadrado-minf)/(paso*T*T); 

            salida6 << T << chiinf << "   " << chisup << endl; 

        }
        if(densidadcalcular==true)
        {
            for(i=0;i<N;i++)
            {
                densidad[i]=densidad[i]*1.0/(pasoMontecarlo); //sumamos los espines igual que hicimos con la magnetización 
                salida3 << i << "  " << densidad[i] << endl;
            }

        }

        if(calcularmagnetizacion==true)
        {
            minf=aux1*1.0/(pasoMontecarlo*N*N);
            msup=aux2*1.0/(pasoMontecarlo*N*N); 
            magnetizacion=magnetizacion*1.0/(pasoMontecarlo*N*N); 

            salida2 << T << "  " << minf << "   " << msup << "  " << magnetizacion << endl;
        }
        if(calcularenergiamedia==true)
        {
            energia=energia/(pasoMontecarlo*2*N); 
            salida4 << T << "    " <<    energia << endl; 
        }

        if(calcularcalorespecifico==true )
        {
            energiacuadrado=energiacuadrado*1.0/(pasoMontecarlo);
            energia=energia/pasoMontecarlo; 
            calor=1.0/(T*T*paso)*((energiacuadrado)-(energia*energia)); 
            salida5 << T << "    " << calor << endl; 

        }
    }

    if(simulacion==true) salida1.close(); 
    if(calcularmagnetizacion=true) salida2.close(); 
    if(densidadcalcular==true) salida3.close(); 
    if(calcularenergiamedia==true) salida4.close(); 
    if(calcularsusceptibilidad==true) salida6.close(); 

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
                uniform_int_distribution<int> i_distribution(2, N-3); //initialize the distribution r_distribution

    double epp = 0.0, p=0.0 , xi=0.0, energia=0.0;

    double suma1=0.0, suma2=0.0; 
    int direction, i1,j1,i2,j2, espin1, espin2;
    
    xi=r_distribution(generator);
    direction=s_distribution(generator4);
    j1=j_distribution(generator2); 
    i1=i_distribution(generator3); 
    i2=i1;
    j2=j1; 

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
