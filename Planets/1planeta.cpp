using namespace std;
#include <math.h>
#include <stdio.h>
#include <iostream>
#define Msol (1.989*pow(10,30))
#define c (1.496*pow(10,11))
#define G (6.67*pow(10,-11))
#include <fstream>

int main ()
{
    double M[9],Rx[9],Ry[9],Ry_aux[9],Vx[9],Vy[9],ax[9],ay[9],Epot[9];
    double reet;
    double wx[9],wy[9];
    double h;
    double T[9],nvuelta[9],cinetica,vtotal,Etot;
    int i,j,k,contador;
    ofstream salida1, salida2, salida3; 

    salida1.open("/home/sara/Escritorio/compu2324/Planetas/datos_1.dat"); 
    salida2.open("/home/sara/Escritorio/compu2324/Planetas/energia_1.dat");
    salida3.open("/home/sara/Escritorio/compu2324/Planetas/periodos_1.dat");


    //reescalado del tiempo
    reet=sqrt((c/(G*Msol)));
    //Inicializo datos del Sol
    M[0]=1;
    Rx[0]=0;
    Ry[0]=0;
    contador=1;

    //Inicializo todos los datos que son 0
    for(i=0;i<9;i++)
    {
        Ry[i]=0;
        Vx[i]=0;
        ax[i]=0;
        ay[i]=0;
        wx[i]=0;
        wy[i]=0;
        nvuelta[i]=0;
    }
    //Datos de Mercurio
    M[1]=(0.33*pow(10,24))/Msol;
    Rx[1]=(57.6*pow(10,9)/c);
    Vy[1]=(47900*reet);

    //Datos de Venus
    M[2]=(4.87*pow(10,24))/Msol;
    Rx[2]=(108.2*pow(10,9))/c;
    Vy[2]=(35000*reet);

    //Inicializo datos de la Tierra
    M[3]=(5.97*pow(10,24))/Msol;
    Rx[3]=1;
    Vy[3]=29800*reet;

    //Datos de Marte
    M[4]=(0.642*pow(10,24))/Msol;
    Rx[4]=(227.9*pow(10,9))/c;
    Vy[4]=24100*reet;

    //Datos Júpiter
    M[5]=1899*pow(10,24)/Msol;
    Rx[5]=(778.6*pow(10,9))/c;
    Vy[5]=(13100*reet);

    //Datos Saturno
    M[6]=(568*pow(10,24))/Msol;
    Rx[6]=(1433.5*pow(10,9))/c;
    Vy[6]=(9700*reet);

    //Datos de Urano
    M[7]=(86.8*pow(10,24))/Msol;
    Rx[7]=(2872.5*pow(10,9))/c;
    Vy[7]=(6800*reet);

    //Datos de Neptuno
    M[8]=(102*pow(10,24))/Msol;
    Rx[8]=(4495.1*pow(10,9))/c;
    Vy[8]=(5400*reet);

    //Pongo una h inicial
    h=0.01;
    //aceleración inicial de los planetas
    for(i=1;i<9;i++)
    {
        for(j=0;j<9;j++)
        {
            ax[i]=ax[i]-M[j]*(Rx[i]-Rx[j])/pow((pow((Rx[i]),2)+pow((Ry[i]),2)),1.5);
            ay[i]=ay[i]-M[j]*(Ry[i]-Ry[j])/pow((pow((Rx[i]),2)+pow((Ry[i]),2)),1.5);
        }
    }

    //Ahora el algoritmo
    for(double t=0; t<5000;t=t+h)
    
    {
        //Primero evalúo las posiciones
        for(i=1;i<9;i++)
        {
            //primero guardo la y anterior para calcular el periodo
            Ry_aux[i]=Ry[i];

            //Ahora calculo las nuevas posiciones
            Rx[i]=Rx[i]+h*Vx[i]+(pow(h,2)/2)*ax[i];
            Ry[i]=Ry[i]+h*Vy[i]+(pow(h,2)/2)*ay[i];
        }

        //Ahora evaúo el vector auxiliar
        for(i=1;i<9;i++)
        {
            wx[i]=Vx[i]+(h/2)*ax[i];
            wy[i]=Vy[i]+(h/2)*ay[i];
        }
       
        //Vuelvo a inicializar las aceleraciones
        for(int i=1;i<9;i++)
        {
            ax[i]=0;
            ay[i]=0;
            //inicializo la energia potencial
            Epot[i]=0;
        }
        //inicializo las variableas de energia totales
        cinetica=0;
        vtotal=0;
        //Vuelvo a evalúar las aceleraciones con la nueva posición
        for(i=1;i<9;i++)
        {
            for(j=0;j<9;j++)
            {
                if(j!=i)
                {
                ax[i]=ax[i]-M[j]*(Rx[i]-Rx[j])/pow((pow((Rx[i]-Rx[j]),2)+pow((Ry[i]-Ry[j]),2)),1.5);
                ay[i]=ay[i]-M[j]*(Ry[i]-Ry[j])/pow((pow((Rx[i]-Rx[j]),2)+pow((Ry[i]-Ry[j]),2)),1.5);
                }
            }
        }
       
        //Evalúo las velocidades nuevas
        for(i=1;i<9;i++)
        {
            Vx[i]=wx[i]+(h/2)*ax[i];
            Vy[i]=wy[i]+(h/2)*ay[i];
        }
        //Voy a sacar los valores por a un archivo
        for(i=0;i<9;i++)
        {
            salida1 << Rx[i]*1.0 << ", " << Ry[i]*1.0 << endl; 
        }
        salida1 << endl; 

        //Calculo el período de cada planeta
        for(k=contador;k<9;k++)
        {
            if(Rx[k]>0)
            {
                if((Ry_aux[k]*Ry[k])<0)
                {
                    nvuelta[k]=nvuelta[k]+1;
                    T[k]=(t/nvuelta[k])*pow((pow(c,3)/(G*Msol)),0.5)/(3600*24);
                    contador=contador+1;
                    salida3 << T[k] << endl; 
                }
            }
        }

        //Calculo la energía del sistema
        //Primero la cinética
        for(i=1;i<9;i++)
        {
            cinetica=cinetica+(0.5)*M[i]*(pow(Vx[i],2)+pow(Vy[i],2));
        }
        //Ahora la potencial
        for(i=1;i<9;i++)        
        {
            for(j=0;j<9;j++)
            {
                if(j!=i)
                {
                    Epot[i]=Epot[i]+((G*M[i]*M[j])/pow((pow((Rx[i]-Rx[j]),2)+pow((Ry[i]-Ry[j]),2)),0.5));
                }
                vtotal=vtotal+Epot[i];
            }
        }
        //las sumo
        Etot=cinetica+vtotal;
        //las saco por f2
        salida2 << Etot << endl; 
    }
    salida1.close();
    salida2.close();
    salida3.close(); 

    return 0;
}