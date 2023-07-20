/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"
#include "../../nifty_lib/lib.h"

using namespace std;


int main (int argc, char *argv[]){

   ofstream write;
   ofstream write2;
   write.open("Random_Walks_d.txt");
   write2.open("Random_Walks_c.txt");
   
   write << "RW,err" <<endl;
   write2 << "RW,err" <<endl;
    
   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

   double a = 1.0;//lunghezza passo
   double x = 0.0;
   double y = 0.0;
   double z = 0.0;
   double mean_prog[100] = {};
   double mean_prog2[100] = {};

   //RW nel discreto

   int M = 10000;
   int t_lim = 100;
   int B = 1;
   int F = M/B;
   double avg[100]={};

for (int k=0; k<B;k++){
  
  for(int j = 0; j<F; j++ ){

    x = 0.0;
    y = 0.0;
    z = 0.0;
   
    double r_x = 0.0;
    double r_y = 0.0;
    double r_z = 0.0;
      
   loading(M,j);

     for(int i = 0; i<t_lim;i++){
         
	double indicator = rnd.Rannyu(0.0,3.0);
    
    if(indicator <1.0){
	r_x = ((int)rnd.Rannyu(0.,2.)*2-1);// questo permette di fare +-1 e non zero!
    x += int(r_x*a);
    }
         
	if(indicator>=1.0 and indicator <2.0){
    r_y = ((int)rnd.Rannyu(0.,2.)*2-1);
    y += int(r_y*a);
    }
         
    if(indicator >=2.0){
	r_z = ((int)rnd.Rannyu(0.,2.)*2-1);
    z += int(r_z*a);
    }

         
	avg[i] += (x*x+y*y+z*z)/double(F); 
    	

     }
    
    }
  }
 
   for(int i = 0; i < t_lim;i++){
    mean_prog[i] += avg[i]/double(B);
    mean_prog2[i] += avg[i]*avg[i];
                    // propagazione errore: normalizzo per 1/2\sqrt{mean}

    write << sqrt(mean_prog[i]) << "," << error(mean_prog[i]/(i+1),mean_prog2[i]/(i+1),i)/(2.*sqrt(mean_prog[i])) << endl;
 
   }
   
   //RW nel continuo
   double mean_prog_bis[100] = {};
   double mean_prog2_bis[100] = {};
   double avg2[100] = {};
  
  
  for (int k=0; k<B;k++){
     loading(B,k);

  for(int j = 0; j<F; j++ ){

    x = 0.0;
    y = 0.0;
    z = 0.0;
   

     for(int i = 0; i<t_lim;i++){
		
	double th = rnd.Rannyu(0.0,M_PI);
	double phi = rnd.Rannyu(0.0,2.0*M_PI);

	x += a*sin(th)*cos(phi);
	y += a*sin(th)*sin(phi);
	z += a*cos(th);
	
	avg2[i] += (x*x+y*y+z*z)/double(F); 
    	

     }
    
    }
  }
 
   for(int i = 0; i<t_lim;i++){
    mean_prog_bis[i] += avg2[i]/double(B);
    mean_prog2_bis[i] += avg2[i]*avg[i];
    
    write2 << sqrt(mean_prog_bis[i]) << "," << error(mean_prog_bis[i]/(i+1),mean_prog2_bis[i]/(i+1),i)/(2.*sqrt(mean_prog[i])) << endl;
 
   }

  write.close();
  write2.close();

   return 0;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
