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


/*double error(double AV, double AV2, int n) 
{
   if(n==0){
      return 0;
   }
   else{
      return sqrt((AV2-AV*AV)/n);
   }
}*/

int main (int argc, char *argv[]){

   ofstream write;
   write.open("test.txt");
   
   write << "mean" << ","<<"error_mean"<<"," << "std"<<","<< "error_std" << endl;

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

   int N = 1000000; //numero di lanci
   double mean = 0;
   double varsig = 0;
   int B=100;  // numero di blocchi

   int E = int(N/B); //numero elementi per blocchi
   
   double mean_progr = 0; // variabili di appoggio per la media
   double mean_progr2 = 0; // variabili di appoggio per la media
   double sig_progr = 0; // variabili di appoggio per la std (sigma)
   double sig_progr2 = 0; // variabili di appoggio per la std (sigma)

   for(int i=0; i<B; i++){
    // devo iterare sui blocchi

     for (int j = 0; j < E; j++){

      double r = rnd.Rannyu();
      mean += r ;
      varsig += pow((r-0.5),2);
      }
   
      mean = mean/E;
      mean_progr += mean; 
      mean_progr2 += mean*mean; 

      varsig = varsig/E;
      sig_progr += varsig;
      sig_progr2 += varsig*varsig;

      if(write.is_open()){
	write << mean_progr/(i+1) <<"," << error(mean_progr/(i+1),mean_progr2/(i+1),i) <<"," <<sig_progr/(i+1) <<"," << error(sig_progr/(i+1),sig_progr2/(i+1),i) << endl;
   }else{cerr << " ERROR: file not opened" <<endl;}

   }
   write.close();
  
   ofstream chi_test;
   chi_test.open("chi.txt");
   chi_test << "chi^2"  <<endl;

   Random rnd2;
   rnd2.SaveSeed();
   
   int M = 100; // Numero di intervalli
   double n = 1e4; // Numeri di lanci per intervallo
   
   for(int i = 0; i<M; i++){// ciclo bocchi

     double chi = 0;
    
     double r = rnd2.Rannyu(0,1);
  
     for(int j =0; j<n; j++){//ciclo dei lanci in un blocco
     
     }
      
      chi = pow((r-n/M),2)/(n/M);
      // Salvo il risultato 
      if (chi_test.is_open()){
         chi_test << chi << endl; //
      } else cerr << "PROBLEM: Unable to open chi.out" << endl;

   // Torno all'inizio del ciclo e ripeto per l'esperimento successivo

   }

  




   chi_test.close();
   rnd.SaveSeed();
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
