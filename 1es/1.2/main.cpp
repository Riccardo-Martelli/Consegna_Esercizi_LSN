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

double inv_exp(double p, double lam = 1){

	return -log(1-p)/lam;
}

double inv_cauchy(double p ,double mu=0, double gamma=1){
     
	return gamma*tan(M_PI*(p-0.5))+mu;
    
}

int main (int argc, char *argv[]){

   ofstream write;
   write.open("test.txt");
   
   write << "unif,exp,cauchy" <<endl;
    
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

  int N = 100;
  //double arr[N] ={};
 
  for (int i=0; i<N; i++){

  double r = rnd.Rannyu(0,1);
  
    write << r <<","<< inv_exp(r) <<","<< inv_cauchy(r) <<endl;
     
  } 
 write.close();
 
 N=100;
 ofstream write2;
 string str = "real"+to_string(N)+".txt";
 write2.open(str);
 
 write2 << "Sum_uniform,Sum_exponential,Sum_cauchy" <<endl;
    
 int realizations = pow(10,4);
    
 for (int i = 0; i <realizations; i++){
  
     double sum_unif = 0;
     double sum_exp = 0;
     double sum_cauchy = 0;
  for (int j=0; j<N ; j++)  { 
     double r = rnd.Rannyu(0.0,1.0);
     sum_unif += r/N;
     sum_exp += inv_exp(r)/N;
     sum_cauchy += inv_cauchy(r)/N;
  }
    write2 << sum_unif<< "," << sum_exp << "," << sum_cauchy << endl;
     
 }
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
