#include <cmath>
#include <iostream>
#include <random>
#include <fstream>
#include <string>
#include "../Random/random.h"
#include "../nifty_lib/lib.h"

using namespace std;

double N(double x)
{
 return 0.5*(1+erf(x/sqrt(2.))) ; 
}

int main(int argc, char *argv[]){
// Parameters
double r = 0.1;
int steps = 100;
double delta_t = 1.0/(double)steps;
float t = 0.0;
float T = 1.0;
double sigma = 0.25;
float K = 100.0;
int S0 = 100.0;
     
ofstream write;
ofstream write2;

   
write.open("res.txt");
write << "C_one,P_one,C_discr,P_discr" << endl;
write2.open("err.txt");
write2 << "err_C_one,err_P_one,err_C_discr,err_P_discr" << endl;

Random rnd;

int seed[4];
int p1,p2;

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


int M = 100000;
int B = 100;//Blocchi
int F = M/B;//Trials per blocco
double S_F = 0.0;
double C_f = 0.0;
double P_f = 0.0;
double C_d = 0.0;
double P_d = 0.0;
double sum [4] = {};
double mean_prog[4] = {};
double mean_prog2[4] = {};  
    
    
for(int i = 0; i<B; i++){
    
  for(int j=0; j<F; j++){
      
    // Direttamente step finale
    double Z = rnd.Gauss(0.0,1.0);
    double S_T = S0*exp((r-sigma*sigma/2.) + sigma*Z);
      
    C_f = exp(-r*T)*max(0., S_T-K);
    P_f = exp(-r*T)*max(0.,-S_T+K);
      
    sum[0] += C_f;
    sum[1] += P_f;
          
    // Discretize approach
    double S[steps];
    S[0] = S0;

    for(int i=1; i <=steps; i++){

     double Z_i = rnd.Gauss(0.0,1.0);
     S[i] = S[i-1]*exp((r-sigma*sigma/2.0)*delta_t + sigma*Z_i*sqrt(delta_t));

    }
      
   S_F = S[steps];
  
   C_d = exp(-r*T)*max(0., S_F-K);
   P_d = exp(-r*T)*max(0., K-S_F);
   sum[2]+= C_d;
   sum[3]+= P_d;

  } 
  
  for (int p = 0; p<=3;p++){
      
    sum[p] /= F;
    mean_prog[p] += sum[p];
    mean_prog2[p]+= sum[p]*sum[p];
      
  }    

  write << mean_prog[0]/(i+1) << "," << mean_prog[1]/(i+1) << ","<< mean_prog[2]/(i+1) << "," << mean_prog[3]/(i+1) << endl;
    
   for (int p = 0; p < 4;p++){
       
    write2 << error(mean_prog[p]/(i+1),mean_prog2[p]/(i+1),i)  << ",";
    
    if (p==3){
        
     write2 << error(mean_prog[p]/(i+1),mean_prog2[p]/(i+1),i) << endl;
    
    }   
   }
      for (int p = 0; p<=3;p++){
      
    sum[p] = 0.0;
      
  }    
 }

write.close();
write2.close();

    
return 0;
}
