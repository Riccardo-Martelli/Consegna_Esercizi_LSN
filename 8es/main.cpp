#include <cmath>
#include <iostream>
#include <random>
#include <fstream>
#include <string>
#include <algorithm>
#include "../Random/random.h"
#include "../nifty_lib/lib.h"

using namespace std;


double psiT(double x, double mu, double sigma){
	return (exp(-(x-mu)*(x-mu)/(2.0*sigma*sigma))+ exp(-(x+mu)*(x+mu)/(2.0*sigma*sigma)))/(2.0*sqrt(2.0*M_PI)*sigma);
 }

double p(double x, double mu, double sigma){
	return pow(psiT(x,mu,sigma),2.0);

 }

double second_derivative(double x, double mu, double sigma){
	return (((x-mu)*(x-mu)/pow(sigma,4)-1.0/pow(sigma,2))*exp(-(x-mu)*(x-mu)/(2.0*sigma*sigma)) + ((x+mu)*(x+mu)/pow(sigma,4)-1.0/pow(sigma,2))*exp(-(x+mu)*(x+mu)/(2.0*sigma*sigma)))/(2.0*sqrt(2.0*M_PI)*sigma);

}

double potential(double x){
	return pow(x,4) - 2.5*pow(x,2);
}

double fun_to_integrate(double x, double mu, double sigma, double I){
	
	double fun =  p(x,mu,sigma)*(-(0.5)*second_derivative(x,mu,sigma)/psiT(x,mu,sigma) +potential(x));
	
	return fun/I;
}

double valore_H(double mu, double sigma){
    
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

    
int B = 100;
int M = pow(10,4);
int F = M/B;

double x = 0.0; // x in (-3,3)
double Tx = 0.0;

double H = 0.0;
double I = 0.0;

double mean_progr = 0.0;
double mean_progr2 = 0.0;
double mean_progr_H = 0.0;
double mean_progr2_H = 0.0;

int count = 0;

for (int i = 0; i<B; i++){

    I = 0.0;
	
    count = 0;

    for (int j = 0; j<F; j++){

       double largh = 2.35;
        
       Tx = x + rnd.Rannyu(-largh,largh);
       
       double r = rnd.Rannyu();
       double T = p(Tx,mu,sigma)/p(x,mu,sigma);

       double A = min(1.0,T);

       if( r <= A ){
	
		x = Tx;
 		count ++;

        }
	
    I += p(x,mu,sigma);
	H += fun_to_integrate(x,mu,sigma,I); 
        
    }

    I /= double(F);
    mean_progr += I;
    mean_progr2 += I*I;

    H /= double(F);
    mean_progr_H += H;
    mean_progr2_H += H*H;

    mean_progr_H = mean_progr_H/double(i+1);

   }    
  
return mean_progr_H;
    
}


int main(int argc, char *argv[]){

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


ofstream write;
ofstream write_H;

write.open("modulus.txt");
write_H.open("H.txt");

write << "I,err" << endl;
write_H << "H,err" << endl;

int B = 100;
int M = pow(10,6);
int F = M/B;

double x = 0.0; // x in (-3,3)
double Tx = 0.0;

double H = 0.0;
double I = 0;
double mu = 0.5;
double sigma = 1.0;

double mean_progr = 0.0;
double mean_progr2 = 0.0;
double mean_progr_H = 0.0;
double mean_progr2_H = 0.0;

int count = 0;

for (int i = 0; i<B; i++){

    I = 0.0;
	
    count = 0;

    for (int j = 0; j<F; j++){

       double largh = 2.35;
       Tx = x + rnd.Rannyu(-largh,largh);
       
       double r = rnd.Rannyu();
       double T = p(Tx,mu,sigma)/p(x,mu,sigma);

       double A = min(1.0,T);

       if(r<=A){

            x = Tx;
            count ++;

        }
        
	I += p(x,mu,sigma);
	H += fun_to_integrate(x,mu,sigma,I); 

    }

    I /= double(F);
    mean_progr += I;
    mean_progr2 += I*I;

    H /= double(F);
    mean_progr_H += H;
    mean_progr2_H += H*H;


    write << mean_progr/(i+1) << "," << error(mean_progr/(i+1),mean_progr2/(i+1),i) << endl;
    write_H << mean_progr_H/(i+1) << "," << error(mean_progr_H/(i+1),mean_progr2_H/(i+1),i) << endl;

   cout << " acc_rate batch "<< i+1 <<": " << ((double)count/F)*100. << "%" << endl;
     
}

write.close();
write_H.close();
    
    
//Ora devo usare la T per trovare il minimo di H

count = 0;
int beta_steps = 5000;
    
double Lmu = 1.0;
double Lsigma = 1.0;
    
double beta = 1.0;
double incr_beta =1.5;
    
ofstream write_mu;
write_mu.open("valori.txt");

write_mu << "beta,mu,sigma,H" << endl;

for(int i = 0; i < beta_steps; i++){
    // devo definire una probabilità con la temperatura(beta) e usare metropolis
    
    //devo definire dei nuovi parametri
    loading(beta_steps,i);
    
    double vecchiaH = valore_H(mu,sigma);
    
    double d_mu = rnd.Rannyu(-Lmu,Lmu)/beta;
    double d_sigma = rnd.Rannyu(-Lsigma,Lsigma)/beta;

    mu += d_mu;
    sigma += d_sigma;
    
    double nuovaH = valore_H(mu,sigma);
    
    double deltaH = vecchiaH - nuovaH; //deve essere la differenza tra le hamiltoniane
    
    double p = exp(beta*deltaH);
    double A = min(1.0,p);
    
    double r = rnd.Rannyu();

    
    if(r <= A)
      {
        
        vecchiaH = nuovaH;
        count++;           
        
      }
    
    else
      {
        mu -= d_mu; 
        sigma -= d_sigma;
      }
    
    //cout  << vecchiaH << endl;
    write_mu << beta << "," << mu << "," << sigma << "," << vecchiaH << endl;
    
    beta += incr_beta;
    
    }
    
write_mu.close();


return 0;
}


