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
	return exp(-pow(x-mu,2.0)/(2.0*pow(sigma,2.0))) + exp(-pow(x+mu,2.0)/(2.0*pow(sigma,2.0)));
 }

double p(double x, double mu, double sigma){
	return psiT(x,mu,sigma)*psiT(x,mu,sigma);
 }

double second_derivative(double x, double mu, double sigma){
	return ((x-mu)*(x-mu)/pow(sigma,4.0)-1.0/pow(sigma,2.0))*exp(-(x-mu)*(x-mu)/(2.0*sigma*sigma)) + ((x+mu)*(x+mu)/pow(sigma,4.0)-1.0/pow(sigma,2.0))*exp(-(x+mu)*(x+mu)/(2.0*sigma*sigma));

}

double potential(double x){
	return pow(x,4.0) - 2.5*pow(x,2.0);
}

double fun_to_integrate(double x, double mu, double sigma){
	
	double fun =  (-(0.5)*second_derivative(x,mu,sigma)/psiT(x,mu,sigma) + potential(x));
	
	return fun;
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

double mean_progr_H = 0.0;
double mean_progr2_H = 0.0;

int count = 0;

for (int i = 0; i<B; i++){

    H = 0.0;
	
    count = 0;

    for (int j = 0; j<F; j++){

       double largh = 2.10;

/*===============Metropolis step=======================*/
       Tx = x + rnd.Rannyu(-largh,largh);
       
       double r = rnd.Rannyu(0.0,1.0);
       double T = p(Tx,mu,sigma)/p(x,mu,sigma);

       double A = min(1.0,T);

       if( r <= A ){
	
		x = Tx;
 		count ++;

        }
        
/*======================  Eval H ====================*/
    H += fun_to_integrate(x,mu,sigma); 
        
    }
    
    H /= (double(F));
    mean_progr_H += H;
    mean_progr2_H += H*H;
    
    //mean_progr_H = mean_progr_H/double(i+1);

   }    
  
return mean_progr_H/B;
    
}

/*=======================================MAIN===============================================*/

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


/*=====================Metro per H e distribuzione=======================================*/

ofstream write;
ofstream write_H;
ofstream points;

write.open("modulus.txt");
write_H.open("H.txt");
points.open("points2.txt");

write << "I,err" << endl;
write_H << "H,err" << endl;
points << "x" << endl;

int B = 100;
int M = pow(10,7);
int F = M/B;

double x = 0.0; // x in (-3,3)
double Tx = 0.0;

double H = 0.0;
double I = 0.0;
double mu = 0.8;
double sigma = 0.6;

double mean_progr = 0.0;
double mean_progr2 = 0.0;
double mean_progr_H = 0.0;
double mean_progr2_H = 0.0;

int count = 0;

for (int i = 0; i<B; i++){

    I = 0.0;
    H = 0.0;
	
    count = 0;

    for (int j = 0; j<F; j++){

       double largh = 2.625;
       Tx = x + rnd.Rannyu(-largh,largh);
       
       double r = rnd.Rannyu(0.0,1.0);
       double T = p(Tx,mu,sigma)/p(x,mu,sigma);

       double A = min(1.0,T);

       if(r <= A){

            x = Tx;
            count ++;
        }
        if(j%(F/1000))
                    points << x << endl;

	I += p(x,mu,sigma);
	H += fun_to_integrate(x,mu,sigma); 

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
points.close();
    
    
/*=======================Simulated annealing===============================*/

//Ora devo usare la "Temperatura" per trovare il minimo di H

count = 0;
int beta_steps = 1000;


//Inizializzazioni

double mean_progr_mu = 0.0;
double mean_progr2_mu = 0.0;
double mean_progr_sigma = 0.0;
double mean_progr2_sigma = 0.0;

mean_progr_H = 0.0;
mean_progr2_H = 0.0;

//valori iniziali
double beta = 1.0;
double incr_beta = 1.5;

mu = 1.0;    
sigma = 0.5;    
double Lmu = 1.0;
double Lsigma = 1.0;
    
    
ofstream write_mu;
write_mu.open("valori.txt");

write_mu << "beta,mu,err_mu,sigma,err_sigma,H,err_H" << endl;

double vecchiaH = valore_H(mu,sigma);

for(int i = 0; i < beta_steps; i++){
    // devo definire una probabilità con la temperatura (cioè la beta) e usare metropolis
    
    //devo definire dei nuovi parametri
    loading(beta_steps,i);
    
    
    double d_mu = rnd.Rannyu(-Lmu,Lmu)/beta;
    double d_sigma = rnd.Rannyu(-Lsigma,Lsigma)/beta;

    mu += d_mu;
    sigma += d_sigma;
    
    double nuovaH = valore_H(mu,sigma);
    
    double deltaH = vecchiaH - nuovaH; //deve essere la differenza tra le hamiltoniane
  
    double q = exp(beta*deltaH);
    double A = min(1.0,q);
    
    double r = rnd.Rannyu();
    
    if(r <= A)
      {
        
        vecchiaH = nuovaH;
        count++;           
        
	H = vecchiaH;
	mean_progr_H += H;
	mean_progr2_H += H*H;
	
	mean_progr_mu += mu;
        mean_progr2_mu += mu*mu;

	mean_progr_sigma += sigma;
	mean_progr2_sigma += sigma*sigma;
        
        
        write_mu << beta << "," << mu << "," << error(mean_progr_mu/(beta+1),mean_progr2_mu/(beta+1),beta) << "," << sigma << "," << error(mean_progr_sigma/(beta+1),mean_progr2_sigma/(beta+1),beta) << "," << vecchiaH << "," << error(mean_progr_H/(beta+1),mean_progr2_H/(beta+1),beta) << endl;
        
        beta += incr_beta;
      }
    
    else
      {
        mu -= d_mu; 
        sigma -= d_sigma;
        i--;
      }
    
  
    
    }
    
write_mu.close();

/*==============================SAMPLING FINAL========================================*/


ofstream write_H_final;

write_H_final.open("H_final.txt");

write_H_final << "H,err" << endl;

B = 100;
M = pow(10,7);
F = M/B;

x = 0.0; // x in (-3,3)
Tx = 0.0;

H = 0.0;
I = 0.0;
mu = 0.8;
sigma = 0.6;

mean_progr = 0.0;
mean_progr2 = 0.0;
mean_progr_H = 0.0;
mean_progr2_H = 0.0;

count = 0;

for (int i = 0; i<B; i++){

    I = 0.0;
    H = 0.0;
	
    count = 0;

    for (int j = 0; j<F; j++){

       double largh = 2.625;
       Tx = x + rnd.Rannyu(-largh,largh);
       
       double r = rnd.Rannyu(0.0,1.0);
       double T = p(Tx,mu,sigma)/p(x,mu,sigma);

       double A = min(1.0,T);

       if(r <= A){

            x = Tx;
            count ++;
        }

	I += p(x,mu,sigma);
	H += fun_to_integrate(x,mu,sigma); 

    }

    I /= double(F);
    mean_progr += I;
    mean_progr2 += I*I;

    H /= double(F);
    mean_progr_H += H;
    mean_progr2_H += H*H;


    write_H_final << mean_progr_H/(i+1) << "," << error(mean_progr_H/(i+1),mean_progr2_H/(i+1),i) << endl;

   cout << " acc_rate batch "<< i+1 <<": " << ((double)count/F)*100. << "%" << endl;
     
}

write_H_final.close();

    
    
    
return 0;
}


