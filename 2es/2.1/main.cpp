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

double sampl(double r){
	
	return 1+sqrt(1-r);
	
	}

int main (int argc, char *argv[]){

   ofstream write;
   ofstream write2;
   write.open("Integral.txt");
   write2.open("Integral_2.txt");
   
   write << "I_1,err_1" << endl;
   write2 << "I_2,err_2" <<endl;
   
    
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


 //Prima parte si fa il sampling con la fuznione uniforme 
 int L = 100000000; //numero di lanci
 int B = 100;//numero di blocchi
 int F = int(L/B);
 double I = 0.0;
 double sum = 0.0;
 double mean_prog = 0.0 ;
 double mean_prog2 = 0.0;
 double err,err_2 = 0.0;

cout << " sampling uniforme" << endl;
 for (int i = 0;i<B; i++){
   loading(B,i); 
     //sum =0.0;
     
    for (int j = 0; j<F;j++){
 
 	double r = rnd.Rannyu(0.0,1.0);

 	sum += M_PI/2 * cos(r*M_PI/2);
	}
     
     sum = sum/F;
     I = sum;
     mean_prog += sum;
     mean_prog2 += sum*sum;
     
     //cout <<I<<","<<error(mean_prog2/(i+1),mean_prog/(i+1),i) << ",";  
     write << mean_prog/(i+1) << "," << error(mean_prog/(i+1),mean_prog2/(i+1),i) << endl;

    }
     
     sum =0.0;
     mean_prog = 0;
     mean_prog2 = 0;
     double prog_svi =0.0;
     double prog_svi2 =0.0;
     double mean_svi =0.0;

 cout << " sampling non uniforme" << endl;

 for (int i = 0;i<B; i++){
   loading(B,i);
     //Ora per la seconda parte si deve utillizzare l'integrale del sampling e si prende lo sviluppo al quart'ordine del coseno
     
     for (int i = 0;i<F; i++){

     double r = rnd.Rannyu(0.0,1.0);
     double u = rnd.Rannyu(0.0,M_PI/2.0);


     sum += M_PI/2 * cos(sampl(r)*M_PI/2)/(2.0*(1-sampl(r)));

 } 
  
     sum = sum/F;
     prog_svi += sum;
     prog_svi2 += sum*sum;
  //cout << I* sum<< ","<< I*sum*sqrt(pow(err/I,2)+pow(err_2/sum,2)) << endl;//devo propagare l'errore poichè moltiplico i due risultati 
  write2 << prog_svi/(i+1)<< ","<< error(prog_svi/(i+1),prog_svi2/(i+1),i) << endl;//devo propagare l'errore poichè moltiplico i due risultati 
 
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
