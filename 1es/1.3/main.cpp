/****************************************************************
*****************************************************************
    _/_/_/_/  _/     _/       Numerical Simulation Laboratory
   _/    _/  _/_/ _/_/      Physics Department
  _/_/_/_/  _/  _/ _/     Universita' degli Studi di Milano
 _/  _/    _/     _/    Riccardo Martelli
_/    _/  _/     _/   email:riccardo.martelli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"
#include "../../nifty_lib/lib.h"

using namespace std;

double s_distr(double r,double x,double L){
	
	double sq = 1.5*sqrt(L*L-r*r)/L;
	double test = x;

	do{
	sq = sqrt(L*L-r*r);
	test = x;
	}while(test > sq);
	  return sq;
}

int main (int argc, char *argv[]){

   ofstream write;
   write.open("Buffon.txt");
   
   write << "PI" <<endl;
   
   ofstream write2;
   write2.open("err.txt");
   
   write2 << "err_pi" <<endl;
    
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

 int B = 100; //Number of blocks
 int M = 1000; //Number of throws
 int F = M/B; // throws per block


 double Pi = 0.0;
 int count = 0;
 
 double r = 0.0;
 double x = 0.0;
 double y = 0.0;
 double pos = 0.0;
 double mean_prog = 0.0;
 double mean_prog2 = 0.0;

 double L = 3.0 ; //lunghezza ago
 double d = 5.0 ; //distanza linee

 //Dalla teoria dipende da dove cade l'ago e chiaramente da con che angolo cade, oltre che dalla sua lunghezza
    
 for (int i=0; i<B; i++){
   
   count = 0;
   
   for (int j=0; j<F; j++){
 	
    r = rnd.Rannyu(0.0,d); //Posizione ago in veritcale
    x = rnd.Rannyu(-L,L); // L*cos(theta)
    y = rnd.Rannyu(-L,L); // L*cos(theta)
    //y = sqrt(L*L-x*x); // uso il fatto che L^2 = x^2 + y^2
    
    if(y*y+x*x <= L*L) 
     pos = r+L*sin(atan(y/x));
    
    else j--;
    //pos = r + s_distr(x,u,L);   
    
    if (pos >= d or pos <= 0.0)
        count ++;//Number of hits
      
   }
    
    if(count != 0){
    
    Pi = 2.0*L*double(F)/double(count);
    Pi /= d;
    
    //mean += Pi;
    mean_prog += Pi;
    mean_prog2 += Pi*Pi;

    write << mean_prog/(i+1) << endl;
    write2 << error(mean_prog/(i+1),mean_prog2/(i+1),i) << endl;//mettere gli errori vari
    
    }
    else
      i--;

   loading(B,i);
	
 }
 cout << endl;

 write.close();
 write2.close();

   return 0;
}

/****************************************************************
*****************************************************************
    _/_/_/_/  _/     _/       Numerical Simulation Laboratory
   _/    _/  _/_/ _/_/      Physics Department
  _/_/_/_/  _/  _/ _/     Universita' degli Studi di Milano
 _/  _/    _/     _/    Riccardo Martelli
_/    _/  _/     _/   email:riccardo.martelli@unimi.it
*****************************************************************
*****************************************************************/
