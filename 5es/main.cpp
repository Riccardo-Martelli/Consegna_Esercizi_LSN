#include <cmath>
#include <iostream>
#include <random>
#include <fstream>
#include <string>
#include <algorithm>
#include "../Random/random.h"
#include "../nifty_lib/lib.h"

using namespace std;

double a0 = 1.0; 
    //0.0529;

double psi1(double x,double y,double z){
    double rag = sqrt(x*x+y*y+z*z);
    
	return exp(-2.*rag)/M_PI;
}

double psi2(double x,double y,double z){
    double rag = sqrt(x*x+y*y+z*z);
	return pow(z,2.0)*exp(-rag)/(32.0*M_PI);
}

double minimo(double A,double B){

	if (A <= B)
	  return A;
	else
	  return B;

}

/*====================================MAIN======================================*/

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

ifstream ReadInput;

ReadInput.open("Input.in");

bool gaus;
ReadInput >> gaus;  
    
ReadInput.close();

string sum = "avg_psi_one_gaus.txt";
string point = "points_psi_one_gaus.txt";
string point2 = "points_psi_two_gaus.txt";
string sum2= "avg_psi_two_gaus.txt";
string equ = "gaus_equilibration.txt";

cout << "gaus=" << gaus <<endl;

if(gaus == false){

sum = "avg_psi_one_unif.txt";
point = "points_psi_one_unif.txt";

point2 ="points_psi_two_unif.txt";
sum2 = "avg_psi_two_unif.txt";

equ = "unif_equilibration.txt";

}


ofstream write2;
ofstream write3;
ofstream sums;
ofstream sums2;
ofstream far;
ofstream far2;
ofstream equil;

write2.open(point);    
write3.open(point2);

    
sums.open(sum);
sums2.open(sum2);
equil.open(equ);

    
sums << "sum,mean,err" << endl;
sums2 << "sum,mean,err" << endl;
equil << "sum,mean,err" << endl;
write2 << "x,y,z" << endl;
write3 << "x,y,z" << endl;
    
int steps= pow(10,5);
int B = 100;
int F = steps/B;


double sumr = 0.0;
double mean_prog = 0.0;
double mean_prog2 = 0.0; 
    
double x0 = 1.;
double y0 = 0.;
double z0 = 0.;

double x = x0;
double y = y0;
double z = z0;
   
double Tx = 2.0;
double Ty = 0.0;
double Tz = 0.0;

int count = 0;
    
int eq = 5.0*pow(10,4);//iterations for equilibration
int eq_blcks = 5.0*pow(10,4);
int rateo = eq/eq_blcks; 
double largh = 0.0;
double r ;

cout << ".------------------------.  " << endl;
cout << "| Equilibration of psi_1 | " << endl;
cout << "'------------------------'  " << endl;
cout << endl;

for (int j=0; j< eq_blcks;j++){
        
   sumr=0.0;
   count = 0;

   for(int i =0; i< rateo;i++){

	if(!gaus){
        largh =1.225;
        
        Tx =  x + rnd.Rannyu( -largh, largh);
        Ty =  y + rnd.Rannyu( -largh, largh);
        Tz =  z + rnd.Rannyu( -largh, largh);
        
        }
        else if(gaus){
        largh = 0.75;         	
	
        Tx =  x + rnd.Gauss(0.0,largh);
        Ty =  y + rnd.Gauss(0.0,largh);
        Tz =  z + rnd.Gauss(0.0,largh);
	
	}
	
	double T = psi1(Tx,Ty,Tz)/psi1(x,y,z);
	double A = minimo(1.0,T);
        double r = rnd.Rannyu(0.0,1.0);

        if( r < A ){

          x = Tx;
          y = Ty;
          z = Tz;

          count ++;
             
        }
        sumr += sqrt(x*x + y*y + z*z);
	
       } 
       
    loading(eq_blcks,j);

    sumr = sumr/double(rateo);
    mean_prog += sumr;
    mean_prog2 += sumr*sumr; 

    equil << sumr << "," << mean_prog/(j+1) << "," << error(mean_prog/(j+1),mean_prog2/(j+1),j)  << endl;
}


cout << endl;
cout << ".----------------------------. " << endl;
cout << "| End equilibration of psi_1 | " << endl;
cout << "'----------------------------' " << endl;
cout << endl;

//####################################psi_{1,0,0}##########################################################

sumr =0.0;
mean_prog = 0.0;
mean_prog2 = 0.0; 

for (int j=0; j< B;j++){
        
   count = 0;

   for(int i =0; i< F;i++){

	if(!gaus){
        largh =1.225;
        
        Tx =  x + rnd.Rannyu( -largh, largh);
        Ty =  y + rnd.Rannyu( -largh, largh);
        Tz =  z + rnd.Rannyu( -largh, largh);
        
        }
        else if(gaus){
        largh = 0.75;         	
	
        Tx =  x + rnd.Gauss(0.0,largh);
        Ty =  y + rnd.Gauss(0.0,largh);
        Tz =  z + rnd.Gauss(0.0,largh);
	
	}
	
	double T = psi1(Tx,Ty,Tz)/psi1(x,y,z);
	double A = minimo(1.0,T);
        double r = rnd.Rannyu(0.0,1.0);

        if( r < A ){

          x = Tx;
          y = Ty;
          z = Tz;
          write2 << x << "," << y << "," << z << endl;

          count ++;
             
        }
        sumr += sqrt(x*x + y*y + z*z);
	
       } 
       
    loading(B,j);

    sumr = sumr/double(F);
    mean_prog += sumr;
    mean_prog2 += sumr*sumr; 

    sums << sumr << "," << mean_prog/(j+1) << "," << error(mean_prog/(j+1),mean_prog2/(j+1),j)  << endl;
}



//####################################psi_{2,1,0}##########################################################

 sumr = 0.0;
 mean_prog= 0.0;
 mean_prog2 = 0.0;
    
 x0 = 5.;
 y0 = 0.;
 z0 = 0.;
    
 x = x0;
 y = y0;
 z = z0;
 
 for (int j= 0; j < B; j++){


    int count = 0;
     
    for(int i = 0; i< F;i++){

	if(!gaus){
        largh =2.85;
        
        Tx =  x + rnd.Rannyu( -largh, largh);
        Ty =  y + rnd.Rannyu( -largh, largh);
        Tz =  z + rnd.Rannyu( -largh, largh);
        
        }
        else if(gaus){
        largh = 1.8;         	
	
        Tx =  x + rnd.Gauss(0.0,largh);
        Ty =  y + rnd.Gauss(0.0,largh);
        Tz =  z + rnd.Gauss(0.0,largh);
	
	}

        double A = min(1., psi2(Tx,Ty,Tz)/psi2(x,y,z));
        double r = rnd.Rannyu(0.0,1.0);

        if( r < A ){

          x = Tx;
          y = Ty;
          z = Tz;

          count ++;
        
          write3 << x << "," << y << "," << z << endl;
          
              }
          sumr += sqrt(x*x+y*y+z*z);
        }
        
    
    sumr /= (double)(F);
    mean_prog += sumr;
    mean_prog2 += sumr*sumr; 
    
    cout << " acc_rate_2 batch "<< j+1 <<": " << ((double)count/F)*100. << "%" << endl;

  sums2 << sumr << "," << mean_prog/(j+1) << "," << error(mean_prog/(j+1),mean_prog2/(j+1),j) << endl;

    }


/*=========================== far from origin  =================================*/

far.open("points_psi_one_far.txt");    
far2.open("points_psi_two_far.txt");

far << "x,y,z" << endl;
far2 << "x,y,z" << endl;

x0 = 20.;
y0 = 0.;
z0 = 0.;

x = x0;
y = y0;
z = z0;
   
/*=========================== psi 1  =================================*/
    
for (int j=0; j< B;j++){
        
   count = 0;

   for(int i =0; i< F;i++){

	if(!gaus){
        largh =1.225;
        
        Tx =  x + rnd.Rannyu( -largh, largh);
        Ty =  y + rnd.Rannyu( -largh, largh);
        Tz =  z + rnd.Rannyu( -largh, largh);
        
        }
        else if(gaus){
        largh = 0.75;         	
	
        Tx =  x + rnd.Gauss(0.0,largh);
        Ty =  y + rnd.Gauss(0.0,largh);
        Tz =  z + rnd.Gauss(0.0,largh);
	
	}
	
	double T = psi1(Tx,Ty,Tz)/psi1(x,y,z);
	double A = minimo(1.0,T);
        double r = rnd.Rannyu(0.0,1.0);

        if( r < A ){

          x = Tx;
          y = Ty;
          z = Tz;
          far << x << "," << y << "," << z << endl;

          count ++;
             
        }
       } 
}

    
 x0 = 50.;
 y0 = 0.;
 z0 = 0.;
    
 x = x0;
 y = y0;
 z = z0;
 
 /*=========================== psi 2  =================================*/
 
 for (int j= 0; j < B; j++){


    int count = 0;
     
    for(int i = 0; i< F;i++){

	if(!gaus){
        largh =2.85;
        
        Tx =  x + rnd.Rannyu( -largh, largh);
        Ty =  y + rnd.Rannyu( -largh, largh);
        Tz =  z + rnd.Rannyu( -largh, largh);
        
        }
        else if(gaus){
        largh = 1.8;         	
	
        Tx =  x + rnd.Gauss(0.0,largh);
        Ty =  y + rnd.Gauss(0.0,largh);
        Tz =  z + rnd.Gauss(0.0,largh);
	
	}

        double A = min(1., psi2(Tx,Ty,Tz)/psi2(x,y,z));
        double r = rnd.Rannyu(0.0,1.0);

        if( r < A ){

          x = Tx;
          y = Ty;
          z = Tz;

          count ++;
        
          far2 << x << "," << y << "," << z << endl;
          
              }
        }
        
    

    }



write2.close();  
write3.close();    
sums.close();
sums2.close();
equil.close();



return 0;
}
