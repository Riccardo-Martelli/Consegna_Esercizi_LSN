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

void exec(bool);

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


 exec(false);
 exec(true);

return 0;
}



/*============================================================*/



void exec(bool gaus = true){

if (gaus) cout << "Simulazione con Gauss" << endl;
else cout << " Simulazione con uniforme " << endl;

Random rnd;

ofstream write;
ofstream write2;
ofstream write3;
ofstream write4;
ofstream sums;

write.open("first.txt",ios::app);
write2.open("second.txt",ios::app);    
write3.open("third.txt",ios::app);
write4.open("first2.txt",ios::app);
    
sums.open("sums.txt");

    
write << "sum,mean,err" << endl;
write4 << "sum,mean,err" << endl;
write2 << "x,y,z" << endl;
write3 << "x,y,z" << endl;
sums << "sums,mean,err" << endl;
    
int steps= pow(10,8);
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
    
int eq = 4.0*pow(10,5);//iterations for equilibration
int eq_blcks = 4.0*pow(10,5);
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

    sums << sumr << "," << mean_prog/(j+1) << "," << error(mean_prog/(j+1),mean_prog2/(j+1),j)  << endl;

   
}
 
 cout << ".--------------------------------.  " << endl;
 cout << "| Equilibration psi_1 terminated | " << endl;
 cout << "'--------------------------------'  " << endl;
 cout << endl;


sumr = 0.0;
mean_prog=0.0;
mean_prog2=0.0;

for (int j= 0; j < B; j++){ 
    
    count = 0;
    
    for(int i = 0; i< F;i++){
        
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
          
       if(i%100 == 0)//riduco il numero di punti presi altrimenti sono troppi
              write2 << x << "," << y << "," << z << endl;
             
        }
          sumr += sqrt(x*x+y*y+z*z);
    
    }
    
    sumr /= double(F);
    mean_prog += sumr;
    mean_prog2 += sumr*sumr; 

    cout << " acc_rate batch "<< j+1 <<": " << ((double)count/F)*100. << "%" << endl;
    write  << sumr << "," << mean_prog/(j+1) << "," << error(mean_prog/(j+1),mean_prog2/(j+1),j)  << endl;
    
    //sumr = 0.0;
    
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

cout << ".------------------------.  " << endl;
cout << "| Equilibration of psi_2 | " << endl;
cout << "'------------------------'  " << endl;
cout << endl;

for (int j=0; j< eq;j++){
        count = 0;


    
	if(!gaus){
        double largh = 2.85;
        
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
	

		//double A = min(1., psi2(Tx,Ty,Tz)/psi2(x,y,z));
        double r = rnd.Rannyu(0.0,1.0);
	double T = psi2(Tx,Ty,Tz)/psi2(x,y,z);
	double A = minimo(1.0,T);

        if( r < A ){

          x = Tx;
          y = Ty;
          z = Tz;

          count ++;
             
        }
        
          sumr += sqrt(x*x+y*y+z*z);
	  loading(eq,j);

}
 
 cout << ".--------------------------------.  " << endl;
 cout << "| Equilibration psi_2 terminated | " << endl;
 cout << "'--------------------------------'  " << endl;
 cout << endl;

sumr = 0.0;
 for (int j= 0; j < B; j++){


    int count = 0;
     
    for(int i = 0; i< F;i++){

        Tx =  x + rnd.Rannyu( -largh, largh);
        Ty =  y + rnd.Rannyu( -largh, largh);
        Tz =  z + rnd.Rannyu( -largh, largh);

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
          if(i%100 == 0)//riduco il numero di punti presi altrimenti sono troppi
              write3 << x << "," << y << "," << z << endl;
              }
          sumr += sqrt(x*x+y*y+z*z);
        }
        
    
    sumr /= (double)(F);
    mean_prog += sumr;
    mean_prog2 += sumr*sumr; 
    
    cout << " acc_rate_2 batch "<< j+1 <<": " << ((double)count/F)*100. << "%" << endl;

  write4 << sumr << "," << mean_prog/(j+1) << "," << error(mean_prog/(j+1),mean_prog2/(j+1),j) << endl;

    }
   
write.close();  
write2.close();  
write3.close();    
write4.close();
sums.close();


}
