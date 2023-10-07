#include <cmath>
#include <iostream>
#include <random>
#include <fstream>
#include <string>
#include <algorithm>
#include <numeric>
#include "../Random/random.h"
#include "../nifty_lib/lib.h"
#include "TSP_GA.h"
#include <map>

using namespace std;

int main(int argc, char **argv){

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

/*==============================================================================*/

int population = 10000;
int num_genes = 34;

chromosome c[population];
problem prob;

for (int i =0; i < population; i++){
	c[i].initialize();// così ogni volta lo inizializza uguale
	c[i].permutation(i);
}

int shape;

ifstream read;

read.open("Input.dat");

read >> shape;

read.close();

ofstream write;
ofstream arrays;
ofstream best_half;

if (shape == 0){
write.open("Best_L2.txt");
arrays.open("best_arrays.txt");
best_half.open("best_half_circ.txt");
}

else if (shape == 1){
write.open("Best_L2_square.txt");
arrays.open("best_arrays_square.txt");
best_half.open("best_half_square.txt");
}

write  << "L2,generation"  << endl;	
best_half << "best_half" << endl;

if (shape == 0){
	prob.circ_cities();
	prob.print_city_file("circ_city.txt");
}
else if (shape == 1){
	prob.square_cities();
	prob.print_city_file("square_city.txt");
}

map <double,int> L2map;
      
int generations = 50000;
double L2arr[population] = {0};
double bha[generations] = {0};

    
/*===========================EVOLUZIONE===============================*/
for(int j = 0; j < generations; j++){
    
map <double,int> L2map;
for(int i = 0; i < population; i++){
	 
	if (rnd.Rannyu() < 0.1) c[i].permutation(i);
	if (rnd.Rannyu() < 0.1) c[i].rev();
	if (rnd.Rannyu() < 0.01){
	
		int pos_shift = int(rnd.Rannyu(2.0,int(num_genes/2)) - 1);
		c[i].shift(pos_shift);

	}
	
	if (rnd.Rannyu() < 0.1){ 

		int pos_perm = int(rnd.Rannyu(2.0,num_genes-2));
		c[i].pair_perm(pos_perm);
		
	 	}

	 
	 L2arr[i] = prob.L2(c[i]);
	 L2map[prob.L2(c[i])] = i;
	 
	 
	 }

 
sort(L2arr,L2arr + population );

for(int i = 0; i < population/2; i++)
	bha[j] += L2arr[i]/double(population/2);
 


write << L2arr[0] << "," << j << endl;
best_half << bha[j] << endl; 

	for(int i = 0; i < population; i++){
	
		if (rnd.Rannyu() < 0.6)
			 prob.crossover_best(c[i],c[L2map[L2arr[0]]],17);
        }
    /*for(int i = 0; i < population-1; i++){
		if (rnd.Rannyu() < 0.1 && i != L2map[L2arr[0]])
			 prob.crossover(c[i],c[i+1],17);
			 
		}*/
			 
	for(int i = 0; i < num_genes-1; i++)
		arrays << c[L2map[L2arr[0]]].get_chrom_pos(i)+1 << "," ; 

	arrays << c[L2map[L2arr[0]]].get_chrom_pos(num_genes-1)+1 << endl;


	loading(generations,j);
}



write.close();
arrays.close();
best_half.close();
	
return 0;
}








