#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include "../Random/random.h"
#include "../nifty_lib/lib.h"
#include "TSP_GA.h"
#include <bits/stdc++.h>
#include <numeric>
#include <map>

using namespace std;

//========================================================//
//       CHROMOSOME                                       //
//========================================================//

chromosome :: chromosome(){} // costruttore
chromosome :: ~chromosome(){} // distruttore

void chromosome :: initialize(){ 
 
 chrom[0] = 1;
 
 for (int i = 1; i< num_genes; i++){
	
	chrom[i] = i+1;
	
 }
 

}

void chromosome :: permutation(int iterations){

for (int i = 0; i < iterations; i++)
	next_permutation(chrom+1,chrom+num_genes-1); // parametro per iterare da aggiungere

}

void chromosome :: print_chromo(){

	cout << "[" << chrom[0];	
	for(int i = 1; i <num_genes; i++)	
		cout <<","<< chrom[i] ;
	cout << "]" << endl;
	cout <<"length: "<< sizeof(chrom)/sizeof(*chrom) << endl;

}

int chromosome :: get_chrom_pos(int i){

	if (i < num_genes && i >=0) return int(chrom[i]-1);
	else {
	
		cerr << "Index out of bounds" << endl;
		cerr << "The value of the index is: " << i << endl;
		exit(1);
	}
}

void chromosome :: set_allele(int i,int value){

	chrom[i] = value;

}

void chromosome :: pair_perm(int pos1) {

	next_permutation(chrom+pos1,chrom+pos1+2);
	
}

void chromosome :: rev(){

	reverse(chrom + 1, chrom +num_genes-1);

}

void chromosome :: shift(int m){

  int copy[num_genes] ={0};
  copy[0] = 1;
  if(m < num_genes/2){	
     for(int i = 1; i < num_genes; i++){	
	if( i + m > num_genes-1)	
		copy[(i+m)%(num_genes-1)] = chrom[i];
		
	else copy[i+m] = chrom[i];
	  }

	 memcpy(chrom, copy, sizeof copy);

	 }
  else cout << " m too big " << endl;

}

//========================================================//
//  	       PROBLEM                                    //
//========================================================//

  problem :: problem(){};

  // destructor
  problem :: ~problem(){};
  

void problem :: circ_cities(){

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

	double theta = rnd.Rannyu(0.0,2*M_PI);
	double x = cos(theta);
	double y = sin(theta);
	
	for(int i= 0; i < N_city;i++){
	
	theta = rnd.Rannyu(0.0,2*M_PI);
	x = cos(theta);
	y = sin(theta);
	city[i][0]= x;
	city[i][1]= y;
	
	}
}

void problem :: square_cities(){

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

	double x;
	double y;
	
	for(int i= 0; i < N_city;i++){
	
		x = rnd.Rannyu(0.0,1.0);
		y = rnd.Rannyu(0.0,1.0);
		city[i][0]= x;
		city[i][1]= y;
	
	}
}


void problem :: print_city_file(string name){

	ofstream write;
	write.open(name);
	
	write  << "x,y"  << endl;
	
	for (int i = 0; i < N_city; i++){
		
		write << city[i][0] << "," << city[i][1] << endl;
		
	}
	write.close();

}


double problem :: L2(chromosome &c) {

	double sum = 0.0; 

	for(int i = 0;i < N_city-1; i++)
		sum += pow(city[c.get_chrom_pos(i)][0] - city[c.get_chrom_pos(i+1)][0] ,2.0) + pow(city[c.get_chrom_pos(i)][1] - city[c.get_chrom_pos(i+1)][1],2.0); 


	return pow(city[c.get_chrom_pos(N_city-1)][0] - city[c.get_chrom_pos(0)][0] ,2.0) + pow(city[c.get_chrom_pos(N_city-1)][1] - city[c.get_chrom_pos(0)][1],2.0) + sum;

}


void problem :: crossover(chromosome &father, chromosome &mother, int it){

	map <int,int> map_ord;
	map <int,int> map_ord2;
	
	int appo[it]={0};
	int appo2[it]={0};
	int position[it]={0};
	int position2[it]={0};

	if(it >= N_city){
		cerr << " Cut too long " << endl;
		exit(0);
		}

	//come father utilizzo il best, e mother gli altri, in verità il crossover è su entrambi, fare un crossover che non intacca il best

	for (int i = 0; i < it; i++){
	
		appo[i] = father.get_chrom_pos(N_city-i-1) + 1;
		appo2[i] = mother.get_chrom_pos(N_city-i-1) + 1;
	}
	
	for (int i = 0; i < it; i++){
		for (int j = 0; j < N_city; j++){ //Qui devo vedere tutto ilcromosoma
			if((mother.get_chrom_pos(j) + 1) == appo[i])
				position[i] = -j+N_city-1;
			if((father.get_chrom_pos(j) + 1) == appo2[i])
				position2[i] = -j+N_city-1;
		}
		
		map_ord[position[i]] = appo[i];
		map_ord2[position2[i]] = appo2[i];
		
	}

	int count = N_city-1;
	for(auto x: map_ord){ // 

		father.set_allele(count,int(x.second));
		count --;
		
		}
	count = N_city-1;
	for(auto x: map_ord2){
	
		mother.set_allele(count,int(x.second));
		count --;
		
		}
}

void problem :: crossover_best(chromosome &father, chromosome &mother, int it){

	map <int,int> map_ord;
	int appo[it]={0};
	int position[it]={0};

	if(it >= N_city){
		cerr << " Cut too long " << endl;
		exit(0);
		}

	//come father utilizzo il best, e mother gli altri, in verità il crossover è su entrambi, fare un crossover che non intacca il best

	for (int i = 0; i < it; i++){
	
		appo[i] = father.get_chrom_pos(N_city-i-1) + 1;

	}

	for (int i = 0; i < it; i++){
		for (int j = 0; j < N_city; j++){ //Qui devo vedere tutto ilcromosoma
			if((mother.get_chrom_pos(j) + 1) == appo[i])
				position[i] = -j+N_city-1;
			
		}
		
		map_ord[position[i]] = appo[i];
		
	}


	int count = N_city-1;
	for(auto x: map_ord){ // 

		father.set_allele(count,int(x.second));
		count --;
		
		}	
	
}






