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
#include <mpi.h> 

using namespace std;

void evolution(int num_genes, int population, int generations,int rank){
	
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

/*==============================EVOLUTION======================================*/

	chromosome c[population];
	problem prob;
	prob.read_cities("American_capitals.dat");

	ofstream arrays;
	ofstream best_half;
	ofstream write;

	string arrays_name = "data_continent/cromosomes_";
	arrays_name+=to_string(rank);
	arrays.open(arrays_name);
	
	string best_half_name = "data_continent/L2_best_half_";
	best_half_name += to_string(rank);
	best_half.open(best_half_name);


	string write_name = "data_continent/L2_best_";
	write_name += to_string(rank);
	write.open(write_name);
	
	for(int i=0; i < population; i++ ){
		c[i].initialize();
		c[i].permutation(i);
	}

	double L2arr[population] = {0};
	double bha[generations] = {0};
	
	for(int j = 0; j < generations; j++){
	    
	map <double,int> L2map;

	for(int i = 0; i < population; i++){
		 
		if (rnd.Rannyu() < 0.05) c[i].permutation(i);
		if (rnd.Rannyu() < 0.05) c[i].rev();
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

	 
	sort(L2arr,L2arr + population);

	for(int i = 0; i < population/2; i++)
		bha[j] += L2arr[i]/double(population/2.0);
	 


	write << L2arr[0] << "," << j << endl;
	best_half << bha[j] << endl; 

		for(int i = 0; i < population; i++){
		
			if (rnd.Rannyu() < 0.6)
				 prob.crossover_best(c[i],c[L2map[L2arr[0]]],27);
				 
			}
				 
		for(int i = 0; i < num_genes-1; i++)
			arrays << c[L2map[L2arr[0]]].get_chrom_pos(i)+1 << "," ; 

		arrays << c[L2map[L2arr[0]]].get_chrom_pos(num_genes-1)+1 << endl;


		loading(generations,j);
	}

	write.close();
	arrays.close();
	best_half.close();
	
}

/*=================================PRINT PARALLEL=================================================*/

void print_parallel(string name,string header,double value,int rank,int flag = -1){

	name += "_" + to_string(rank);

	ofstream write;
	write.open(name,ofstream::app);

	//write << header << endl;
	if (flag == -1)
		write << value << endl;
	else if (flag != -1)
		write << value <<","<< flag << endl;
		

	write.close();
}

/*=================================PRINT PARALLEL CHROM============================================*/

void print_parallel_chrom(string name,chromosome c[],double L2arr[],map <double,int> L2map, int rank,int num_genes){

	ofstream arrays;
	name += "_" + to_string(rank);
	arrays.open(name,ofstream::app);

	for(int i = 0; i < num_genes-1; i++)
		arrays << c[L2map[L2arr[0]]].get_chrom_pos(i)+1 << "," ; 

	arrays << c[L2map[L2arr[0]]].get_chrom_pos(num_genes-1)+1 << endl;


	arrays.close();
}

/*==================================================================================================*/

int main(int argc, char **argv){


int population = 15000;
int num_genes = 50;
int N_migr = 50;
int generations = 15000;

int size,rank; 

MPI_Init(&argc,&argv);// as on the guide
MPI_Comm_size(MPI_COMM_WORLD, &size); //to know how many process
MPI_Comm_rank(MPI_COMM_WORLD, &rank); // to know who I am
MPI_Status stat;

//devo risettare il seed per ogni nuova simulazione altrimenti ho la rnd.SetPrimesCouple(rank);
// voglio settare io qunado migra o meno e devo anche fare un send-give.
bool mig;

ifstream migration;

migration.open("Input.dat");

 migration >> mig;

migration.close();

if (mig==0){
	
	evolution(num_genes, population, generations,rank);
	
	
	}
/*==================================================EVOLUTION WITH MIGRATION============================================================*/
else if (mig==1){


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

	// rnd.SetPrimesCouple(rank);
/*==============================EVOLUTION WITH MIGRATION======================================*/

	chromosome c[population];
	problem prob;
	prob.read_cities("American_capitals.dat");
	
	for(int i=0; i < population; i++ ){
		c[i].initialize();
		c[i].permutation(i);
	}


	double L2arr[population] = {0};
	double bha[generations] = {0};
	
	for(int j = 0; j < generations; j++){
	    
	map <double,int> L2map;
	
	int b = 0; //selects the best
	int sender = 0;
	int receiver = 1;
	int migrator[num_genes];
	int itag=1;
	
	for(int i = 0; i < population; i++){
		 
		 if(i%N_migr == 0){
		 
			int b = (int)(population*(1-pow(rnd.Rannyu(),20)))-1; //selects the best
			int sender = (int)rnd.Rannyu(0,size); ;
			int receiver = (int)rnd.Rannyu(0,size); ;
			
		   	if(rank == 0){
				
				sender = (int)rnd.Rannyu(0,size);            
				receiver = (int)rnd.Rannyu(0,size);           
				while(receiver==sender){                      
				    receiver = (int)rnd.Rannyu(0,size);      
				}
	//			cout << endl << "Migration (b, sender, receiver): " << b << ", " << sender << ", " << receiver << endl;
			}
		 
		 }
		 
		 MPI_Bcast(&b, 1, MPI_INT, 0, MPI_COMM_WORLD); // il core 0 manda a tutti gli altri
		 MPI_Bcast(&sender, 1, MPI_INT, 0, MPI_COMM_WORLD); 
		 MPI_Bcast(&receiver, 1, MPI_INT, 0, MPI_COMM_WORLD); 
		 
            // se sono il sender, salvo su migrator e lo spedisco al receiver
		if(rank == sender){
                	for(int i = 0; i< num_genes; i++){
                    		migrator[i] = c[b].get_chrom_pos(i);//pop.Chr[h].GetGen(i);
                	}
               		MPI_Send(migrator,num_genes,MPI_INT,receiver,itag,MPI_COMM_WORLD);

            	}
            // se sono il ricevitore lo salvo
		if(rank == receiver){
	           
	            MPI_Recv(migrator,num_genes,MPI_INT,sender,itag,MPI_COMM_WORLD,&stat);
	            for(int i = 0; i< num_genes; i++)
	           	 c[b].set_allele(i,migrator[i]);
	           	 //COSA SUCCEDE SE METTO DUE CROMOSOMI UGUALI??
            	}
         
		 
		 /*=========MUTATIONS=================*/
		 
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

	 
	sort(L2arr,L2arr + population);

	for(int i = 0; i < population/2; i++)
		bha[j] += L2arr[i]/double(population/2.0);
	 
\

	print_parallel("data_parallel/L2_best","L2,genetion",L2arr[0],rank,j);
	print_parallel("data_parallel/L2_best_half","L2BH,",bha[j],rank);

		for(int i = 0; i < population; i++){
		
			if (rnd.Rannyu() < 0.6)
				 prob.crossover_best(c[i],c[L2map[L2arr[0]]],26);
			if (rnd.Rannyu() < -0.1)
				 prob.crossover_best(c[i],c[L2map[L2arr[1]]],26);
				 
			}
				 
		print_parallel_chrom("data_parallel/cromosomes",c,L2arr,L2map,rank,num_genes);


		loading(generations,j);
	}

}	
	
	
MPI_Finalize();

return 0;
}

