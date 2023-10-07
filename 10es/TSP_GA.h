#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include "../Random/random.h"
#include "../nifty_lib/lib.h"

using namespace std;

  

//========================================================//
//       CHROMOSOME                                       //
//========================================================//

#ifndef __CHROMOSOME__
#define __CHROMOSOME__

class chromosome {

public:
  // constructors
  chromosome();

  // destructor
  ~chromosome();
  
  void initialize();
  void print_chromo();
  void permutation(int);
  int get_chrom_pos(int);
  void set_allele(int,int);
  void pair_perm(int);
  void rev();
  void shift(int);

private:
 int num_genes = 50;
 static const int N = 50; 
 int fitness;
 
 int chrom[N] = {0};

 
};
#endif


//========================================================//
//     		  PROBLEM                                 //
//========================================================//

#ifndef __PROBLEM__
#define __PROBLEM__

class problem {

public:

  // constructors
  problem();

  // destructor
  ~problem();
  
   void read_cities(string);
   void print_city_file(string);
   double L2(chromosome &c);
   void crossover(chromosome &father, chromosome &mother, int);
   void crossover_best(chromosome &father, chromosome &mother, int);


private:

   static const int N_city = 50;
   double city[N_city][2]={{0,0}};
   
 
};
#endif


