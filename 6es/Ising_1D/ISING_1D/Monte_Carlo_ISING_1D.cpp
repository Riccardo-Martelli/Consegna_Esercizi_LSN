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
#include <ostream>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <sstream>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main()
{ 
  
  Input(); //Inizialization
  
  if(metro ==1){
    cout << ".--------------------------------------. "<< endl;
    cout << "| Simulazione con algoritmo Metropolis |"<< endl;
    cout << "| h ="<<h<<"                               |" << endl;
    cout << "'______________________________________'"<< endl;
    cout << endl;
  }

  else{
    cout << ".----------------------------------. "<< endl;
    cout << "| Simulazione con Gibbs' sampling  |"<< endl;
    cout << "| h = "<<h<<"                          |" << endl;
    cout << "'__________________________________'"<< endl;
    cout << endl;
    }

  int max = 2000;
  if(!mix){
    ofstream Ene;
    Ene.open("output.eq",ios::app);

  cout << " equilibrando... " << endl;
   for (int i = 0; i< max; i++){
    Move(metro);
    Measure();
    Accumulate();
    
    stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy
    glob_av[iu]  += stima_u;
    glob_av2[iu] += stima_u*stima_u;

    err_u=Error(glob_av[iu],glob_av2[iu],i);
    Ene << i <<  "," << stima_u << "," << glob_av[iu]/(double)i  << endl;

     }

    Ene.close();
   }


  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move(metro);
      Measure();
      Accumulate(); //Update block averages
    }
    Averages(iblk);   //Print results for current block
  }
  Stampa();
  ConfFinal(); //Write final configuration

  return 0;
}


void Input(void)
{
  ifstream ReadInput;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;
    
  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;

  ReadInput >> mix;

  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();


//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
 
  n_props = 4; //Number of observables
//bisogna ri inizializzare gli spin
//initial configuration
 if(!mix){
  for (int i=0; i<nspin; ++i)
  {
    if(rnd.Rannyu() >= 0.5) s[i] = 1;
    else s[i] = -1;
   }
  }
 else {
	 ifstream ReadConfig;
    ReadConfig.open("config.out");

    for (int i=0; i<nspin; ++i)
    {
      ReadConfig >> s[i];
    } // leggo configurazione >> temp equilibrata
    }
//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
}


void Move(int metro)
{
  int o;
  //double p, energy_old, energy_new, sm;
  //double energy_up, energy_down;

  for(int i=0; i<nspin; ++i)
  {
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);

    if(metro==1) //Metropolis
    {
  // propongo il flip
      int s_try = -s[o];
      attempted++;

      // calcolo variazione energia DeltaE
      double DeltaE = 2*Boltzmann(s_try, o); // <<< 

      // cuore algoritmo metropolis
      double q = exp(-beta*DeltaE);
      double A = min(1.,q);

      if (A==1) // accetto direttamente
      {
          s[o] = s_try; 
          accepted++;
      }
      else // accetto con probabilità A
      {
        if(rnd.Rannyu() < A)
        {
          s[o] = s_try; 
          accepted++;
        }
        // << qui è la differenza da gibbs: potrei non accettare
      }
        
    }
    else //Gibbs sampling
    {

      attempted++;
      
      
      //----------------------------------
      // estraggo spin random
      int snew = ((int)rnd.Rannyu(0,2))*2-1;

      // calcolo variazione energia DeltaE e probabilità associata
      double DeltaE = -2*Boltzmann(snew, o); 
      double p = 1.0/(1.0+exp(-beta*DeltaE));

      // impongo il flip con probabilità p
      if(rnd.Rannyu() < p)
      {
        s[o] = snew;
      } else s[o] = -snew;
      //-------------------------------/

      accepted++;
    } 
  }    
}

double Boltzmann(int sm, int ip)
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}

void Measure()
{
  int bin;
  double u = 0.0, m = 0.0;
  double H = 0.0;

//cycle over spins
  for (int i=0; i<nspin; ++i)
  {
     u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
// INCLUDE YOUR CODE HERE
    H = -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
   // u += H; 
    //c += H*H; // prima era sbagliato!
    m += s[i];
    // x +=  s[i]*s[i]; // non serve, perchè h = 0
  
  } 
      
  
  walker[iu] = u;
// INCLUDE YOUR CODE HERE
  walker[ic] = u*u; // 
   //walker[ic] = beta*beta*(c/(double)nspin-u*u/pow((double)nspin,2)); // prima era sbagliato!
   // magnetization
  walker[im] = m;                 
   // susceptivity
  walker[ix] = beta*(m*m); // non beta*(x - m*m), perché h = 0    
   
}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{
    
   ofstream Ene, Heat, Mag, Chi;
    
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    Ene.open("output.ene.0",ios::app);
    stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy
    glob_av[iu]  += stima_u;
    glob_av2[iu] += stima_u*stima_u;

    err_u=Error(glob_av[iu],glob_av2[iu],iblk);
    Ene << iblk <<  "," << stima_u << "," << glob_av[iu]/(double)iblk << "," << err_u << endl;
    Ene.close();
    
        // CAPACITÀ TERMICA
    Heat.open("output.heat.0",ios::app);
    stima_c = beta*beta * (blk_av[ic]/blk_norm - pow(blk_av[iu]/blk_norm, 2))/(double)nspin; 
    //stima_c = blk_av[ic]/blk_norm; // prima era sbagliato!
    glob_av[ic]  += stima_c; // accu globale
    glob_av2[ic] += stima_c*stima_c; // accu quadratico

    err_c=Error(glob_av[ic],glob_av2[ic],iblk);
    Heat << iblk  <<  "," << setprecision(9)<< stima_c << "," << setprecision(9)<< glob_av[ic]/(double)iblk << "," << setprecision(9)<< err_c << endl;
    Heat.close();

    // MAGNETIZZAZIONE
    Mag.open("output.mag.0",ios::app);
    stima_m = blk_av[im]/blk_norm/(double)nspin; 
    glob_av[im]  += stima_m; // accu globale
    glob_av2[im] += stima_m*stima_m; // accu quadratico

    err_m=Error(glob_av[im],glob_av2[im],iblk);
    Mag << iblk <<  "," << stima_m << "," << glob_av[im]/(double)iblk << "," << err_m << endl;
    Mag.close();

    // SUSCETTIVITÀ
    Chi.open("output.chi.0",ios::app);
    stima_x = blk_av[ix]/blk_norm/(double)nspin;
    glob_av[ix]  += stima_x; // accu globale
    glob_av2[ix] += stima_x*stima_x; // accu quadratico

    err_x=Error(glob_av[ix],glob_av2[ix],iblk);
    Chi << iblk <<  "," << stima_x << "," << glob_av[ix]/(double)iblk << "," << err_x << endl;
    Chi.close();

    cout << "----------------------------" << endl << endl;

// INCLUDE YOUR CODE HERE
   }
void Stampa(){
   ofstream Ene, Heat, Mag, Chi;
  if(h==0)
  {
    // ENERGY
    Ene.open("ene.0",ios::app);
    Ene << temp << "," << glob_av[iu]/(double)nblk << "," << err_u << endl;
    Ene.close();

    // HEAT CAPACITY
    Heat.open("heat.0",ios::app);
    Heat << temp << "," << glob_av[ic]/(double)nblk << "," << err_c << endl;
    Heat.close();

    // MAGN SUSCEPTIBILITY
    Chi.open("chi.0",ios::app);
    Chi << temp << "," << glob_av[ix]/(double)nblk << "," << err_x << endl;
    Chi.close();
  }

  else if(h!=0)
  {
    // MAGNETIZATION
    stringstream stream;
    stream << fixed << std::setprecision(3) << h;
    string hh = stream.str();
    //string hh = to_string(h);
    Mag.open("mag."+hh,ios::app);
    Mag << temp << "," << glob_av[im]/(double)nblk << "," << err_m << "," << h << endl;
    Mag.close();
  }

}


void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    if(iblk==1) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
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
