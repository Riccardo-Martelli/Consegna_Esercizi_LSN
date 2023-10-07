#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <string>
#include "lib.h"

using namespace std;

char st[11] = {' ',' ',' ',' ',' ',' ',' ',' ',' ',' ','\0'};

void loading(double max, double j){

    
       if((max-j)/(max/10.) - int((max-j)/(max/10.)) == 0){
          st[int(10 - int((max-j)/(max/10.)))-1] = '='; 
          cout <<'\r' <<"["<< st << "]  " << j <<" out of " << max << " " << double(j/max*100) << "%" << flush;
       
       }  
        
       else if(j==max-1){
	st[9]= '=';
        cout <<'\r' <<"["<< st << "]  " << j+1 <<" out of " << max << " " << double((j+1)/max*100) << "%" << flush;
        for(int i =0; i<11; i++)
          st[i] = ' ';
       	
	cout << endl; 
       }
}


double error(double AV, double AV2, int n){
   if(n==0){
      return 0;
   }
   else{
      return sqrt((AV2-AV*AV)/n);
   }

}

