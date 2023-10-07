 sed -i '8s/.*/1/' input.gas
 sed -i '8s/.*/1/' input.liquid
 sed -i '8s/.*/1/' input.solid

 sed -i '9s/.*/30000/' input.gas
 sed -i '9s/.*/30000/' input.liquid
 sed -i '9s/.*/30000/' input.solid

 make

 ./equilibration.sh
 
 ./over.exe

 cp equlibration_pot.dat pot_equilibration.dat
 cp equlibration_temp.dat temp_equilibration.dat

 sed -i '8s/.*/50000/' input.gas
 sed -i '8s/.*/50000/' input.liquid
 sed -i '8s/.*/50000/' input.solid

 sed -i '9s/.*/10/' input.gas
 sed -i '9s/.*/10/' input.liquid
 sed -i '9s/.*/10/' input.solid

 ./equilibration.sh
