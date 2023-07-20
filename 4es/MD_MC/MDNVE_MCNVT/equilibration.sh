 sed -i '8s/.*/1/' input.gas
 sed -i '8s/.*/1/' input.liquid
 sed -i '8s/.*/1/' input.solid

 sed -i '9s/.*/30000/' input.gas
 sed -i '9s/.*/30000/' input.liquid
 sed -i '9s/.*/30000/' input.solid

 ./execute.sh
 
 ./over.exe

 cp equlibration_temp.dat temp_equilibration.dat

 sed -i '8s/.*/60/' input.gas
 sed -i '8s/.*/60/' input.liquid
 sed -i '8s/.*/60/' input.solid

 sed -i '9s/.*/2000/' input.gas
 sed -i '9s/.*/2000/' input.liquid
 sed -i '9s/.*/2000/' input.solid

 ./execute.sh
