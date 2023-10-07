<<com
	comments
com

rm -r data_parallel && mkdir data_parallel

sed -i '1s/.*/0/' Input.dat

make && mpiexec -np 4 ./main.exe


sed -i '1s/.*/1/' Input.dat

mpiexec -np 4 ./main.exe 


