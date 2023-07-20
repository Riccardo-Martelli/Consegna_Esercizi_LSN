rm output.* *.0 *.0.*

sed -i '5s/.*/0/' input.dat
for i in $(seq 0.5 0.25 3.0)
do
	t=$i
	sed -i '1s/.*/'${t}'/' input.dat
	sed -i '4s/.*/0.0/' input.dat

	make && ./Monte_Carlo_ISING_1D.exe

	sed -i '4s/.*/0.02/' input.dat

	./Monte_Carlo_ISING_1D.exe
done

sed -i '5s/.*/1/' input.dat
for i in $(seq 0.5 0.25 3.0)
do
	t=$i
	sed -i '1s/.*/'${t}'/' input.dat
	sed -i '4s/.*/0.0/' input.dat

	./Monte_Carlo_ISING_1D.exe

	sed -i '4s/.*/0.02/' input.dat

	./Monte_Carlo_ISING_1D.exe

done



