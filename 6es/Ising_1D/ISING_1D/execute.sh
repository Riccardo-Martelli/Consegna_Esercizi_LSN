rm seed.out
rm config.final

sed -i '8s/.*/1' input.dat
./run.sh

rm seed.out
rm config.final

sed -i '8s/.*/0' input.dat
./run.sh


