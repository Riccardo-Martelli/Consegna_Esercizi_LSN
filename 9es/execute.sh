sed -i '1s/.*/0/' Input.dat

make && echo " Cities on a circle" && ./main.exe

sed -i '1s/.*/1/' Input.dat

echo "Cities in a square"

./main.exe
