make && ./main.exe

sed -i '1s/.*/1/' Input.in

./main.exe

sed -i '1s/.*/0/' Input.in
