rm equlibration_temp.dat output_etot.dat output_press.dat output_temp.dat output_ekin.dat output_epot.dat
touch equlibration_temp.dat output_etot.dat output_press.dat output_temp.dat output_ekin.dat output_epot.dat
./NVE_NVT.exe input.gas && ./NVE_NVT.exe input.liquid && ./NVE_NVT.exe input.solid
