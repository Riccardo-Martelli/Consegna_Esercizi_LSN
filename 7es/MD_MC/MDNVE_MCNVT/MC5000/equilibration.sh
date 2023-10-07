rm equlibration_temp.dat equlibration_pot.dat output_etot.dat output_press.dat output_temp.dat output_ekin.dat output_epot.dat output_gdr_MC.dat output_gdr_MD.dat output_gdr_final.dat

touch equlibration_temp.dat equlibration_pot.dat output_etot.dat output_press.dat output_temp.dat output_ekin.dat output_epot.dat output_gdr_MC.dat output_gdr_MD.dat 

./NVE_NVT.exe input.gas && ./NVE_NVT.exe input.liquid && ./NVE_NVT.exe input.solid
