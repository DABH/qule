#! /usr/bin/bash

#Script to automate running tests, tests for plot scripts included, uncomment and edit as needed
#cd Build/Release
#chmod a+x ./qule

#Number of iterations for specific examples
#./qule --matrix_type Had --method GDP --n 2 --num_tests 100 --alpha 0.1 
#./qule --matrix_type Had --method "GDP*" --n 2 --num_tests 100 --alpha 0.1
#./qule --matrix_type Had --method GDLM --n 2 --num_tests 100 --alpha 0.1
#./qule --matrix_type Had --method "GDLM*" --n 2 --num_tests 100 --alpha 0.1
#./qule --matrix_type Had --method DLR --n 2 --num_tests 100 --alpha 0.1
#./qule --matrix_type Had --method "DLR*" --n 2 --num_tests 100 --alpha 0.1
#./qule --matrix_type Had --method NM --n 2 --num_tests 100 --alpha 0.1
#./qule --matrix_type QFT --method GDP --n 4 --num_tests 100 --alpha 0.1 
#./qule --matrix_type QFT --method "GDP*" --n 4 --num_tests 100 --alpha 0.1
#./qule --matrix_type QFT --method GDLM --n 4 --num_tests 100 --alpha 0.1
#./qule --matrix_type QFT --method "GDLM*" --n 4 --num_tests 100 --alpha 0.1
#./qule --matrix_type QFT --method DLR --n 4 --num_tests 100 --alpha 0.1
#./qule --matrix_type QFT --method "DLR*" --n 4 --num_tests 100 --alpha 0.1
#./qule --matrix_type QFT --method NM --n 4 --num_tests 100 --alpha 0.1
#./qule --matrix_type QFT --method GDP --n 8 --num_tests 100 --alpha 0.1 
#./qule --matrix_type QFT --method "GDP*" --n 8 --num_tests 100 --alpha 0.1
#./qule --matrix_type QFT --method GDLM --n 8 --num_tests 100 --alpha 0.1
#./qule --matrix_type QFT --method "GDLM*" --n 8 --num_tests 100 --alpha 0.1
#./qule --matrix_type QFT --method DLR --n 8 --num_tests 100 --alpha 0.1
#./qule --matrix_type QFT --method "DLR*" --n 8 --num_tests 100 --alpha 0.1
#./qule --matrix_type QFT --method NM --n 8 --num_tests 100 --alpha 0.1
#./qule --matrix_type grover --method GDP --n 8 --num_tests 100 --alpha 0.1 
#./qule --matrix_type grover --method "GDP*" --n 8 --num_tests 100 --alpha 0.1
#./qule --matrix_type grover --method GDLM --n 8 --num_tests 100 --alpha 0.1
#./qule --matrix_type grover --method "GDLM*" --n 8 --num_tests 100 --alpha 0.1
#./qule --matrix_type grover --method DLR --n 8 --num_tests 100 --alpha 0.1
#./qule --matrix_type grover --method "DLR*" --n 8 --num_tests 100 --alpha 0.1
#./qule --matrix_type grover --method NM --n 8 --num_tests 100 --alpha 0.1
#./qule --matrix_type grover --method GDP --n 16 --num_tests 100 --alpha 0.1 
#./qule --matrix_type grover --method "GDP*" --n 16 --num_tests 100 --alpha 0.1
#./qule --matrix_type grover --method GDLM --n 16 --num_tests 100 --alpha 0.1
#./qule --matrix_type grover --method "GDLM*" --n 16 --num_tests 100 --alpha 0.1
#./qule --matrix_type grover --method DLR --n 16 --num_tests 100 --alpha 0.1
#./qule --matrix_type grover --method "DLR*" --n 16 --num_tests 100 --alpha 0.1
#./qule --matrix_type grover --method NM --n 16 --num_tests 100 --alpha 0.1

#Random matrix, wall time vs size
./qule --matrix_type rand --method GDP --n 4 --num_tests 100 --alpha 0.25
./qule --matrix_type rand --method GDP --n 8 --num_tests 100 --alpha 0.25 
./qule --matrix_type rand --method GDP --n 12 --num_tests 100 --alpha 0.25 
./qule --matrix_type rand --method GDP --n 16 --num_tests 100 --alpha 0.25
./qule --matrix_type rand --method GDP --n 20 --num_tests 100 --alpha 0.25  
./qule --matrix_type rand --method GDP --n 24 --num_tests 100 --alpha 0.25
./qule --matrix_type rand --method GDP --n 28 --num_tests 100 --alpha 0.25 
./qule --matrix_type rand --method GDP --n 32 --num_tests 100 --alpha 0.25
./qule --matrix_type rand --method GDP* --n 4 --num_tests 100 --alpha 0.25 
./qule --matrix_type rand --method GDP* --n 8 --num_tests 100 --alpha 0.25 
./qule --matrix_type rand --method GDP* --n 12 --num_tests 100 --alpha 0.25 
./qule --matrix_type rand --method GDP* --n 16 --num_tests 100 --alpha 0.25
./qule --matrix_type rand --method GDP* --n 20 --num_tests 100 --alpha 0.25  
./qule --matrix_type rand --method GDP* --n 24 --num_tests 100 --alpha 0.25
./qule --matrix_type rand --method GDP* --n 28 --num_tests 100 --alpha 0.25 
./qule --matrix_type rand --method GDP* --n 32 --num_tests 100 --alpha 0.25  
./qule --matrix_type rand --method GDLM --n 4 --num_tests 100 --alpha 0.25
./qule --matrix_type rand --method GDLM --n 8 --num_tests 100 --alpha 0.25
./qule --matrix_type rand --method GDLM --n 12 --num_tests 100 --alpha 0.25 
./qule --matrix_type rand --method GDLM --n 16 --num_tests 100 --alpha 0.25
./qule --matrix_type rand --method GDLM --n 20 --num_tests 100 --alpha 0.25  
./qule --matrix_type rand --method GDLM --n 24 --num_tests 100 --alpha 0.25
./qule --matrix_type rand --method GDLM --n 28 --num_tests 100 --alpha 0.25 
./qule --matrix_type rand --method GDLM --n 32 --num_tests 100 --alpha 0.25 
./qule --matrix_type rand --method GDLM* --n 4 --num_tests 100 --alpha 0.125
./qule --matrix_type rand --method GDLM* --n 8 --num_tests 100 --alpha 0.125
./qule --matrix_type rand --method GDLM* --n 12 --num_tests 100 --alpha 0.125
./qule --matrix_type rand --method GDLM* --n 16 --num_tests 100 --alpha 0.125
./qule --matrix_type rand --method GDLM* --n 20 --num_tests 100 --alpha 0.125  
./qule --matrix_type rand --method GDLM* --n 24 --num_tests 100 --alpha 0.125
./qule --matrix_type rand --method GDLM* --n 28 --num_tests 100 --alpha 0.125 
./qule --matrix_type rand --method GDLM* --n 32 --num_tests 100 --alpha 0.125 
./qule --matrix_type rand --method DLR --n 4 --num_tests 100 --alpha 1 
./qule --matrix_type rand --method DLR --n 8 --num_tests 100 --alpha 1 
./qule --matrix_type rand --method DLR --n 12 --num_tests 100 --alpha 1 
./qule --matrix_type rand --method DLR --n 16 --num_tests 100 --alpha 1
./qule --matrix_type rand --method DLR --n 20 --num_tests 100 --alpha 1  
./qule --matrix_type rand --method DLR --n 24 --num_tests 100 --alpha 1
./qule --matrix_type rand --method DLR --n 28 --num_tests 100 --alpha 1 
./qule --matrix_type rand --method DLR --n 32 --num_tests 100 --alpha 1
./qule --matrix_type rand --method DLR* --n 4 --num_tests 100 --alpha 0.5
./qule --matrix_type rand --method DLR* --n 8 --num_tests 100 --alpha 0.5
./qule --matrix_type rand --method DLR* --n 12 --num_tests 100 --alpha 0.5
./qule --matrix_type rand --method DLR* --n 16 --num_tests 100 --alpha 0.5
./qule --matrix_type rand --method DLR* --n 20 --num_tests 100 --alpha 0.5
./qule --matrix_type rand --method DLR* --n 24 --num_tests 100 --alpha 0.5
./qule --matrix_type rand --method DLR* --n 28 --num_tests 100 --alpha 0.5
./qule --matrix_type rand --method DLR* --n 32 --num_tests 100 --alpha 0.5
./qule --matrix_type rand --method NM --n 4 --num_tests 100 --alpha 1 
./qule --matrix_type rand --method NM --n 8 --num_tests 100 --alpha 1 
./qule --matrix_type rand --method NM --n 12 --num_tests 100 --alpha 1 
./qule --matrix_type rand --method NM --n 16 --num_tests 100 --alpha 1
./qule --matrix_type rand --method NM --n 20 --num_tests 100 --alpha 1  
./qule --matrix_type rand --method NM --n 24 --num_tests 100 --alpha 1
./qule --matrix_type rand --method NM --n 28 --num_tests 100 --alpha 1 
./qule --matrix_type rand --method NM --n 32 --num_tests 100 --alpha 1

#Number of iterations vs num training examples; 4 qubits, ie n=16
#./qule --matrix_type rand --method GDP* --n 16 --num_x 8 --num_tests 100 --alpha 0.1
#./qule --matrix_type rand --method GDP* --n 16 --num_x 12 --num_tests 100 --alpha 0.1
#./qule --matrix_type rand --method GDP* --n 16 --num_x 16 --num_tests 100 --alpha 0.1
#./qule --matrix_type rand --method GDP* --n 16 --num_x 20 --num_tests 100 --alpha 0.1
#./qule --matrix_type rand --method GDP* --n 16 --num_x 24 --num_tests 100 --alpha 0.1
#./qule --matrix_type rand --method GDLM* --n 16 --num_x 8 --num_tests 100 --alpha 0.1
#./qule --matrix_type rand --method GDLM* --n 16 --num_x 12 --num_tests 100 --alpha 0.1
#./qule --matrix_type rand --method GDLM* --n 16 --num_x 16 --num_tests 100 --alpha 0.1
#./qule --matrix_type rand --method GDLM* --n 16 --num_x 20 --num_tests 100 --alpha 0.1
#./qule --matrix_type rand --method GDLM* --n 16 --num_x 24 --num_tests 100 --alpha 0.1
#./qule --matrix_type rand --method DLR* --n 16 --num_x 8 --num_tests 100 --alpha 0.1
#./qule --matrix_type rand --method DLR* --n 16 --num_x 12 --num_tests 100 --alpha 0.1
#./qule --matrix_type rand --method DLR* --n 16 --num_x 16 --num_tests 100 --alpha 0.1
#./qule --matrix_type rand --method DLR* --n 16 --num_x 20 --num_tests 100 --alpha 0.1
#./qule --matrix_type rand --method DLR* --n 16 --num_x 24 --num_tests 100 --alpha 0.1

#sequential vs batch DLR*
#./qule --matrix_type rand --method DLR* --n 8 --num_tests 100 --alpha 0.1 --each_mat --sequential
#./qule --matrix_type rand --method DLR* --n 16 --num_tests 100 --alpha 0.1 --each_mat --sequential
#./qule --matrix_type rand --method DLR* --n 32 --num_tests 100 --alpha 0.1 --each_mat --sequential
#./qule --matrix_type rand --method DLR* --n 8 --num_tests 100 --alpha 0.1 --each_mat 
#./qule --matrix_type rand --method DLR* --n 16 --num_tests 100 --alpha 0.1 --each_mat 
#./qule --matrix_type rand --method DLR* --n 32 --num_tests 100 --alpha 0.1 --each_mat 
