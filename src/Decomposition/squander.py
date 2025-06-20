# !pip install qiskit
# !pip install git+https://github.com/LNoorl/qiskit-terra.git@feature/sk-pass
import sys
import os
import numpy as np
import re
import math
from qiskit.circuit import QuantumCircuit
from qiskit.quantum_info import Operator

import sys  
sys.path.append('/home/rhea_huang/anaconda3/envs/qgd/lib/python3.10/site-packages')

#decomposition classes
from qgd_python.decomposition.qgd_N_Qubit_Decomposition_adaptive import qgd_N_Qubit_Decomposition_adaptive as N_Qubit_Decomposition_adaptive
from qgd_python.decomposition.qgd_N_Qubit_Decomposition_custom import qgd_N_Qubit_Decomposition_custom as N_Qubit_Decomposition_custom
from qgd_python.decomposition.qgd_N_Qubit_Decomposition import qgd_N_Qubit_Decomposition as N_Qubit_Decomposition


matrix_str = sys.argv[1:][0]
matrix_dim = int(sys.argv[1:][1])
qubit = math.log2(matrix_dim)

# check the dimension, exit if not power of 2
if (not qubit.is_integer()) :
    print("input matrix is not divisible into qubit subsystems")
    sys.exit(1)


# Define a pattern to match each complex number
pattern = r"\(([+-]?\d+\.\d+)[eE]([+-]?\d+),([+-]?\d+\.\d+)[eE]([+-]?\d+)\)"
# Extract the real and imaginary parts 
matches = re.findall(pattern, matrix_str)
# convert to numpy array
arr = np.array([complex(float(m[0]+'e'+m[1]), float(m[2]+'e'+m[3])) for m in matches]).reshape(-matrix_dim, matrix_dim)
print("Matrix to be factorized:\n",arr)
op = Operator(arr)

qc = QuantumCircuit(qubit)
qc.unitary(op, list(range(int(qubit))), label="original gate")

print("\nOriginal circuit:")
print(qc.draw())
print("\n")

# start decomposing....

# creating a class to decompose the 
cDecompose = N_Qubit_Decomposition(np.conj(arr).T )

# setting the verbosity of the decomposition
cDecompose.set_Verbose( 3 )

# setting the debugfile name. If it is not set, the program will not debug. 
# cDecompose.set_Debugfile("debug.txt")

# setting the tolerance of the optimization process. The final error of the decomposition would scale with the square root of this value.
cDecompose.set_Optimization_Tolerance( 1e-12 )

# set the number of block to be optimized in one shot
cDecompose.set_Optimization_Blocks( 20 )

# starting the decomposition
cDecompose.Start_Decomposition()

# list the decomposing operations
cDecompose.List_Gates()

# get the decomposing operations
quantum_circuit = cDecompose.get_Quantum_Circuit()

print("\nDiscretized circuit:")
# print the quantum circuit
print(quantum_circuit)
print("Error:", np.linalg.norm(Operator(quantum_circuit).data - Operator(arr).data))
