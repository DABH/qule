import sys
sys.path.append("/home/yuxin/anaconda3/envs/qiskit_env/lib/python3.11/site-packages")
import numpy as np
import re
import math
import time
import qiskit
from qiskit.circuit import QuantumCircuit
from qiskit.quantum_info import Operator, process_fidelity
from qiskit import transpile
from decimal import Decimal, getcontext



matrix_str = sys.argv[1:][0]
matrix_dim = int(sys.argv[1:][1])
nqubits = math.log2(matrix_dim)

# check the dimension, exit if not power of 2
if (not nqubits.is_integer()) :
    print("input matrix is not divisible into qubit subsystems")
    sys.exit(1)

nqubits = int(nqubits)

# Define a pattern to match each complex number
pattern = r"\(([+-]?\d+\.\d+)[eE]([+-]?\d+),([+-]?\d+\.\d+)[eE]([+-]?\d+)\)"
# Extract the real and imaginary parts 
matches = re.findall(pattern, matrix_str)
# convert to numpy array
arr = np.array([complex(float(m[0]+'e'+m[1]), float(m[2]+'e'+m[3])) for m in matches]).reshape(-matrix_dim, matrix_dim)
print("Matrix to be factorized:\n",arr)
op = Operator(arr)

qc = QuantumCircuit(nqubits)
qc.unitary(op, list(range(int(nqubits))), label="original gate")

print("\nOriginal circuit:")
print(qc.draw())
print("\n")


# ------------start decomposing------------------------

trans_qc = transpile(qc, basis_gates=['u3', 'cx'], optimization_level=3)
print(trans_qc.draw())