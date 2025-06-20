# !pip install qiskit
# !pip install git+https://github.com/LNoorl/qiskit-terra.git@feature/sk-pass
import sys
import numpy as np
import re
import math
from qiskit.circuit import QuantumCircuit
from qiskit.circuit.library import SGate, HGate, SdgGate, TGate, TdgGate, ZGate, XGate, YGate
from qiskit.transpiler.passes.synthesis import SolovayKitaev
from qiskit.quantum_info import Operator
from qiskit.synthesis.discrete_basis.generate_basis_approximations import (
    generate_basic_approximations,
)

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

basic_approximations = generate_basic_approximations(basis_gates=[HGate(), SGate(), TGate(), ZGate(), TdgGate(), SdgGate(), XGate(), YGate()], depth=10)
skd = SolovayKitaev(recursion_degree=2, basic_approximations = basic_approximations)
discretized = skd(qc)

print("\nDiscretized circuit:")
print(discretized.draw())

print("Error:", np.linalg.norm(Operator(qc).data - Operator(discretized).data))

sys.exit(0)
