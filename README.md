# QuLe
Finding Quantum Operators Through Machine Learning


## About

QuLe (rhymes with "mule") is a program for generating quantum operators (in the form of unitary matrices) that (as closely as possible) map given input states to given output states of a quantum
computer.  QuLe is introduced and explained in the article:

## Setup

```
cd ./src
cmake .
make
```

## Usage

Sample Usage: `./qule --matrix_type rand --alpha 0.1 --method GDLM --n 6 --num_tests 5`

- Matrix type can be of type: ex_const, rand, Had (Hadamard Gate), QFT (Quantum Fourier Transform), grover
- Method can be of type: GDP (Gradient Descent Penalty), GDP*, GDLM (Gradient Descent Lagrange Multipler), DLR (Delta Learning Rule), DLR*, and New (Newtons Method). (*'s represent unitary constraint).
- alpha is the learning rate
- n is the size of the nxn matrix
- num_tests is the number of iterations to run

#### Optional Flags

- sequential: if set, learns in sequential mode instead of batch
- matrix_analysis: if set, prints SVD and eig info about U and X
- cond_x: set condition bound of x, default is 1000
- num_x: set the number of training examples, default is n
- each_mat: will log the data for every matrix in a test (not just the average)
- condx_upper: sets an upper bound on the condition number of randomly generated X
- condx_lower: sets a lower bound on the condition number of randomly generated X

## Examples


## Feedback


## License



