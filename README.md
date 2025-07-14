# QuLe

Finding Quantum Operators Through Machine Learning

## About

QuLe (rhymes with "mule") is a program for generating quantum operators (in the form of unitary matrices) that map given input states to given output states of a quantum computer (as closely as possible).
When combined with unitary factorization algorithms (as included in this repository and described in our paper), QuLe enables end-to-end synthesis of quantum algorithms: starting with a set of desired input/output pairs, the end result of our pipeline can be a quantum circuit diagram that closely approximates the desired behavior.

QuLe is introduced and explained in the research paper: (Coming Soon).
Please [cite](#citation) our paper if you find our work useful.

## Prerequsites and Setup

QuLe has been test on Windows and Linux.  Compilation on Mac should be possible, but will require (1) installing a suitable compiler that supports OpenMP (e.g., `brew install gcc`) and using that compiler for CMake (e.g., `export CC=gcc-15 && export CXX=g++-15`), and (2) ensuring paths for libraries like Armadillo are specified correctly.  A PR is welcomed to make Mac work smoothly out of the box.

QuLe depends on:
- [Armadillo](https://arma.sourceforge.net/download.html)
- [Boost](https://www.boost.org/) (Note, this can often be installed more easily via package managers like `apt` or `brew`)

Optionally, on Windows, you may build QuLe using [Intel MKL](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html).

One prerequisites are installed, you can build QuLe simply via:

```
cd ./src
cmake .
make -j
```

## Usage

Sample Usage: `./qule --matrix_type rand --alpha 0.1 --method GDLM --n 6 --num_tests 5`

- Matrix type can be of type: ex_const, rand, Had (Hadamard Gate), QFT (Quantum Fourier Transform), grover
- Method can be of type: GDP (Gradient Descent Penalty), GDP*, GDLM (Gradient Descent Lagrange Multipler), DLR (Delta Learning Rule), DLR*, and New (Newtons Method). (*'s represent unitary constraint).
- alpha is the learning rate
- n is the size of the nxn matrix
- num_tests is the number of iterations to run

#### Optional Flags

A few flags are listed here.  We recommend browsing the source code to see the most up-to-date list of flags and their meanings.

- sequential: if set, learns in sequential mode instead of batch
- matrix_analysis: if set, prints SVD and eig info about U and X
- cond_x: set condition bound of x, default is 1000
- num_x: set the number of training examples, default is n
- each_mat: will log the data for every matrix in a test (not just the average)
- condx_upper: sets an upper bound on the condition number of randomly generated X
- condx_lower: sets a lower bound on the condition number of randomly generated X

## Feedback

Please open an [issue](https://github.com/DABH/qule/issues) or [pull request](https://github.com/DABH/qule/pulls) to provide feedback on the project.

## License

QuLe is licensed under the University of Illinois / NCSA Open Source License.  See [LICENSE.md](LICENSE.md).

## Citation

If you use QuLe in your research or would otherwise like to reference it, please use the following citation:

Yuxin Huang, Benjamin Grossman-Ponemon, David Hyde.  "Automated Synthesis of Quantum Algorithms via Classical Numerical Techniques."  ACM Transactions on Quantum Computing.  In press (2025).

## Sponsorhsip

Please consider [sponsorship](https://github.com/sponsors/DABH) to support further development and maintenance of QuLe.
