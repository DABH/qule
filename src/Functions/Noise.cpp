#include "Noise.h"

cx_mat GaussianNoise::addNoise(const cx_mat& X, double magnitude = 1) const
{
	cx_mat noiseMatrix = magnitude*randn<cx_mat>(X.n_rows, X.n_cols); //  distr_param(10, 20)
	return X + noiseMatrix;
}