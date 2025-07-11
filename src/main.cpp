#include <stdio.h>
#include <armadillo>
#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <cstdlib>
#include <unistd.h>
#include <omp.h>
#include <boost/throw_exception.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/complex.hpp>
#include <boost/program_options.hpp>
#include "LearningAlgorithms/GradientDescentPenalty.h"
#include "LearningAlgorithms/DeltaLearningPenalty.h"
#include "LearningAlgorithms/GradientDescentLagrange.h"
#include "LearningAlgorithms/NewtonsMethod.h"
#include "Functions/OrthogonalProcrustes.h"
#include "Functions/Orthogonalization.h"
#include "LearningRates/Constant.h"
#include "LearningRates/Wolfe.h"
#include "LearningRates/Adam.h"
#include "Functions/Noise.h"
#include "Globals.h"

using namespace std;
using namespace arma;
namespace po = boost::program_options;
using namespace std::chrono;

arma::cx_mat BuildMatrix(const vector<vector<cx_double>> &a);
vector<vector<cx_double>> ExportMatrix(const cx_mat &input);
arma::cx_mat LoadMatrix(const string &name);

bool sequential_mode = false;
bool matrix_analysis = false;
bool each_mat = false;
int n;
int num_x;
string output_name = "test.txt";
string matrix_type = "";
string method = "";
string decomposition = "";
int user_def_cond_x = 100;
int condx_upper = 100;
int condx_lower = 0;
int num_tests;
double alpha;
bool add_noise = false;
double noise_mag;
bool DEBUG_FLAG = false;

arma::cx_mat Utarget;

void print_mat(arma::cx_mat my_matrix)
{

	uint cols = my_matrix.n_cols;
	uint rows = my_matrix.n_rows;

	std::cout << "--------\n";
	for (uint rX = 0; rX < rows; rX++)
	{
		std::cout << " " << rX << ": ";
		for (uint cX = 0; cX < cols; cX++)
		{
			std::cout << my_matrix(rX, cX) << " ";
		}
		std::cout << "\n";
	}
	std::cout << "--------\n";
}

int main(int argc, const char *argv[])
{

	// Parse command-line options
	po::options_description desc("Allowed options");
	desc.add_options()("help", "produce help message")("sequential", "learn in sequential mode")("matrix_analysis", "print SVD and eig info about U and X")
		// TODO: Add matrix possibilities
		("matrix_type", po::value<string>(&matrix_type)->required(), "Matrix example")("n", po::value<int>(&n)->required(), "number of rows/cols of U")("method", po::value<string>(&method)->required(), "Learning method")("num_x", po::value<int>(&num_x), "Number of training examples, defaults to n")("alpha", po::value<double>(&alpha)->required(), "Learning rate")("each_mat", "Log all info associated with a test, not just average.")("cond_x", po::value<int>(&user_def_cond_x), "Condition number bound of X, default = 1000")("condx_upper", po::value<int>(&condx_upper), "Condition number bound of X, default = 100")("condx_lower", po::value<int>(&condx_lower), "Condition number bound of X, default = 0")("output_name", po::value<string>(&output_name), "name of output file, default = test.txt, outputs one test per line in format: iterations, condition of x")("num_tests", po::value<int>(&num_tests)->required(), "number of tests")("decomposition", po::value<string>(&decomposition), "Decompose the matrix to a circuit of unitary gates")("add_noise", "Add Gaussian Noise to Y vector")("noise_mag", po::value<double>(&noise_mag), "Specify magnitude of the noise")("debug", "Set debug statements");

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);

	if (vm.count("help"))
	{
		std::cout << desc << "\n";
		return 1;
	}
	try
	{
		po::notify(vm);
	}
	catch (const po::error &ex)
	{
		std::cerr << ex.what() << '\n';
		return 1;
	}

	std::cout << "       Welcome to        \n";
	std::cout << "   ____        __       \n";
	std::cout << "  / __ \\__  __/ /   ___ \n";
	std::cout << " / / / / / / / /   / _ \\\n";
	std::cout << "/ /_/ / /_/ / /___/  __/\n";
	std::cout << "\\___\\_\\__,_/_____/\\___/ \n";
	std::cout << "       Version 1.0\n";

	sequential_mode = vm.count("sequential");
	matrix_analysis = vm.count("matrix_analysis");
	each_mat = vm.count("each_mat");
	add_noise = vm.count("add_noise");
	DEBUG_FLAG = vm.count("debug");

	arma::cx_mat X, Y, U0, U;
	// arma::cx_mat X, Y, U0; // only to test NM, change back to the above declaration when finished
	// TODO: Add the other ventura matrices
	//  Do underconstrained test...
	if (matrix_type == "ex_const")
	{
		X = LoadMatrix("toronto_ventura_qft_exactly_constrained_X");
		Y = LoadMatrix("toronto_ventura_qft_exactly_constrained_Y");
		U0 = LoadMatrix("toronto_ventura_qft_exactly_constrained_U0");
		Utarget = U0;

	}

	std::cout << "Matrix size, n = " << n << endl;
	std::cout << "Number of tests = " << num_tests << endl;
	std::cout << "Printing data in " << output_name << endl;

	// If the user didn't specify a number of training examples, set size of X to n
	if (num_x == 0)
	{
		num_x = n;
	}

	// Generate a random nxn unitary U
	if (matrix_type == "rand")
	{
		// arma_rng::set_seed_random();
		U0 = randn<cx_mat>(n, n) + i1 * randn<cx_mat>(n, n);
		for (int k = 0; k < n; ++k)
			U0.col(k) /= norm(U0.col(k), 2);
		arma::cx_mat Q, R;
		if (!qr(Q, R, U0))
		{
			std::cout << "QR FAIL!!!" << endl;
		}
		R = diagmat(diagmat(R) / arma::abs(diagmat(R)));
		U0 = Q * R;
		Utarget = U0;
	}

	// Generates random x
	//  generate 10,000 random training examples
	/*X=randn<cx_mat>(n,n)+i1*randn<cx_mat>(n,n);
	for(int k=0;k<n;k++){
		X.col(k)/=norm(X.col(k),2);
	}*/
	// Y=U0*X;

	// Single-qubit Hadamard gate
	if (matrix_type == "Had")
	{
		U0 = BuildMatrix({{1.0, 1.0}, {1.0, -1.0}});
		U0 *= 1.0 / sqrt(2.0);
		Utarget = U0;
	}
	// QFT
	if (matrix_type == "QFT")
	{
		cx_double omega = exp(2.0 * datum::pi * cx_double(0., 1.) / (double)n);

		U0 = zeros<cx_mat>(n, n);
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < n; ++j)
				U0(i, j) = pow(omega, (double)i * j);
		U0 /= sqrt((double)n);
		Utarget = U0;
	}

	// TODO: Each matrix needs own X?
	//  QUESTION: Implemented correctly?
	if (matrix_type == "grover")
	{
		// X=ones<cx_mat>(n,n);
		// for(int k=0;k<n;k++)
		//     X(k,k)*=-1;
		// X/=sqrt(n);

		U0 = zeros<cx_mat>(n, n);
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < n; ++j)
				U0(i, j) = i == j ? 2.0 / (double)n - 1.0 : 2.0 / (double)n;
		Utarget = U0;
	}

	vector<int> counts(num_tests, 0);
    vector<int> stalls(num_tests, 0);
	vector<unsigned long long int> wall_times(num_tests, 0);
	double totalFNormYUX = 0;
	double totalFNormUHUI = 0;

	double froDiffBetweenNoiseAndClean = 0;

	FILE *outputptr;
	outputptr = fopen(output_name.c_str(), "w");

	FILE *logptr;
	string log_name = "log.csv";
	logptr = fopen(log_name.c_str(), "a");

    FILE *logptr_indv;
	string log_name_indv = method + "_" + matrix_type + "_" + std::to_string(n) + ".csv";
    while(log_name_indv.find("*") != string::npos) log_name_indv.replace(log_name_indv.find("*"), 1, "-");
    if (matrix_type != "rand") // only using this for the first bar chart, and so far ended up not using it even there
	    logptr_indv = fopen(log_name_indv.c_str(), "a");

	
	// Learning rate log
	// std::string alphaLogName = "AdamAlphaLog_" + method + "_" + std::to_string(n) + ".csv";
	// std::ofstream alphaLog(alphaLogName, std::ios::out | std::ios::trunc);

	// Check if file is empty
	fseek(logptr, 0, SEEK_END);
	if (ftell(logptr) == 0)
	{
		// Add headers
		fprintf(logptr, "method,matrix_type,alpha,sequential_mode,condx_lower,condx_upper,condx,n,num_x,num_tests,mat_num,fro_Y-UX,fro_UHU-I,avg_wall_time,avg_iters,n_stalls,noise_mag,noise_nonoise_frodiff");
	}

	for (int test = 0; test < num_tests; ++test)
	{
		if (matrix_analysis)
		{
			cx_mat svd_U;
			cx_mat svd_V;
			vec svd_s;
			cx_vec eigval_U;
			cx_mat eigvec_U;

			svd(svd_U, svd_s, svd_V, U0);
			eig_gen(eigval_U, eigvec_U, U0);
			std::cout << "singular values of U0: " << svd_s << endl;
			std::cout << "eigenvalues of U0: " << eigval_U << endl;
			std::cout << "cond(U0) = " << cond(U0) << endl;
			std::cout << "rank(U0) = " << arma::rank(U0) << endl;
		}

		// Generate a random X and enforce that cond(X) <= 1000
		// Skip random X generation for ex_const case since X is loaded from file
		int condx = 1000000;  // Declare condx outside the if block so it's accessible later
		if (matrix_type != "ex_const")
		{
			// while(condx>user_def_cond_x){
			// TODO: Is there an edge case here?
			while (!(condx >= condx_lower && condx <= condx_upper))
			{ // condx>condx_upper && condx > condx_lower){
				// arma_rng::set_seed_random();
				X = randn<cx_mat>(n, num_x) + i1 * randn<cx_mat>(n, num_x);
				for (int k = 0; k < num_x; ++k)
					X.col(k) /= norm(X.col(k), 2);
				condx = cond(X);
			}
		}
		else
		{
			// For ex_const case, get the condition number of the loaded X matrix
			condx = cond(X);
		}
		// std::cout << X << endl;
		//  } //else {//if (matrix_type != "grover") {
		//     arma_rng::set_seed_random();
		//     X=randn<cx_mat>(n,n)+i1*randn<cx_mat>(n,n);
		//     for(int k=0;k<n;++k)
		//         X.col(k)/=norm(X.col(k),2);
		// }
		// build sigma
		// create an identity matrix
		cx_mat sigma = eye<cx_mat>(n, num_x);
		sigma(n - 1, num_x - 1) = randu<cx_double>();
		sigma(0, 0) *= condx;
		// generate an array of random numbers inbetween the value of the first element of sigma and the last
		vec sigma_rand = randu<vec>(n - 2);
		// sort sigma_rand in descending order
		sigma_rand = sort(sigma_rand, "descend");
		// multiply the rest of the diagnal of sigma by corresponding entries of sigma_rand
		for (int k = 1; k < n - 1; ++k)
		{
			sigma(k, k) *= sigma_rand(k - 1);
		}
		// std::cout<< "sigma: \n" << sigma <<endl;
		// get a random orthogonal matrix U where the norm of each column is 1
		cx_mat Ux = randn<cx_mat>(n, num_x) + i1 * randn<cx_mat>(n, num_x);
		for (int k = 0; k < num_x; ++k)
			Ux.col(k) /= norm(Ux.col(k), 2);
		// get a VT which is the transpose of a nxn matrix containing the orthonormal eigenvectors of sigma
		cx_mat VT = randn<cx_mat>(n, n) + i1 * randn<cx_mat>(n, n);
		for (int k = 0; k < n; ++k)
			VT.col(k) /= norm(VT.col(k), 2);

		// get a random matrix VT
		Y = Utarget * X;
		arma::cx_mat Y_noise(Y);

		if (add_noise)
		{
			if (noise_mag == 0)
			{
				noise_mag = 0.000001;
			}
			GaussianNoise noiseMaker = GaussianNoise();
			// print_mat(Y);
			Y_noise = noiseMaker.addNoise(Y, noise_mag);
			// print_mat(Y);
		}

		if (matrix_analysis)
		{
			cx_mat Xsvd_U;
			cx_mat Xsvd_V;
			vec Xsvd_s;
			cx_vec Xeigval_U;
			cx_mat Xeigvec_U;

			svd(Xsvd_U, Xsvd_s, Xsvd_V, X);
			eig_gen(Xeigval_U, Xeigvec_U, X);
			std::cout << "singular values of X: " << Xsvd_s << endl;
			std::cout << "eigenvalues of X: " << Xeigval_U << endl;
			std::cout << "cond(X) = " << cond(X) << endl;
			std::cout << "rank(X) = " << arma::rank(X) << endl;
		}

		// forget U0 and input an initial guess
		U0 = eye<cx_mat>(n, n);

		Function *DescentObjective = new OrthogonalProcrustes();
		vector<Function *> DescentConstraints = {new Orthogonalization()};
		vector<Function *> NoDescentConstraints = {};
		LearningRate *DescentLearningRate = new Constant(alpha);

		// LearningRate *DescentLearningRate = new Wolfe(alphaLog);

		// LearningRate *DescentLearningRate = new Adam(alphaLog, alpha);
		// LearningRate *DescentLearningRate = new Adam(alpha);

		std::vector<double> DescentPenalties = {1.0};
		LearningAlgorithm *DescentLearningAlgorithm;

		if (method == "GDP")
		{
			DescentLearningAlgorithm = new GradientDescentPenalty(DescentObjective, NoDescentConstraints, DescentLearningRate, DescentPenalties);
		}
		else if (method == "GDP*")
		{
			DescentLearningAlgorithm = new GradientDescentPenalty(DescentObjective, DescentConstraints, DescentLearningRate, DescentPenalties);
		}
		else if (method == "GDLM")
		{
			DescentLearningAlgorithm = new GradientDescentLagrange(DescentObjective, NoDescentConstraints, DescentLearningRate);
		}
		else if (method == "GDLM*")
		{
			DescentLearningAlgorithm = new GradientDescentLagrange(DescentObjective, DescentConstraints, DescentLearningRate);
		}
		else if (method == "DLR")
		{
			DescentLearningAlgorithm = new DeltaLearningPenalty(DescentObjective, NoDescentConstraints, DescentLearningRate, DescentPenalties);
		}
		else if (method == "DLR*")
		{
			DescentLearningAlgorithm = new DeltaLearningPenalty(DescentObjective, DescentConstraints, DescentLearningRate, DescentPenalties);
		}
		else if (method == "NM")
		{
			// TODO: Deal with this warning
			//  \warning: Newton's Method is designed solely for unconstrained Procrustes Problem. Will need to work out Hessian of Unitarization constraint
			// Add in a check for unconstrained Procrustes Problem
			DescentLearningAlgorithm = new NewtonsMethod(DescentObjective, DescentConstraints, DescentLearningRate);
		}
		else
		{
			std::cout << method << " is not a valid method" << std::endl;
			std::cout << "Valid methods are: GDP, GDP*, GDLM, GDLM*, DLR, DLR*, NM where * means a unitary constraint is applied" << std::endl;
			return 1;
		}

		cx_mat U = zeros<cx_mat>(n, n);
		cx_mat U_nonoise = zeros<cx_mat>(n, n); // only relavent if add_noise is true
		if (sequential_mode)
		{
			for (int i = 0; i < n; ++i)
			{
				high_resolution_clock::time_point t1 = high_resolution_clock::now();
				pair<cx_mat, int> r = (*DescentLearningAlgorithm).Learn(U0.row(i), X, Y.row(i));
				U.row(i) = r.first;
				if(r.second == -1)
					stalls[test] = 1;
				else
					// counts[test] += r.second;
					// change to max to account for the fact that some rows may take longer to converge than others
					counts[test] = max(counts[test], r.second);
				high_resolution_clock::time_point t2 = high_resolution_clock::now();
				wall_times[test] += duration_cast<microseconds>(t2 - t1).count();
			}
		}
		else if (add_noise)
		{
			high_resolution_clock::time_point t1 = high_resolution_clock::now();
			pair<cx_mat, int> r = (*DescentLearningAlgorithm).Learn(U0, X, Y_noise);
			U = r.first;
			if(r.second == -1)
				stalls[test] = 1;
			else
				counts[test] = r.second;
			high_resolution_clock::time_point t2 = high_resolution_clock::now();
			wall_times[test] = duration_cast<microseconds>(t2 - t1).count();

			r = (*DescentLearningAlgorithm).Learn(U0, X, Y);
			U_nonoise = r.first;
		}
		else
		{
			high_resolution_clock::time_point t1 = high_resolution_clock::now();
			pair<cx_mat, int> r = (*DescentLearningAlgorithm).Learn(U0, X, Y);
			U = r.first;
			if(r.second == -1)
				stalls[test] = 1;
			else
				counts[test] = r.second;
			high_resolution_clock::time_point t2 = high_resolution_clock::now();
			wall_times[test] = duration_cast<microseconds>(t2 - t1).count();
		}
		std::cout << "Finished! " << test << endl;
		bool print_U = false;
		if (print_U)
			std::cout << "Learned U:" << endl
					  << U << endl;

		// QUESTION: How are these frobenius norms different?
		/*for(int m=0;m<n;m++)
			cout<<"Fro norm of row "<<m<<" = "<<norm(Y.row(m)-U.row(m)*X,"fro")<<endl;
		cout<<"Fro norm of U-GOODU = "<<norm(U-GOODU,"fro")<<endl;
		cout<<"Fro norm of Y-UX = "<<norm(Y-U*X,"fro")<<endl;
		for(int m=0;m<n;m++){
			cx_mat u=U.row(m);
			cx_mat uX=u*X;
			cx_mat p=U*X;
			cx_mat pr=p.row(m);
			cout<<"Y-UX row "<<m<<" = "<<uX-pr<<endl;
		}

		cout<<"Y-UX= "<<U*X<<endl;*/
		std::cout << "Frobenius norm of Y-UX = " << norm(Y - U * X, "fro") << endl;
		std::cout << "Frobenius norm of U^H U - I = " << norm(U.t() * U - eye<cx_mat>(n, n), "fro") << endl;

		totalFNormYUX += norm(Y - U * X, "fro");
		totalFNormUHUI += norm(U.t() * U - eye<cx_mat>(n, n), "fro");

		if (add_noise)
		{
			// froDiffBetweenNoiseAndClean += norm(Y - U * X, "fro") + norm(U.t() * U - eye<cx_mat>(n, n), "fro") - norm(Y - U_nonoise * X, "fro") + norm(U_nonoise.t() * U_nonoise - eye<cx_mat>(n, n), "fro");
			froDiffBetweenNoiseAndClean += norm(U - U_nonoise, "fro");
		}

		double Xcond = cond(X);
		std::cout << "Cond(X) = " << Xcond << endl;
		std::cout << "Iterations to convergence = " << counts[test] << endl;
		std::cout << endl;
		fprintf(outputptr, "%d,%.20f\n", counts[test], Xcond);

		if (each_mat)
		{
			fprintf(logptr, "\n%s,", method.c_str());
			fprintf(logptr, "%s,", matrix_type.c_str());
			fprintf(logptr, "%f,", alpha);
			fprintf(logptr, "%d,", sequential_mode);
			fprintf(logptr, "%d,", condx_lower);
			fprintf(logptr, "%d,", condx_upper);
			fprintf(logptr, "%f,", Xcond);
			// fprintf(logptr,"%d,", user_def_cond_x);
			fprintf(logptr, "%d,", n);
			fprintf(logptr, "%d,", num_x);
			fprintf(logptr, "%d,", num_tests);
			// TODO: Can I outprint wall time for an individual matrix?
			// represents matrix number
			fprintf(logptr, "%d,", test);
			fprintf(logptr, "%.10f,", norm(Y - U * X, "fro"));
			fprintf(logptr, "%.10f,", norm(U.t() * U - eye<cx_mat>(n, n), "fro"));
			// set to -1 for wtsum. TODO: Check that that makes sense
			fprintf(logptr, "%-d,", -1);
			fprintf(logptr, "%d,", counts[test]);
			fprintf(logptr, "%d,", stalls[test]); // this is just 0 or 1, 1 if this test stalled
			fprintf(logptr, "%f,", noise_mag);
			fprintf(logptr, "%f",  0.0);
		}

        if (matrix_type != "rand")
        {
            if (stalls[test] == 0)
                fprintf(logptr_indv, "%d\n", counts[test]);
        }

		// std::cout<<"partial unitarity = "<<norm(U.cols(0,n-9).t()*U.cols(0,n-9)-eye<cx_mat>(n-8,n-8),"fro")<<endl;

		// QUESTION: Is unitarization a constraint or done as a post-processing step?
		//  Unitarization step
		/*Function* OrthoObjective=new Orthogonalization();
		vector<Function*> OrthoConstraints={};
		LearningRate* OrthoLearningRate=new Constant(0.0001);
		std::vector<double> OrthoPenalties={};
		GradientDescentPenalty OrthoLearningAlgorithm(OrthoObjective,OrthoConstraints,OrthoLearningRate,OrthoPenalties);
		cx_mat U2=U.cols(n-8,n-1);
		U2=OrthoLearningAlgorithm.Learn(U2,X,Y).first;
		U.cols(n-8,n-1)=U2;
		cout << "Finished Orthogonalization via Gradient Descent" << endl;
		cout<<"Fro norm of Y-UX = "<<norm(Y-U*X,"fro")<<endl;
		cout<<"Fro norm of U^H U - I = "<<norm(U.t()*U-eye<cx_mat>(n,n),"fro")<<endl<<endl;*/
	}

	std::fclose(outputptr);

	int sum = 0;
	for (int s : counts)
		sum += s;

    int n_stalls = 0;
	for (int s : stalls)
		n_stalls += s;

	double avg_iters = num_tests == n_stalls ? 0 : (double)sum / ((double)num_tests - (double)n_stalls);
	std::cout << "Average number of iterations (when convergence didn't stall): " << avg_iters << endl;

	unsigned long long int wtsum = 0;
    if (num_tests != n_stalls){
        for (unsigned long long int s : wall_times)
            wtsum += s;
        std::cout << "Average wall time across " << (num_tests-n_stalls) << " trials (microseconds): " << wtsum / (unsigned long long int)(num_tests - n_stalls) << endl;
        std::cout << (double)wtsum / (double)(num_tests - n_stalls) << std::endl;
        std::cout << "Average wall time across " << (num_tests - n_stalls) << " trials (seconds): " << wtsum / (unsigned long long int)(num_tests - n_stalls) / (unsigned long long int)1000000 << endl;
        std::cout << (double)wtsum / (double)(num_tests - n_stalls) / (double)1000000 << std::endl;
        std::cout << "Average Frobenius norm of Y-UX = " << (double)totalFNormYUX / (double)(num_tests-n_stalls) << std::endl;
        std::cout << "Average Frobenius norm of U^H U - I = " << (double)totalFNormUHUI / (double)(num_tests-n_stalls) << std::endl;
        std::cout << "wtsum = " << wtsum << endl;
        std::cout << endl;
    }

	// Outprint infomation to log file so that python script can plot data
	// method,matrix_type,alpha,sequential_mode,condx_lower,condx_upper,condx,n,num_x,num_tests,mat_num,fro_Y-UX,fro_UHU-I,avg_wall_time,avg_iters
	fprintf(logptr, "\n%s,", method.c_str());
	fprintf(logptr, "%s,", matrix_type.c_str());
	fprintf(logptr, "%f,", alpha);
	fprintf(logptr, "%d,", sequential_mode);
	fprintf(logptr, "%d,", condx_lower);
	fprintf(logptr, "%d,", condx_upper);
	fprintf(logptr, "%-d,", -1);
	// fprintf(logptr,"%d,", user_def_cond_x);
	fprintf(logptr, "%d,", n);
	fprintf(logptr, "%d,", num_x);
	fprintf(logptr, "%d,", num_tests);
	// represents matrix number (python script will check if this is equal to the num_tests and then
	// will know this entry is the average of all the trials)
	fprintf(logptr, "%d,", num_tests);
	fprintf(logptr, "%.10f,", num_tests == n_stalls ? 0. : (double)totalFNormYUX / (double)(num_tests-n_stalls));
	fprintf(logptr, "%.10f,", num_tests == n_stalls ? 0. :(double)totalFNormUHUI / (double)(num_tests-n_stalls));
	fprintf(logptr, "%f,", num_tests == n_stalls ? 0. :(double)wtsum / (double)(num_tests-n_stalls) / (double)1000000);
	fprintf(logptr, "%f,", avg_iters);
	fprintf(logptr, "%d,", n_stalls);
	fprintf(logptr, "%f,", noise_mag);
	fprintf(logptr, "%f",  froDiffBetweenNoiseAndClean);

	std::fclose(logptr);
	if (matrix_type != "rand")
		std::fclose(logptr_indv);
	// BEN: A test of Newtons Method
	/*Function* MyObjective=new OrthogonalProcrustes();
	vector<Function*> MyConstraints={};
	LearningRate* MyLearningRate=new Wolfe();//Constant(0.01);
	NewtonsMethod MyAlgorithm(MyObjective,MyConstraints,MyLearningRate);
	cx_mat U=MyAlgorithm.Learn(U0,X,Y);
	cout << "Finished Newton's Method" << endl;
	cout << "Learned U = " << endl;
	cout << U << endl;
	cout << endl;*/

	// Unitarization step
	//  cx_mat UOrth=zeros<cx_mat>(n,n);
	// Function* OrthoObjective=new Orthogonalization();
	// vector<Function*> OrthoConstraints={};
	// LearningRate* OrthoLearningRate=new Constant(0.01);
	// std::vector<double> OrthoPenalties={};
	// GradientDescentPenalty OrthoLearningAlgorithm(OrthoObjective,OrthoConstraints,OrthoLearningRate,OrthoPenalties);
	// cx_mat UOrth = OrthoLearningAlgorithm.Learn(UOrth,X,Y);
	// cout << "Finished Orthogonalization via Gradient Descent" << endl;
	// cout<<"Fro norm of Y-UX = "<<norm(Y-UOrth*X,"fro")<<endl;
	// cout<<"Fro norm of U^H U - I = "<< norm(UOrth.t()*UOrth-eye<cx_mat>(n,n),"fro") <<endl<<endl;
	// cout << "Learned U = " << endl;
	// cout << U << endl;

	// Solovay-Kitaev factorization, user should have python 3 installed
	if (decomposition == "qiskit")
	{
		const char *python = "/usr/bin/python3";
		const char *script = "./Decomposition/qiskit_decomposition.py";
		char *args[4] = {0};

		// Convert the matrix to a string representation
		std::stringstream ss;
		Utarget.save(ss, arma::raw_ascii);
		std::string str = ss.str();
		// Convert the string to a char* buffer
		const char *buffer = str.c_str();

		std::string n_str = std::to_string(n);
		const char *matrix_size = n_str.c_str();

		args[0] = (char *)python;
		args[1] = (char *)script;
		args[2] = (char *)buffer;
		args[3] = (char *)matrix_size;

		execvp(python, args);
		// If execvp returns, it means there was an error
		std::cerr << "Error executing factorization" << std::endl;
		return 1;
	}

	if (decomposition == "SQUANDER")
	{
		const char *python = "/usr/bin/python3";
		const char *script = "./Decomposition/squander.py";
		char *args[4] = {0};

		// Convert the matrix to a string representation
		std::stringstream ss;
		Utarget.save(ss, arma::raw_ascii);
		std::string str = ss.str();
		// Convert the string to a char* buffer
		const char *buffer = str.c_str();

		std::string n_str = std::to_string(n);
		const char *matrix_size = n_str.c_str();

		args[0] = (char *)python;
		args[1] = (char *)script;
		args[2] = (char *)buffer;
		args[3] = (char *)matrix_size;

		execvp(python, args);
		// If execvp returns, it means there was an error
		std::cerr << "Error executing factorization" << std::endl;
		return 1;
	}

	if (decomposition == "OPENQL")
	{
		const char *python = "/usr/bin/python3";
		const char *script = "./Decomposition/openQLs.py";
		char *args[4] = {0};

		// Convert the matrix to a string representation
		std::stringstream ss;
		Utarget.save(ss, arma::raw_ascii);
		std::string str = ss.str();
		// Convert the string to a char* buffer
		const char *buffer = str.c_str();

		std::string n_str = std::to_string(n);
		const char *matrix_size = n_str.c_str();

		args[0] = (char *)python;
		args[1] = (char *)script;
		args[2] = (char *)buffer;
		args[3] = (char *)matrix_size;

		execvp(python, args);
		// If execvp returns, it means there was an error
		std::cerr << "Error executing factorization" << std::endl;
		return 1;
	}

	return 0;
}

cx_mat BuildMatrix(const vector<vector<cx_double>> &a)
{
	cx_mat A(a.size(), a[0].size());
	for (unsigned int i = 0; i < a.size(); ++i)
		for (unsigned int j = 0; j < a[i].size(); ++j)
			A(i, j) = a[i][j];
	return A;
}

vector<vector<cx_double>> ExportMatrix(const cx_mat &input)
{
	vector<vector<cx_double>> A(input.n_rows);
	for (unsigned int i = 0; i < input.n_rows; ++i)
	{
		A[i].resize(input.n_cols);
		for (unsigned int j = 0; j < input.n_cols; ++j)
			A[i][j] = input(i, j);
	}
	return A;
}

cx_mat LoadMatrix(const string &name)
{
	vector<vector<cx_double>> input;
	std::ifstream ifs("../data/" + name + ".txt");
	boost::archive::text_iarchive ia(ifs);
	ia &input;
	return BuildMatrix(input);
}