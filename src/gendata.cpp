#include <stdio.h>
#include <armadillo>
#include <iostream>
#include <fstream>
#include <string>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/complex.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

using namespace std;
using namespace arma;
namespace po=boost::program_options;

cx_mat BuildMatrix(const vector<vector<cx_double> >& a);
void WriteMatrix(const cx_mat& input,const string& name);
vector<vector<cx_double> > ExportMatrix(const cx_mat& input);

bool overwrite=false;

int main(int argc,const char* argv[])
{
    // Parse command-line options
    po::options_description desc("Allowed options");
desc.add_options()
    ("help", "produce help message")
    ("overwrite", "overwrite existing data files");

    po::variables_map vm;
    po::store(po::parse_command_line(argc,argv,desc),vm);
    po::notify(vm);
    
    if(vm.count("help")){
        cout<<desc<<"\n";
        return 1;
    }

    if(vm.count("overwrite"))
        overwrite=true;

    // ADD YOUR DATA HERE!
    cx_mat X; // use only one matrix

    // Ventura 2000: learn Hadamard from 2 examples
    // appears this ex converges in 1 or 2 steps of SGD for Ventura...
    X=BuildMatrix({{2.0/sqrt(5.0),-2.0/sqrt(20.0)},{1.0/sqrt(5.0),4.0/sqrt(20.0)}});
    WriteMatrix(X,"ventura_2000_X");
    X=BuildMatrix({{3.0/sqrt(10.0),2/sqrt(40.0)},{1.0/sqrt(10.0),-6.0/sqrt(40.0)}});
    WriteMatrix(X,"ventura_2000_Y");
    X=BuildMatrix({{0,1},{-1,0}});
    WriteMatrix(X,"ventura_2000_U0");

    // Toronto/Ventura 2006: Quantum Fourier Transform; exactly constrained
    X=BuildMatrix({{1.,sqrt(2.),sqrt(2.),2.},{1.,0.,0.,0.},{1.,sqrt(2.),0.,0.},{1.,0.,sqrt(2.),0.}});
    X*=0.5;
    WriteMatrix(X,"toronto_ventura_qft_exactly_constrained_X");
    X=BuildMatrix({{2.,sqrt(2.),sqrt(2.),1.},{0.,0.,exp(7.*datum::pi*cx_double(0.,1.)/4.),1.},
                   {0.,sqrt(2.),0.,1.},{0.,0.,exp(datum::pi*cx_double(0.,1.)/4.),1.}});
    X*=0.5;
    WriteMatrix(X,"toronto_ventura_qft_exactly_constrained_Y");
    X=eye<cx_mat>(X.n_rows,X.n_rows);
    WriteMatrix(X,"toronto_ventura_qft_exactly_constrained_U0");

    // Toronto/Ventura 2006: underconstrained operator learning
    X=BuildMatrix({{1,0},{0,1},{0,0},{0,0}});
    WriteMatrix(X,"toronto_ventura_underconstrained_X");
    X=BuildMatrix({{sqrt(2),0},{sqrt(2),0},{0,sqrt(2)},{0,sqrt(2)}});
    X*=0.5;
    WriteMatrix(X,"toronto_ventura_underconstrained_Y");
    X=eye<cx_mat>(4,4);
    WriteMatrix(X,"toronto_ventura_underconstrained_U0");

    // Toronto/Ventura 2006: Grover's iterate learning (non-unitary)
    X=ones<cx_mat>(8,8);
    for(int i=0;i<8;++i)
        X(i,i) *= -1.;
    X*=sqrt(2.)/4.;
    WriteMatrix(X,"toronto_ventura_grover_X");
    X=eye<cx_mat>(X.n_rows,X.n_rows);
    WriteMatrix(X,"toronto_ventura_grover_Y");
    WriteMatrix(X,"toronto_ventura_grover_U0");

    return 0;
}

void WriteMatrix(const cx_mat& input,const string& name){
    if(boost::filesystem::exists("../../../data"+name+".txt")&&!overwrite)return;
    vector<vector<cx_double> > data=ExportMatrix(input);
    std::ofstream ofs("../../../data/"+name+".txt");
    boost::archive::text_oarchive oa(ofs);
    oa & data;
}

cx_mat BuildMatrix(const vector<vector<cx_double> >& a) {
    cx_mat A(a.size(),a[0].size());
    for(unsigned int i=0;i<a.size();i++)
        for(unsigned int j=0;j<a[i].size();j++)
            A(i,j)=a[i][j];
    return A;
}

vector<vector<cx_double> > ExportMatrix(const cx_mat& input){
    vector<vector<cx_double> > A(input.n_rows);
    for(unsigned int i=0;i<input.n_rows;i++){
        A[i].resize(input.n_cols);
        for(unsigned int j=0;j<input.n_cols;j++)
            A[i][j]=input(i,j);
    }
    return A;
}
