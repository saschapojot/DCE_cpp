//
// Created by polya on 2/26/24.
//

#ifndef DCE_CPP_TIMEEVOLUTION_HPP
#define DCE_CPP_TIMEEVOLUTION_HPP
#include <string>
#include <iostream>
#include <memory>
#include <cstdio>
#include <array>
#include <regex>
#include <cmath>
#include <complex>
#include <vector>
#include <Eigen/Sparse>
#include <unsupported/Eigen/KroneckerProduct>

using namespace std::complex_literals;
const auto PI=std::numbers::pi;
class DCE_Evolution {
public:
    /// This constructor initializes all parameters
    /// @param group group number of parameters
    /// @param row row number of parameters
    DCE_Evolution(const int &group, const int &row) {
        //
        this->groupNum = group;
        this->rowNum = row;

        this->parseCSV(group, row);
        for (int n1 =0;n1<N1;n1++){
            this->x1ValsAll.push_back(-L1+dx1*n1);
        }
        for (int n2=0;n2<N2;n2++){
            this->x2ValsAll.push_back(-L2+dx2*n2);
        }

//        std::cout<<"dt="<<dt<<std::endl;
//        std::cout<<"dx1="<<dx1<<std::endl;
//        std::cout<<"dx2="<<dx2<<std::endl;
//        std::cout<<"M="<<M<<std::endl;
//        std::cout<<"x1vec has size "<<x1ValsAll.size()<<std::endl;
//        std::cout<<"x2vec has size "<<x2ValsAll.size()<<std::endl;


            

    }

public:
    int jH1 = -1;
    int jH2 = -1;
    double g0 = 0;
    double omegam = 0;
    double omegap=0;
    double omegac = 0;
    double er = 0;
    double thetaCoef = 0;
    int groupNum = -1;
    int rowNum = -1;
    double theta=0;
    double lmd=0;
    double Deltam=0;

     const int N1=2;//500;
     const int N2=2;//2048;
    double L1=5;
    double L2=40;

    double dx1=2*L1/(static_cast<double>(N1));
    double dx2=2*L2/(static_cast<double >(N2));
    double dtEst=0.002;
    double tFlushStart=0;
    double tFlushStop=5;
    int flushNum=1;
    double tTotPerFlush=tFlushStop-tFlushStart;
    int M=static_cast<int>(std::ceil(tTotPerFlush/dtEst));
    double dt=tTotPerFlush/(static_cast<double>(M) );

    std::vector<double> x1ValsAll;
    std::vector<double> x2ValsAll;

    Eigen::SparseMatrix<std::complex<double>> IN1N2;
//    Eigen::SparseMatrix<std::complex<double>>IN2;
//    Eigen::SparseMatrix<std::complex<double>> D2;
    Eigen::SparseMatrix<std::complex<double>> H0;
    Eigen::SparseMatrix<std::complex<double>> H2;
//    Eigen::SparseMatrix<std::complex<double>> S2;
    Eigen::SparseMatrix<std::complex<double>> H3;
    Eigen::SparseMatrix<std::complex<double>> H6;


    ///
    /// @param group group number
    /// @param row row number
    ///parse csv with group number group
    void parseCSV(const int &group, const int &row);

    ///
    /// @param cmd python execution string
    /// @return signal from the python
    static std::string execPython(const char *cmd);

    ///This function initializes the sparse matrices
    void populatedMatrices();

};


#endif //DCE_CPP_TIMEEVOLUTION_HPP
