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
#include <boost/filesystem.hpp>
#include <memory>
#include <msgpack.hpp>
#include <fstream>

namespace fs = boost::filesystem;
using namespace std::complex_literals;
const auto PI=std::numbers::pi;

using wvVec=Eigen::VectorXcd;

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

    }//end of constructor

    ///
    /// @param j index of time
    /// @return H1j
    Eigen::SparseMatrix<std::complex<double>> H1Val(const int & j);

    ///
    /// @param j  index of time
    /// @return H4j
    Eigen::SparseMatrix<std::complex<double>> H4Val(const int & j);

    ///
    /// @param j index of time
    /// @return H5j
    Eigen::SparseMatrix<std::complex<double>>H5Val(const int& j);

    ///initialize wavefunction
    void initPsi();


    template<class T>
    static void printVec(const std::vector <T> &vec) {
        for (int i = 0; i < vec.size() - 1; i++) {
            std::cout << vec[i] << ",";
        }
        std::cout << vec[vec.size() - 1] << std::endl;
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


     static const int N1=500;
     static const int N2=4096;
    double L1=5;
    double L2=80;


    double dx1=2*L1/(static_cast<double>(N1));
    double dx2=2*L2/(static_cast<double >(N2));
    double dtEst=0.002;
    double tFlushStart=0;
    double tFlushStop=0.01;
    int flushNum=3;
    std::vector<int> jIndsAll;//index for time steps
    double tTotPerFlush=0;
    int stepsPerFlush=0;//total time steps
    double dt=0;
    std::vector<int> timeIndsAll;//all indices of time

    std::string outDir;

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
    Eigen::SparseMatrix<std::complex<double>> H7;
    Eigen::SparseMatrix<std::complex<double>> H8;

    Eigen::SparseMatrix<std::complex<double>> HSumStatic;//sum of H0,H2,H3,H6,H7,H8
//    using wvVec=Eigen::Vector<std::complex<double>,N1*N2>;
    std::vector<wvVec> PsiAll;
    wvVec Psi0;


public:
    ///
    /// @param group group number
    /// @param row row number
    ///parse csv with group number group
    void parseCSV(const int &group, const int &row);

    ///
    /// @param cmd python execution string
    /// @return signal from the python
    static std::string execPython(const char *cmd);

    ///This function initializes the sparse matrices H0,H2,H3,H6,H7,H8
    void populatedMatrices();

    ///
    /// @param j index of time
    /// @return HDj
    Eigen::SparseMatrix<std::complex<double>> HDj(const int &j);

    ///create output directories
    void createOutDir();



    ///initialize time indices and dt
    void initTimeInds();

    ///
    /// @param j time step
    /// @param PsiCurr wavefunction before evolution
    /// @return wavefunction after evolution
    wvVec oneStepEvolution(const int& j, const wvVec& PsiCurr);

    //evolution and write to file by flush
    wvVec evolutionPerFlush(const int &fls, const wvVec& initVec);

    ///evolution
    void evolution();

    ///
    /// @param solutions solutions per flush
    /// @return converted to cpp data type
    std::vector<std::vector<std::complex<double>>> eigen2cppType(const std::vector<wvVec > & solutions);



};
namespace msgpack {
    MSGPACK_API_VERSION_NAMESPACE(MSGPACK_DEFAULT_API_NS) {
        namespace adaptor {

            template<>
            struct convert<std::complex<double>> {
                msgpack::object const& operator()(msgpack::object const& o, std::complex<double>& v) const {
                    if(o.type != msgpack::type::ARRAY || o.via.array.size != 2) {
                        throw msgpack::type_error();
                    }
                    v.real(o.via.array.ptr[0].as<double>());
                    v.imag(o.via.array.ptr[1].as<double>());
                    return o;
                }
            };

            template<>
            struct pack<std::complex<double>> {
                template <typename Stream>
                msgpack::packer<Stream>& operator()(msgpack::packer<Stream>& o, std::complex<double> const& v) const {
                    // The complex number is packed as an array of two doubles.
                    o.pack_array(2);
                    o.pack(v.real());
                    o.pack(v.imag());
                    return o;
                }
            };

        } // namespace adaptor
    } // MSGPACK_API_VERSION_NAMESPACE(MSGPACK_DEFAULT_API_NS)
} // namespace msgpack


#endif //DCE_CPP_TIMEEVOLUTION_HPP
