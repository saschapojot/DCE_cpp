//
// Created by polya on 3/12/24.
//

#ifndef DCE_CPP_COMBINESEGMENTS_HPP
#define DCE_CPP_COMBINESEGMENTS_HPP
#include <iostream>
#include <regex>
#include <cmath>
#include <complex>
#include <boost/filesystem.hpp>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <msgpack.hpp>
#include <Eigen/Sparse>
#include <boost/json.hpp>
#include <sstream>
#include <chrono>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/complex.hpp>
//#include <future>



namespace fs = boost::filesystem;
using namespace std::complex_literals;
const auto PI=std::numbers::pi;
using wvVec=Eigen::VectorXcd;


class combineSegments {


public:
    /// This constructor initializes all parameters
    /// @param group group number of parameters
    /// @param row row number of parameters
    combineSegments(const int &group, const int &row){
        //
        this->groupNum = group;
        this->rowNum = row;
        this->outDir="./groupNew"+std::to_string(groupNum)+"/row"+std::to_string(rowNum)+"/";
        this->parseCSV(group, row);
        for (int n1 =0;n1<N1;n1++){
            this->x1ValsAll.push_back(-L1+dx1*n1);
        }
        for (int n2=0;n2<N2;n2++){
            this->x2ValsAll.push_back(-L2+dx2*n2);
        }
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
    int N1 = 0;
    int N2 = 0;
    double L1 = 0;
    double L2 = 0;
    double dx1=0;
    double dx2=0;
    double dtEst = 0.002;
    double tFlushStart = 0;
    double tFlushStop = 1;
    std::vector<int> jIndsAll;//index for time steps
    double tTotPerFlush=0;
    int stepsPerFlush=0;//total time steps
    double dt=0;
    std::vector<int> timeIndsAll;//all indices of time
    std::vector<double> x1ValsAll;
    std::vector<double> x2ValsAll;
    std::string outDir;
    int flushNum=0;

    std::vector<std::vector<std::complex<double>>> wvFunctions;//wavefunctions

    Eigen::SparseMatrix<std::complex<double>> NcMat;
    Eigen::SparseMatrix<std::complex<double>> NmMat;

    std::vector<wvVec> solutions;





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


    ///
    /// @param path path containing computational results of wavefunctions
    /// @return sorted files containing wavefunctions
    std::vector<std::string> listFiles(const std::string & path);

    ///
    /// @param binFiles files found by listFiles()
    /// @return sorted files by flush number
    std::vector<std::string> sortFiles(const std::vector<std::string>& binFiles);

    ///
    /// @param v vector
    /// @return argsort
    std::vector<size_t> argsort(const std::vector<int>& v) {
        // Initialize a vector of indices (0 to N-1)
        std::vector<size_t> indices(v.size());
        std::iota(indices.begin(), indices.end(), 0);

        // Sort indices based on comparing values in v
        std::sort(indices.begin(), indices.end(),
                  [&v](size_t i1, size_t i2) { return v[i1] < v[i2]; });

        return indices;
    }

    template<class T>
    static void printVec(const std::vector<T>& vec){
        for(int i=0;i<vec.size()-1;i++){
            std::cout<<vec[i]<<",";
        }
        std::cout<<vec[vec.size()-1]<<std::endl;
    }

    /// parse N1, N2, L1, L2
    /// @param binFileName file containing wavefunctions
    void catchParameters(const std::string& binFileName);

    ///initialize  dt
    void initTimeInds();

    /// read wavefunctions
    /// @param sortedFiles sorted bin files by flush number
    void readBinFiles(const std::vector<std::string>& sortedFiles);

    ///
    /// @param oneFileName one bin file name
    /// @return wavefunctions in this bin file
    std::vector<std::vector<std::complex<double>>> readOneBinFile(const std::string& oneFileName);



    ///fill in NcMat, NmMat
    void popolateMatrices();

    ///convert vector to eigen's vector
    std::vector<wvVec >  cppType2Eigen();

    ///
    /// @param solutionsInOneFile vectors containing in one file
    /// @return Eigen's vector
    std::vector<wvVec >cppType2EigenOneFile(const std::vector<std::vector<std::complex<double>>>& solutionsInOneFile);


    ///
    /// @param vec wavefunction
    /// @return photon number
    double numOfPhoton(const wvVec &vec);

    ///
    /// @return photon numbers at each time
    std::vector<double> photonAllParallel();

    ///
    /// @return photon numbers at each time
    std::vector<double> photonAllSerial();

    ///
    /// @param solutionsPerFlush solutions in one flush
    /// @return photon numbers in one flush
    std::vector<double> photonPerFlushSerial(const std::vector<wvVec >& solutionsPerFlush);


    ///
    /// @param vec wavefunction
    /// @return phonon number
    double numOfPhonon(const wvVec &vec);

    ///
    /// @return phonon numbers at each time
    std::vector<double> phononAllParallel();

    ///
    /// @return phonon numbers at each time
    std::vector<double> phononAllSerial();

    ///
    /// @param solutionsPerFlush solutions in one flush
    /// @return phonon numbers in one flush
    std::vector<double> phononPerFlushSerial(const std::vector<wvVec >& solutionsPerFlush);

    ///
    /// @param numVecvec photon/phonon numbers
    /// @return duplicated entries removed
    std::vector<double> removeHeadTail(const std::vector<std::vector<double>> &numVecvec);

    /// write to json
    /// @param photonNumAll all photon numbers
    /// @param phononNumAll all phonon numbers
    void to_json(const std::vector<double> &photonNumAll,const std::vector<double> &phononNumAll);

};


#endif //DCE_CPP_COMBINESEGMENTS_HPP
