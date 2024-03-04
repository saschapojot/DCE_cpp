//
// Created by polya on 2/22/24.
//

#include "timeEvolution.hpp"

///
/// @param cmd python execution string
/// @return signal from the python
void DCE_Evolution::parseCSV(const int &group, const int &row) {

std::string commandToReadCSV="python3 readCSV.py "+std::to_string(group)+" "+std::to_string(row);

std::string result=this->execPython(commandToReadCSV.c_str());
//    std::cout<<result<<std::endl;

    std::regex pattern_j1H("j1H(\\d+)j2H");
    std::smatch  match_j1H;
    if (std::regex_search(result,match_j1H,pattern_j1H)){
        this->jH1=std::stoi(match_j1H[1].str());
    }

    std::regex pattern_j2H("j2H(\\d+)g0");
    std::smatch match_j2H;
    if (std::regex_search(result,match_j2H,pattern_j2H)){
        this->jH2=std::stoi(match_j2H[1].str());
    }

    std::regex pattern_g0("g0([+-]?\\d+(\\.\\d+)?)omegam");
    std::smatch match_g0;
    if (std::regex_search(result,match_g0,pattern_g0)){
        this->g0=std::stod(match_g0[1].str());
    }


    std::regex pattern_omegam("omegam([+-]?\\d+(\\.\\d+)?)omegap");
    std::smatch match_omegam;
    if (std::regex_search(result,match_omegam,pattern_omegam)){
        this->omegam=std::stod(match_omegam[1].str());
    }


    std::regex pattern_omegap("omegap([+-]?\\d+(\\.\\d+)?)omegac");
    std::smatch match_omegap;
    if (std::regex_search(result,match_omegap,pattern_omegap)){
        this->omegap=std::stod(match_omegap[1].str());
    }

    std::regex pattern_omegac("omegac([+-]?\\d+(\\.\\d+)?)er");
    std::smatch match_omegac;
    if (std::regex_search(result,match_omegac,pattern_omegac)){
        this->omegac=std::stod(match_omegac[1].str());
    }

    std::regex pattern_er("er([+-]?\\d+(\\.\\d+)?)thetaCoef");
    std::smatch match_er;
    if(std::regex_search(result,match_er,pattern_er)){
        this->er=std::stod(match_er[1].str());
    }

    std::regex pattern_thetaCoef("thetaCoef([+-]?\\d+(\\.\\d+)?)");
    std::smatch  match_thetaCoef;
    if (std::regex_search(result,match_thetaCoef,pattern_thetaCoef)){
        this->thetaCoef=std::stod(match_thetaCoef[1].str());
    }
//    std::cout<<"jH1="<<jH1<<std::endl;
//
//    std::cout<<"jH2="<<jH2<<std::endl;
//
//    std::cout<<"g0="<<g0<<std::endl;
//
//    std::cout<<"omegam="<<omegam<<std::endl;
//
//    std::cout<<"omegac="<<omegac<<std::endl;
//
//    std::cout<<"omegap="<<omegap<<std::endl;
//
//    std::cout<<"er="<<er<<std::endl;
//
//    std::cout<<"thetaCoef="<<thetaCoef<<std::endl;
      double e2r=std::pow(er,2);
      double eM2r=1/e2r;
      this->Deltam=this->omegam-this->omegap;
      this->lmd=(e2r-eM2r)/(e2r+eM2r)*Deltam;
      this->theta=thetaCoef*PI;
//      std::cout<<"lambda="<<lmd<<std::endl;
//      std::cout<<"theta="<<theta<<std::endl;
//      std::cout<<"Deltam="<<Deltam<<std::endl;



}


std::string DCE_Evolution::execPython(const char *cmd) {
    std::array<char, 4096> buffer; // Buffer to store command output
    std::string result; // String to accumulate output

    // Open a pipe to read the output of the executed command
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }

    // Read the output a chunk at a time and append it to the result string
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }

    return result; // Return the accumulated output

}



///This function initializes the sparse matrices
void DCE_Evolution::populatedMatrices(){
    //initialize identity matrix
    this->IN1N2=Eigen::SparseMatrix<std::complex<double>>(N1*N2,N1*N2);
    for (int i=0;i<N1*N2;i++){
        this->IN1N2.insert(i,i)=1.0;
    }

//    this->D2=Eigen::SparseMatrix<std::complex<double>>(N2,N2);
//    for(int n2=0;n2<N2;n2++){
//        D2.insert(n2,n2)=this->x2ValsAll[n2];
//    }

    //initialize H0
    this->H0=(-0.5*omegac-0.5*Deltam)*this->IN1N2;

    //initialize H2
    this->H2=Eigen::SparseMatrix<std::complex<double>>(N1*N2,N1*N2);
    for(int n1=0;n1<N1;n1++) {

        for (int n2 = 0; n2 < N2; n2++) {
            this->H2.insert(n1 * N2 + n2, n1 * N2 + n2) = std::pow(this->x1ValsAll[n1], 2);

        }
    }
    H2*=0.5*std::pow(omegac,2);

    //initialize H3
    this->H3=Eigen::SparseMatrix<std::complex<double>>(N1*N2,N1*N2);

    std::vector<std::complex<double>> S2Diag;
    for (int n2=0;n2<N2;n2++){
        S2Diag.push_back(std::pow(this->x2ValsAll[n2],2));
    }

    for(int n1=0;n1<N1;n1++){
        for(int n2=0;n2<N2;n2++){

            this->H3.insert(n1*N2+n2,n1*N2+n2)=S2Diag[n2];
        }

    }
    H3*=(0.5*lmd*std::cos(theta)+0.5*Deltam*omegam);

    //initialize H6
    this->H6=Eigen::SparseMatrix<std::complex<double>>(N1*N2,N1*N2);
    //fill in diagonal values
    for (int i=0;i<N1*N2;i++){
        this->H6.insert(i,i)=-2.0;
    }
    //fill in super-block-diagonal
    Eigen::SparseMatrix<std::complex<double>> H6Upper(N1*N2,N1*N2);
    for (int i=0;i<(N1-1)*N2;i++){
        H6Upper.insert(i,i+N2)=1.0;
    }
    Eigen::SparseMatrix<std::complex<double>> H6Lower=H6Upper.transpose();
    H6+=H6Upper;
    H6+=H6Lower;
//    std::cout<<"H6 block="<<H6<<std::endl;
    H6*=-1.0/(std::pow(dx1,2));
//        std::cout<<"H6 block="<<H6<<std::endl;


}//end of populate matrices

