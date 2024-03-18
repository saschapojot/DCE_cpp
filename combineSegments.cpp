//
// Created by polya on 3/12/24.
//

#include "combineSegments.hpp"


///
/// @param group group number
/// @param row row number
///parse csv with group number group
void combineSegments::parseCSV(const int &group, const int &row){
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




///
/// @param cmd python execution string
/// @return signal from the python
std::string combineSegments::execPython(const char *cmd){
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



///
/// @param path path containing computational results of wavefunctions
/// @return sorted files containing wavefunctions
std::vector<std::string> combineSegments::listFiles(const std::string & path){
    std::vector<std::string> retVec;
    try {
        // Check if the given path is a directory
        if (fs::is_directory(path)) {
            // Iterate over the directory entries
            for (const auto& entry : fs::directory_iterator(path)) {
                // Check if the directory entry is a file
                if (fs::is_regular_file(entry.status()) and entry.path().extension()==".bin") {
//                    std::cout << entry.path() << std::endl; // Print the file path
                    retVec.push_back(entry.path().filename().string());
                }
            }
        } else {
            std::cout << path << " is not a directory." << std::endl;
        }
    } catch (const fs::filesystem_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    this->flushNum=retVec.size();
    return retVec;


}


///
/// @param binFiles files found by listFiles()
/// @return sorted files by flush number
std::vector<std::string> combineSegments::sortFiles(const std::vector<std::string>& binFiles){

std::vector<int> flushVec;

for(const auto&fileName:binFiles){

    std::regex flushRegex("flush(\\d+)N1");
    std::smatch  matchFlush;
    std::regex_search(fileName,matchFlush,flushRegex);
    int fls=std::stoi(matchFlush.str(1));
    flushVec.push_back(fls);

}

std::vector<size_t> inds= argsort(flushVec);
std::vector<std::string> sortedFiles;
for(const auto& idx:inds){
    sortedFiles.push_back(this->outDir+ binFiles[idx]);
}
//    printVec(sortedFiles);
    return sortedFiles;


}

/// parse N1, N2, L1, L2
/// @param binFileName file containing wavefunctions
void combineSegments::catchParameters(const std::string& binFileName){
    std::regex N1Regex("N1(\\d+)N2");
    std::regex N2Regex("N2(\\d+)L1");

    std::regex L1Regex("L1(\\d*(\\.\\d+)?)L2");
    std::regex L2Regex("L2(\\d*(\\.\\d+)?)solution");

    std::smatch matchN1;
    std::smatch matchN2;
    std::smatch matchL1;
    std::smatch matchL2;

    std::regex_search(binFileName,matchN1,N1Regex);
    this->N1=std::stoi(matchN1.str(1));

    std::regex_search(binFileName,matchN2,N2Regex);
    this->N2=std::stoi(matchN2.str(1));

    std::regex_search(binFileName,matchL1,L1Regex);
    this->L1=std::stod(matchL1.str(1));

    std::regex_search(binFileName,matchL2,L2Regex);
    this->L2=std::stod(matchL2.str(1));
    this->dx1=2*L1/(static_cast<double>(N1));
     this->dx2=2*L2/(static_cast<double >(N2));

//    std::cout<<"N1="<<N1<<", N2="<<N2<<", L1="<<L1<<", L2="<<L2<<std::endl;

    this->parseCSV(groupNum, rowNum);
    for (int n1 =0;n1<N1;n1++){
        this->x1ValsAll.push_back(-L1+dx1*n1);
    }
    for (int n2=0;n2<N2;n2++){
        this->x2ValsAll.push_back(-L2+dx2*n2);
    }

}

///initialize time indices and dt
void combineSegments::initTimeInds(){

    this->tTotPerFlush=this->tFlushStop-this->tFlushStart;
    this->stepsPerFlush=static_cast<int>(std::ceil(tTotPerFlush/dtEst));
    this->dt=tTotPerFlush/static_cast<double >(stepsPerFlush);

//    for(int fls=0;fls<flushNum;fls++){
//        int startingInd=fls*stepsPerFlush;
//        for(int j=0;j<stepsPerFlush;j++){
//            this->timeIndsAll.push_back(startingInd+j);
//        }
//    }

}


/// read wavefunctions
/// @param sortedFiles sorted bin files by flush number
void combineSegments::readBinFiles(const std::vector<std::string>& sortedFiles) {
    for (const auto &filePath: sortedFiles) {
        std::ifstream file(filePath, std::ios::binary);
        if (!file.is_open()) {
            std::cerr << "Failed to open file: " << filePath << std::endl;
            exit(1);
        }

        // Read the entire file into a std::vector<char>
        std::vector<char> buffer((std::istreambuf_iterator<char>(file)),
                                 std::istreambuf_iterator<char>());

        // Now, use msgpack to unpack the buffer
        msgpack::object_handle oh =
                msgpack::unpack(buffer.data(), buffer.size());
        // Convert the unpacked object to your desired C++ data structure
        std::vector<std::vector<std::complex<double>>> deserializedData;
        oh.get().convert(deserializedData);

        for (const auto &vec: deserializedData) {
            this->wvFunctions.push_back(vec);
        }


    }
//    std::cout << wvFunctions.size() << std::endl;

}


///fill in NcMat, NmMat
void combineSegments::popolateMatrices() {
    //NcMat
    Eigen::SparseMatrix<std::complex<double>> NcPart0 = Eigen::SparseMatrix<std::complex<double>>(N1 * N2, N1 * N2);

//    for (int n1 = 0; n1 < N1; n1++) {
//        double x1n1Tmp2 = std::pow(x1ValsAll[n1], 2);
//        int startingPos = n1 * N2;
//        for (int n2 = 0; n2 < N2; n2++) {
//            NcPart0.insert(startingPos + n2, startingPos + n2) = x1n1Tmp2;
//        }
//    }
std::vector<double> x1Squared(x1ValsAll.size(),0);
for(int n1=0;n1<x1ValsAll.size();n1++){
    x1Squared[n1]=std::pow(x1ValsAll[n1],2);
}

        for(int i=0;i<N1*N2;i++){
            int n1= static_cast<int>(std::floor(static_cast<double >(i)/static_cast<double >(N2)));
            NcPart0.insert(i,i)=x1Squared[n1];
        }

    NcPart0 *= 0.5 * omegac;

    for (int i = 0; i < N1 * N2; i++) {
        NcPart0.coeffRef(i, i) += -0.5;
    }


    Eigen::SparseMatrix<std::complex<double>> H6 = Eigen::SparseMatrix<std::complex<double>>(N1 * N2, N1 * N2);
    H6 = Eigen::SparseMatrix<std::complex<double>>(N1 * N2, N1 * N2);
    //fill in diagonal values
    for (int i = 0; i < N1 * N2; i++) {
        H6.insert(i, i) = -2.0;
    }
    //fill in super-block-diagonal
    Eigen::SparseMatrix<std::complex<double>> H6Upper(N1 * N2, N1 * N2);
    for (int i = 0; i < (N1 - 1) * N2; i++) {
        H6Upper.insert(i, i + N2) = 1.0;
    }
    Eigen::SparseMatrix<std::complex<double>> H6Lower = H6Upper.transpose();
    H6 += H6Upper;
    H6 += H6Lower;
//    std::cout<<"H6 block="<<H6<<std::endl;
    H6 *= -1.0 / (2.0 * std::pow(dx1, 2));
//        std::cout<<"H6 block="<<H6<<std::endl;
    this->NcMat = NcPart0 + H6 / omegac;

//NmMat
    std::vector<double> x2ValsSquared;
    for (const auto &val: x2ValsAll) {
        x2ValsSquared.push_back(std::pow(val, 2));
    }
    Eigen::SparseMatrix<std::complex<double>> NmPart0 = Eigen::SparseMatrix<std::complex<double>>(N1 * N2, N1 * N2);
//    for (int n1 = 0; n1 < N1; n1++) {
//        int startingPos = n1 * N2;
//        for (int n2 = 0; n2 < N2; n2++) {
//            NmPart0.insert(startingPos + n2, startingPos + n2) = x2ValsSquared[n2];
//        }
//    }
    //S2 in diagonal
    for(int i=0;i<N1*N2;i++){
        int n2=i%N2;
        NmPart0.insert(i,i)=x2ValsSquared[n2];
    }


    NmPart0 *= 0.5 * omegam;

    for(int i=0;i<N1*N2;i++){
        NmPart0.coeffRef(i,i)+=-0.5;
    }

    Eigen::SparseMatrix<std::complex<double>> NmPart1 = Eigen::SparseMatrix<std::complex<double>>(N1 * N2, N1 * N2);
//    for(int n1=0;n1<N1;n1++){
//        int startingPos = n1 * N2;
//        //diagonal
//        for(int n2=0;n2<N2;n2++){
//            NmPart1.insert(startingPos+n2,startingPos+n2)=-2.0;
//        }
//        //superdiagonal
//        for(int n2=0;n2<N2-1;n2++){
//            NmPart1.insert(startingPos+n2,startingPos+n2+1)=1.0;
//        }
//        //subdiagonal
//        for(int n2=0;n2<N2-1;n2++){
//            NmPart1.insert(startingPos+n2+1,startingPos+n2)=1.0;
//        }
//    }
std::vector<int> upperRowInds;
upperRowInds.reserve((N2-1)*N1);
    for(int n2=0;n2<N1*N2;n2++){
        if ((n2+1)%N2==0){
            continue;
        }else{
            upperRowInds.push_back(n2);
        }
    }
    //diagonal
    for(int i=0;i<N1*N2;i++){
        NmPart1.insert(i,i)=-2.0;
    }
    //upper and lower diagonal
    for(const auto&n2:upperRowInds){
        NmPart1.insert(n2,n2+1)=1.0;
        NmPart1.insert(n2+1,n2)=1.0;
    }

    NmPart1*=-1/(2.0*omegam*std::pow(dx2,2));

    this->NmMat=NmPart0+NmPart1;





}



///convert vector to eigen's vector
std::vector<wvVec >  combineSegments::cppType2Eigen() {

    std::vector<wvVec> retVec;
    for (const auto &vec: this->wvFunctions) {
        int length = vec.size();
        wvVec vecEigenTmp = wvVec(length);
       for(int i=0;i<length;i++){
           vecEigenTmp(i)=vec[i];
       }
//        std::cout<<"deserialized: "<<vecEigenTmp<<std::endl;
        retVec.push_back(vecEigenTmp);
    }

    return retVec;

}

///
/// @param solutionsInOneFile vectors containing in one file
/// @return Eigen's vector
std::vector<wvVec >combineSegments::cppType2EigenOneFile(const std::vector<std::vector<std::complex<double>>>& solutionsInOneFile){

    std::vector<wvVec> retVec;
    for(const auto&vec:solutionsInOneFile){
        int length=vec.size();
        wvVec vecEigenTmp = wvVec(length);
        for(int i=0;i<length;i++){
            vecEigenTmp(i)=vec[i];
        }
        retVec.push_back(vecEigenTmp);
    }

    return retVec;


}



///
/// @param vec wavefunction
/// @return photon number
double combineSegments::numOfPhoton(const wvVec &vec){

    double numPht=std::abs(vec.dot(this->NcMat*vec));
    return numPht;


}

///
/// @param vec wavefunction
/// @return phonon number
double combineSegments::numOfPhonon(const wvVec &vec){

    double numPhn=std::abs(vec.dot(this->NmMat*vec));
    return numPhn;
}



///
/// @return photon numbers at each time
//std::vector<double> combineSegments::photonAllParallel(){
//
//    std::vector<std::future<double>> futures;
//    for(const auto&vec:this->solutions){
//        std::future<double> fut=std::async(std::launch::async,&combineSegments::numOfPhoton,this,vec);
//        futures.push_back(std::move(fut));
//    }
//    std::vector<double>retVec;
//    for(auto& fut:futures){
//        retVec.push_back(fut.get());
//    }
//    return retVec;
//
//
//}



///
/// @return phonon numbers at each time
//std::vector<double> combineSegments::phononAllParallel(){
//    std::vector<std::future<double>> futures;
//    for(const auto&vec:this->solutions){
//        std::future<double> fut=std::async(std::launch::async,&combineSegments::numOfPhonon,this,vec);
//        futures.push_back(std::move(fut));
//    }
//
//
//    std::vector<double>retVec;
//    for(auto& fut:futures){
//        retVec.push_back(fut.get());
//    }
//    return retVec;
//
//
//
//}

///
/// @return photon numbers at each time
std::vector<double> combineSegments::photonAllSerial(){
    std::vector<double>retVec;
    for(const auto&vec:this->solutions){
        retVec.push_back(this->numOfPhoton(vec));
    }
    return retVec;

}

///
/// @return phonon numbers at each time
std::vector<double> combineSegments::phononAllSerial(){
    std::vector<double>retVec;
    for(const auto&vec:this->solutions){
        retVec.push_back(this->numOfPhonon(vec));
    }
    return retVec;

}



/// write to csv
/// @param photonNumAll all photon numbers
/// @param phononNumAll all phonon numbers
void combineSegments::to_json(const std::vector<double> &photonNumAll,const std::vector<double> &phononNumAll){

    std::string suffix="N1"+std::to_string(N1)+"N2"+std::to_string(N2)+"L1"+std::to_string(L1)+"L2"+std::to_string(L2);
std::string outPhotonFileName=this->outDir+"photon"+suffix+".json";

std::string outPhononFileName=this->outDir+"phonon"+suffix+".json";

boost::json::object objPhoton;
boost::json::array arrPhoton;
for(const auto &val:photonNumAll){
    arrPhoton.push_back(val);
}

    boost::json::array timeAll;
for(int i=0;i<photonNumAll.size();i++){
    double ti=static_cast<double >(i)*dt;
    timeAll.push_back(ti);
}



objPhoton["photonNum"]=arrPhoton;
objPhoton["time"]=timeAll;


    boost::json::object objPhonon;
    boost::json::array arrPhonon;
    for(const auto&val:phononNumAll){
        arrPhonon.push_back(val);
    }

    objPhonon["phononNum"]=arrPhonon;
    objPhonon["time"]=timeAll;

    std::ofstream ofsPhoton(outPhotonFileName);
    std::string photon_str=boost::json::serialize(objPhoton);
    ofsPhoton<<photon_str<<std::endl;
    ofsPhoton.close();


    std::ofstream  ofsPhonon(outPhononFileName);
    std::string phonon_str=boost::json::serialize(objPhonon);
    ofsPhonon<<phonon_str<<std::endl;
    ofsPhonon.close();

}


///
/// @param oneFileName one bin file name
/// @return wavefunctions in this bin file
std::vector<std::vector<std::complex<double>>> combineSegments::readOneBinFile(const std::string& oneFileName){

std::ifstream file(oneFileName,std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << oneFileName << std::endl;
        exit(1);
    }
    // Read the entire file into a std::vector<char>

    std::vector<char> buffer((std::istreambuf_iterator<char>(file)),
                             std::istreambuf_iterator<char>());
    // Now, use msgpack to unpack the buffer

    msgpack::object_handle oh =
            msgpack::unpack(buffer.data(), buffer.size());
    // Convert the unpacked object to your desired C++ data structure
    std::vector<std::vector<std::complex<double>>> deserializedData;
    oh.get().convert(deserializedData);

    std::vector<std::vector<std::complex<double>>> retVec;
    for(const auto& vec:deserializedData){
        retVec.push_back(vec);
    }
    return retVec;



}


///
/// @param solutionsPerFlush solutions in one flush
/// @return photon numbers in one flush
std::vector<double> combineSegments::photonPerFlushSerial(const std::vector<wvVec >& solutionsPerFlush) {

    std::vector<double> retVec;
    for (const auto &vec: solutionsPerFlush) {
        retVec.push_back(this->numOfPhoton(vec));
    }
    return retVec;


}



///
/// @param solutionsPerFlush solutions in one flush
/// @return phonon numbers in one flush
std::vector<double> combineSegments::phononPerFlushSerial(const std::vector<wvVec >& solutionsPerFlush){
    std::vector<double> retVec;

    for(const auto&vec:solutionsPerFlush){
        retVec.push_back(this->numOfPhonon(vec));
    }
    return retVec;


}
///
/// @param numVecvec photon/phonon numbers
/// @return duplicated entries removed
std::vector<double> combineSegments::removeHeadTail(const std::vector<std::vector<double>> &numVecvec){
    std::vector<double> retVec;
    for(int i=0;i<numVecvec.size()-1;i++){
        for(int j=0;j<numVecvec[i].size()-1;j++){
            retVec.push_back(numVecvec[i][j]);
        }
    }

    for(int j=0;j<numVecvec[numVecvec.size()-1].size();j++){
        retVec.push_back(numVecvec[numVecvec.size()-1][j]);
    }

    return retVec;
}


///write the init wavefunction and last wavefunction to xml file
///
/// @param vec initial or final wavefunction
void combineSegments::writeInitLast(const std::string &filename,std::vector<std::complex<double>>& vec){

    std::ofstream ofs(filename);
    boost::archive::xml_oarchive oa(ofs);
    oa & BOOST_SERIALIZATION_NVP(vec);

}