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



///This function initializes the sparse matrices H0,H2,H3,H6,H7,H8
void DCE_Evolution::populatedMatrices(){
    const auto tPopStart{std::chrono::steady_clock::now()};

//    const auto tIN1N2Start{std::chrono::steady_clock::now()};
    //initialize identity matrix
    this->IN1N2=Eigen::SparseMatrix<std::complex<double>>(N1*N2,N1*N2);
    for (int i=0;i<N1*N2;i++){
        this->IN1N2.insert(i,i)=1.0;
    }
//    const auto tIN1N2End{std::chrono::steady_clock::now()};
//    const std::chrono::duration<double> elapsed_IN1N2{tIN1N2End- tIN1N2Start};
//    std::cout<<"IN1N2 time: "<< elapsed_IN1N2.count() / 3600.0 << " h" << std::endl;


//    const auto tH0Start{std::chrono::steady_clock::now()};
    //initialize H0
    this->H0=(-0.5*omegac-0.5*Deltam)*this->IN1N2;
//    const auto tH0End{std::chrono::steady_clock::now()};
//    const std::chrono::duration<double> elapsed_H0{tH0End- tH0Start};
//    std::cout<<"H0 time: "<< elapsed_H0.count() / 3600.0 << " h" << std::endl;

    //initialize H2
//    const auto tH2Start{std::chrono::steady_clock::now()};
    this->H2=Eigen::SparseMatrix<std::complex<double>>(N1*N2,N1*N2);
    for(int n1=0;n1<N1;n1++) {

        for (int n2 = 0; n2 < N2; n2++) {
            this->H2.insert(n1 * N2 + n2, n1 * N2 + n2) = std::pow(this->x1ValsAll[n1], 2);

        }
    }
    H2*=0.5*std::pow(omegac,2);
//    const auto tH2End{std::chrono::steady_clock::now()};
//    const std::chrono::duration<double> elapsed_H2{tH2End- tH2Start};
//    std::cout<<"H2 time: "<< elapsed_H2.count() / 3600.0 << " h" << std::endl;

    //initialize H3
//    const auto tH3Start{std::chrono::steady_clock::now()};
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
    H3*=(0.5*lmd*omegam*std::cos(theta)+0.5*Deltam*omegam);
//    const auto tH3End{std::chrono::steady_clock::now()};
//    const std::chrono::duration<double> elapsed_H3{tH3End- tH3Start};
//    std::cout<<"H3 time: "<< elapsed_H3.count() / 3600.0 << " h" << std::endl;

    //initialize H6
//    const auto tH6Start{std::chrono::steady_clock::now()};
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
    H6*=-1.0/(2.0*std::pow(dx1,2));
//        std::cout<<"H6 block="<<H6<<std::endl;
//    const auto tH6End{std::chrono::steady_clock::now()};
//    const std::chrono::duration<double> elapsed_H6{tH6End- tH6Start};
//    std::cout<<"H6 time: "<< elapsed_H6.count() / 3600.0 << " h" << std::endl;

    // initialize H7
//    const auto tH7Start{std::chrono::steady_clock::now()};
    this->H7=Eigen::SparseMatrix<std::complex<double>>(N1*N2,N1*N2);

    //H7's upper and lower diagonal indices
    std::vector<int> H7UpperRowInds;
    H7UpperRowInds.reserve((N2-1)*N1);
    for(int n2=0;n2<N1*N2;n2++){
        if ((n2+1)%N2==0){
            continue;
        }else{
            H7UpperRowInds.push_back(n2);
        }
    }
    //H7's diagonal
    for(int i=0;i<N1*N2;i++){
        H7.insert(i,i)=-2.0;
    }
    //H7's upper and lower diagonal
    for(const auto &n2:H7UpperRowInds){
        H7.insert(n2,n2+1)=1.0;
        H7.insert(n2+1,n2)=1.0;
    }
    H7*=(-Deltam/(2*omegam)+lmd*std::cos(theta)/(2*omegam))/(std::pow(dx2,2));
//    const auto tH7End{std::chrono::steady_clock::now()};
//    const std::chrono::duration<double> elapsed_H7{tH7End- tH7Start};
//    std::cout<<"H7 time: "<< elapsed_H7.count() / 3600.0 << " h" << std::endl;

    //initialze H8
//    const auto tH8Start{std::chrono::steady_clock::now()};
    this->H8=Eigen::SparseMatrix<std::complex<double>>(N1*N2,N1*N2);
    std::vector<double> upperDiagVals;
    for(int n2=0;n2<N2-1;n2++){
        upperDiagVals.push_back(x2ValsAll[n2]+x2ValsAll[n2+1]);
    }

    for(int n1=0;n1<N1;n1++) {
        int startingRowPos = n1 * N2;
        int startingColPos = n1 * N2;
        //upper superdiagonal
        for(int n2=0;n2<N2-1;n2++){
            H8.insert(startingRowPos+n2,startingColPos+n2+1)=upperDiagVals[n2];
            H8.insert(startingRowPos+n2+1,startingColPos+n2)=-upperDiagVals[n2];
        }

    }
    H8*=1i*lmd*std::sin(theta)/(4*dx2);
//    const auto tH8End{std::chrono::steady_clock::now()};
//    const std::chrono::duration<double> elapsed_H8{tH8End- tH8Start};
//    std::cout<<"H8 time: "<< elapsed_H8.count() / 3600.0 << " h" << std::endl;

//    const auto tHSumStart{std::chrono::steady_clock::now()};
    this->HSumStatic=H0+H2+H3+H6+H7+H8;
//    const auto tHSumEnd{std::chrono::steady_clock::now()};
//    const std::chrono::duration<double> elapsed_HSum{tHSumEnd- tHSumStart};
//    std::cout<<"HSum time: "<< elapsed_HSum.count() / 3600.0 << " h" << std::endl;

    const auto tPopEnd{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_secondsAll{tPopEnd - tPopStart};
    std::cout<<"populate matrices time: "<< elapsed_secondsAll.count() / 3600.0 << " h" << std::endl;

//    std::cout<<"HS="<<HSumStatic<<std::endl;


}//end of populate matrices



///
/// @param j index of time
/// @return H1j
Eigen::SparseMatrix<std::complex<double>> DCE_Evolution::H1Val(const int & j){

    Eigen::SparseMatrix<std::complex<double>> retMat=Eigen::SparseMatrix<std::complex<double>>(N1*N2,N1*N2);
    //fill in diagonal
    for (int i=0;i<N1*N2;i++){
        int n2=i%N2;
        retMat.insert(i,i)=x2ValsAll[n2];
    }

//    std::cout<<retMat<<std::endl;

    double tj=static_cast<double >(j)*dt;

    retMat*=-0.5*g0*std::sqrt(2.0*omegam)*std::cos(omegap*(tj+0.5*dt));

    return retMat;



}




///
/// @param j j index of time
/// @return H4j
Eigen::SparseMatrix<std::complex<double>> DCE_Evolution::H4Val(const int & j){
    Eigen::SparseMatrix<std::complex<double>> retMat=Eigen::SparseMatrix<std::complex<double>>(N1*N2,N1*N2);

    for(int n1=0;n1<N1;n1++){
        int startingPos=n1*N2;
        double x1Tmp=x1ValsAll[n1];
        double x1Tmp2=std::pow(x1Tmp,2);
        for(int n2=0;n2<N2;n2++){
            retMat.insert(startingPos+n2,startingPos+n2)=x1Tmp2*x2ValsAll[n2];
        }
    }
    double tj=static_cast<double >(j)*dt;

    retMat*=g0*omegac*std::sqrt(2.0*omegam)*std::cos(omegap*(tj+0.5*dt));

    return retMat;



}



///
/// @param j index of time
/// @return H5j
Eigen::SparseMatrix<std::complex<double>>DCE_Evolution::H5Val(const int& j){
    //construct A5j
    Eigen::SparseMatrix<std::complex<double>> A5jPart0=Eigen::SparseMatrix<std::complex<double>>(N1*N2,N1*N2);
    Eigen::SparseMatrix<std::complex<double>> A5jPart1=Eigen::SparseMatrix<std::complex<double>>(N1*N2,N1*N2);

    for(int n1=0;n1<N1;n1++){
        int startingPos=n1*N2;
        double x1Tmp=x1ValsAll[n1];
        double x1Tmp2=std::pow(x1Tmp,2);
        for(int n2=0;n2<N2;n2++){
            A5jPart0.insert(startingPos+n2,startingPos+n2)=x1Tmp2;

        }
    }
    double tj=static_cast<double >(j)*dt;

    A5jPart0*=-1i*g0*omegac*std::sqrt(2.0/omegam)*std::sin(omegap*(tj+0.5*dt));

    for(int i=0;i<N1*N2;i++){
        A5jPart1.insert(i,i)=1.0;

    }
    A5jPart1*=1j*0.5*g0*std::sqrt(2.0/omegam)*std::sin(omegap*(tj+0.5*dt));

    Eigen::SparseMatrix<std::complex<double>> A5j=A5jPart0+A5jPart1;


    //construct diag P2
    Eigen::SparseMatrix<std::complex<double>> diagPart=Eigen::SparseMatrix<std::complex<double>>(N1*N2,N1*N2);

    for (int n1=0;n1<N1;n1++){
        int startingPos=n1*N2;
        for(int n2=0;n2<N2-1;n2++){
            diagPart.insert(startingPos+n2,startingPos+n2+1)=1.0;
            diagPart.insert(startingPos+n2+1,startingPos+n2)=-1.0;
        }
    }


    Eigen::SparseMatrix<std::complex<double>> retMat=A5j*diagPart/(2*dx2);

    return retMat;
}


///
/// @param j index of time
/// @return HDj
Eigen::SparseMatrix<std::complex<double>> DCE_Evolution::HDj(const int &j){


    return HSumStatic+ H1Val(j)+ H4Val(j)+ H5Val(j);


}

///initialize wavefunction
void DCE_Evolution::initPsiSerial(){
    const auto tInitPsiStart{std::chrono::steady_clock::now()};

    std::vector<double> valVec;
    for(int i=0;i<N1*N2;i++){
        valVec.push_back(this->fillOneValInPsi(i));
    }

    this->Psi0=wvVec (N1*N2);
    for(int i=0;i<N1*N2;i++){
        this->Psi0(i)=valVec[i];
    }

    double nm0=Psi0.norm();
    Psi0/=nm0;
    const auto tInitPsiEnd{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_secondsAll{tInitPsiEnd - tInitPsiStart};
    std::cout<<"init Psi time: "<< elapsed_secondsAll.count() / 3600.0 << " h" << std::endl;





}


///initialize wavefunction in parallel
void DCE_Evolution::initPsiParallel() {
    const auto tInitPsiStart{std::chrono::steady_clock::now()};
    std::vector<std::future<double>> futures;
    for (int i = 0; i < N1 * N2; i++) {
        std::future<double> fut = std::async(std::launch::async, &DCE_Evolution::fillOneValInPsi, this, i);
        futures.push_back(std::move(fut));
    }

    std::vector<double> retVec;
    for (auto &fut: futures) {
        retVec.push_back(fut.get());
    }


    this->Psi0 = wvVec(N1 * N2);
    for (int i = 0; i < N1 * N2; i++) {
        this->Psi0(i) = retVec[i];
    }
    double nm0=Psi0.norm();
    Psi0/=nm0;
    const auto tInitPsiEnd{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_secondsAll{tInitPsiEnd - tInitPsiStart};
    std::cout<<"init Psi time: "<< elapsed_secondsAll.count() / 3600.0 << " h" << std::endl;
}

///create output directories
void DCE_Evolution::createOutDir(){
    const auto tCreatDirStart{std::chrono::steady_clock::now()};
    this->outDir="./groupNew"+std::to_string(this->groupNum)+"/row"+std::to_string(this->rowNum)+"/";

    if(!fs::is_directory(outDir) || !fs::exists(outDir)){
        fs::create_directories(outDir);
    }

    const auto tCreatDirEnd{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_secondsAll{tCreatDirEnd - tCreatDirStart};
    std::cout<<"create dir time: "<< elapsed_secondsAll.count() / 3600.0 << " h" << std::endl;



}



///initialize time indices and dt
void DCE_Evolution::initTimeInds(){
    const auto tInitTimeStart{std::chrono::steady_clock::now()};
    this->tTotPerFlush=this->tFlushStop-this->tFlushStart;
    this->stepsPerFlush=static_cast<int>(std::ceil(tTotPerFlush/dtEst));
    this->dt=tTotPerFlush/static_cast<double >(stepsPerFlush);

    for(int fls=0;fls<flushNum;fls++){
        int startingInd=fls*stepsPerFlush;
        for(int j=0;j<stepsPerFlush;j++){
            this->timeIndsAll.push_back(startingInd+j);
        }
    }

//std::vector<double> tAll;
//    for (const auto&val:timeIndsAll){
//        tAll.push_back(double (val)*dt);
//    }
//
//    printVec(tAll);

    const auto tInitTimeEnd{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_secondsAll{tInitTimeEnd - tInitTimeStart};
    std::cout<<"init time inds time: "<< elapsed_secondsAll.count() / 3600.0 << " h" << std::endl;


}



//evolution and write to file by flush
wvVec DCE_Evolution::evolutionPerFlush(const int &fls, const wvVec& initVec){

    int startingInd=fls*this->stepsPerFlush;

    int nextStartingInd=startingInd+stepsPerFlush;
    std::cout<<"starting loop: "<<startingInd<<", ending loop: "<<nextStartingInd-1<<std::endl;

    const auto tStart{std::chrono::steady_clock::now()};
    std::unique_ptr<std::vector<wvVec>> PsiPerFlushPtr=std::make_unique<std::vector<wvVec>>();
    PsiPerFlushPtr->resize(stepsPerFlush+1);

   for(auto &vec:*PsiPerFlushPtr){
       vec=wvVec (N1*N2);
   }
    (*PsiPerFlushPtr)[0]=initVec;
    wvVec PsiCurr=initVec;
    int counter=1;
    for(int j=startingInd;j<nextStartingInd;j++){

        wvVec PsiNext= oneStepEvolution(j,PsiCurr);
        (*PsiPerFlushPtr)[counter]=(PsiNext);
        counter++;
        PsiCurr=PsiNext;


    }

    const auto tEnd{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_secondsAll{tEnd - tStart};
    std::cout<<"flush "<<fls<<" time: "<< elapsed_secondsAll.count() / 3600.0 << " h" << std::endl;

    const auto tOutStart{std::chrono::steady_clock::now()};
    std::vector<std::vector<std::complex<double>>> outData=this->eigen2cppType(PsiPerFlushPtr);


    std::string outFile=this->outDir+"flush"+std::to_string(fls)+"N1"+std::to_string(N1)
            +"N2"+std::to_string(N2)+"L1"+std::to_string(L1)
            +"L2"+std::to_string(L2)+"solution.bin";
    std::ofstream ofs(outFile,std::ios::binary);
    msgpack::pack(ofs,outData);
    ofs.close();

    const auto tOutEnd{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_Out{tOutEnd - tOutStart};
    std::cout<<"Out time: "<< elapsed_Out.count() / 3600.0 << " h" << std::endl;

    int length=(*PsiPerFlushPtr).size();


    return (*PsiPerFlushPtr)[length-1];

}



///
/// @param j time step
/// @param PsiCurr wavefunction before evolution
/// @return wavefunction after evolution
wvVec DCE_Evolution::oneStepEvolution(const int& j, const wvVec& PsiCurr){

//    const auto tInitHStart{std::chrono::steady_clock::now()};
    Eigen::SparseMatrix<std::complex<double>> HDjMat= HDj(j);
//    const auto tInitHEnd{std::chrono::steady_clock::now()};
//    const std::chrono::duration<double> elapsed_secondsAll{tInitHEnd - tInitHStart};
//    std::cout<<"init H: "<< elapsed_secondsAll.count() / 3600.0 << " h" << std::endl;


//    const auto tSolve_yStart{std::chrono::steady_clock::now()};
    Eigen::SparseMatrix<std::complex<double>> mat0Tmp=0.5*1i*dt*HDjMat;
    //add scalar 1 to matTmp, i.e., add identity matrix to matTmp
    for(int i=0;i<N1*N2;i++){
        mat0Tmp.coeffRef(i,i)+=1.0;
    }

    Eigen::BiCGSTAB<Eigen::SparseMatrix<std::complex<double>>> solver;
    solver.compute(mat0Tmp);
    wvVec y=solver.solve(PsiCurr);
//    const auto tSolve_yEnd{std::chrono::steady_clock::now()};
//    const std::chrono::duration<double> elapsed_Solve_y{tSolve_yEnd- tSolve_yStart};
//    std::cout<<"solve y: "<< elapsed_Solve_y.count() / 3600.0 << " h" << std::endl;

//    const auto tProdStart{std::chrono::steady_clock::now()};
    Eigen::SparseMatrix<std::complex<double>> mat1Tmp=-0.5*1i*dt*HDjMat;
    for(int i=0;i<N1*N2;i++){
        mat1Tmp.coeffRef(i,i)+=1.0;
    }

    wvVec PsiNext=mat1Tmp*y;
//    const auto tProdEnd{std::chrono::steady_clock::now()};
//    const std::chrono::duration<double> elapsedProd{tProdEnd- tProdStart};
//    std::cout<<"prod: "<< elapsedProd.count() / 3600.0 << " h" << std::endl;
    return PsiNext;


}


///
/// @param solutions solutions per flush
/// @return converted to cpp data type
std::vector<std::vector<std::complex<double>>> DCE_Evolution::eigen2cppType(const std::vector<wvVec > & solutions){

    std::vector<std::vector<std::complex<double>>> retData;
    for(const auto& eigenTypeVec:solutions){
        std::vector<std::complex<double>> oneVec;
        for(int i=0;i<eigenTypeVec.size();i++){
            oneVec.push_back(eigenTypeVec(i));
        }
        retData.push_back(oneVec);
    }

    return retData;


}


///evolution
void DCE_Evolution::evolution(){
    wvVec PsiInit=this->Psi0;
    for(int fls=0;fls<this->flushNum;fls++){
        wvVec PsiFinal=this->evolutionPerFlush(fls,PsiInit);
        PsiInit=PsiFinal;
//        std::cout<<"flush "<<fls<<std::endl;


    }

}


///
/// @param ind index in Psi
/// @return Psi[ind]
double DCE_Evolution::fillOneValInPsi(const int& ind) {
    int n2 = ind % N2;
    int n1 = static_cast<int>((ind - n2) / N2);
    double x1Tmp = this->x1ValsAll[n1];
    double x2Tmp = this->x2ValsAll[n2];
    double valTmp = std::exp(-0.5 * omegac * std::pow(x1Tmp, 2))
            * std::hermite(this->jH1, std::sqrt(omegac) * x1Tmp)
                    * std::exp(-0.5 * omegam * std::pow(x2Tmp, 2))
                    *std::hermite(this->jH2,std::sqrt(omegam)*x2Tmp);
    return valTmp;
}


///
/// @param solutionsPtr unique ptr of solutions per flush
/// @return converted to cpp data type
std::vector<std::vector<std::complex<double>>> DCE_Evolution::eigen2cppType(const std::unique_ptr<std::vector<wvVec >> &solutionsPtr){
    std::vector<std::vector<std::complex<double>>> retData;
    for(const auto& eigenTypeVec:*solutionsPtr){
        std::vector<std::complex<double>> oneVec;
        for(int i=0;i<eigenTypeVec.size();i++){
            oneVec.push_back(eigenTypeVec(i));
        }
        retData.push_back(oneVec);
    }
    return retData;

}
