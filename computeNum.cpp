//
// Created by polya on 3/12/24.
//
#include "combineSegments.hpp"

int main(int argc, char *argv[]) {
//if(argc!=3){
//    std::cout << "wrong arguments" << std::endl;
//    std::exit(2);
//}
    int group = 6;
    int row = 0;
    auto loader = combineSegments(group, row);
    auto strVec = loader.listFiles(loader.outDir);
    auto sortedFiles = loader.sortFiles(strVec);
    loader.catchParameters(sortedFiles[0]);
    loader.initTimeInds();


//    loader.readBinFiles(sortedFiles);
    loader.popolateMatrices();

    std::vector<std::vector<double>> photonNumFromAllFlushes;
    std::vector<std::vector<double>> phononNumFromAllFlushes;
    for(const auto&name:sortedFiles){
        const auto tOneFlushStart{std::chrono::steady_clock::now()};
        std::vector<std::vector<std::complex<double>>> inVec=loader.readOneBinFile(name);

        std::vector<wvVec > solutionsInOneFlush=loader.cppType2EigenOneFile(inVec);



        std::vector<double> photonsPerFlush=loader.photonPerFlushSerial(solutionsInOneFlush);


        std::vector<double> phononsPerFlush=loader.phononPerFlushSerial(solutionsInOneFlush);



        photonNumFromAllFlushes.push_back(photonsPerFlush);
        phononNumFromAllFlushes.push_back(phononsPerFlush);
        const auto tOneFlushEnd{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_secondsAll{tOneFlushEnd - tOneFlushStart};
    std::cout << "One flush time: " << elapsed_secondsAll.count() / 3600.0 << " h" << std::endl;
    }

    std::vector<double> photonVals=loader.removeHeadTail(photonNumFromAllFlushes);
    std::cout<<photonVals.size()<<std::endl;
//    loader.solutions = loader.cppType2Eigen();

//    std::vector<double> pn = loader.photonAllSerial();
//    std::vector<double> ph = loader.phononAllSerial();
//    loader.printVec(pn);
//    loader.printVec(ph);
//loader.to_json(pn,ph);



}


