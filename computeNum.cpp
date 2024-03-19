//
// Created by polya on 3/12/24.
//
#include "combineSegments.hpp"

int main(int argc, char *argv[]) {
if(argc!=3){
    std::cout << "wrong arguments" << std::endl;
    std::exit(2);
}
    int group = std::stoi(argv[1]);
    int row = std::stoi(argv[2]);
    auto loader = combineSegments(group, row);
    auto strVec = loader.listFiles(loader.outDir);
    auto sortedFiles = loader.sortFiles(strVec);
    loader.catchParameters(sortedFiles[0]);
    loader.initTimeInds();


//    loader.readBinFiles(sortedFiles);
    loader.popolateMatrices();

    std::vector<std::vector<double>> photonNumFromAllFlushes;
    std::vector<std::vector<double>> phononNumFromAllFlushes;
    const auto tStart{std::chrono::steady_clock::now()};
    for (int j = 0; j < sortedFiles.size(); j++) {
        if(j%20==0) {
            std::cout<<"flush "<<j<<std::endl;
        }
        std::string name = sortedFiles[j];

        std::vector<std::vector<std::complex<double>>> inVec = loader.readOneBinFile(name);
        for(const auto& vec:inVec){
            std::cout<<"norm is "<<loader.checkNorm(vec)<<std::endl;
        }
        std::string suffix = "N1" + std::to_string(loader.N1) + "N2" + std::to_string(loader.N2)
                             + "L1" + std::to_string(loader.L1) + "L2" + std::to_string(loader.L2);

        if (j == 0) {
            std::string outInitName = loader.outDir + suffix + "init.xml";
            loader.writeInitLast(outInitName, inVec[0]);
        }
        if (j == sortedFiles.size() - 1) {
            std::string outLastName = loader.outDir + suffix + "final.xml";
            int lengthTmp = inVec.size();
            loader.writeInitLast(outLastName, inVec[lengthTmp - 1]);
        }
        std::vector<wvVec> solutionsInOneFlush = loader.cppType2EigenOneFile(inVec);


        std::vector<double> photonsPerFlush = loader.photonPerFlushSerial(solutionsInOneFlush);


        std::vector<double> phononsPerFlush = loader.phononPerFlushSerial(solutionsInOneFlush);


        photonNumFromAllFlushes.push_back(photonsPerFlush);
        phononNumFromAllFlushes.push_back(phononsPerFlush);


        std::vector<double> photonVals = loader.removeHeadTail(photonNumFromAllFlushes);
        std::vector<double> phononVals = loader.removeHeadTail(phononNumFromAllFlushes);

        loader.to_json(photonVals, phononVals);


    }
    const auto tEnd{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_secondsAll{tEnd - tStart};
    std::cout << "All flushes time: " << elapsed_secondsAll.count() / 3600.0 << " h" << std::endl;
}




