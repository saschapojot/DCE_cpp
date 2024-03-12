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

    loader.readBinFiles(sortedFiles);
    loader.popolateMatrices();

    loader.solutions = loader.cppType2Eigen();

    std::vector<double> pn = loader.photonAll();
    std::vector<double> ph = loader.phononAll();
//    loader.printVec(pn);
//    loader.printVec(ph);
loader.to_json(pn,ph);



}
