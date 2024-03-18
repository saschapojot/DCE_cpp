// MyProgram.cpp


#include "timeEvolution.hpp"

int main(int argc, char *argv[]){
    if(argc!=3){
        std::cerr<<"wrong number of arguments"<<std::endl;
        exit(1);
    }
    int groupNum=std::stoi(argv[1]);
    int rowNum=std::stoi(argv[2]);
    auto evo=DCE_Evolution(groupNum,rowNum);
    evo.populatedMatrices();
    evo.initPsiSerial();

    evo.createOutDir();
    evo.initTimeInds();

    evo.evolution();













}
