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
    const auto tEvoStart{std::chrono::steady_clock::now()};
    evo.evolution();
    const auto tEvoEnd{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_secondsAll{tEvoEnd - tEvoStart};
    std::cout<<"Evolution time: "<< elapsed_secondsAll.count() / 3600.0 << " h" << std::endl;













}
