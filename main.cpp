// MyProgram.cpp


#include "timeEvolution.hpp"

int main(int argc, char *argv[]){
    auto evo=DCE_Evolution(6,0);
    evo.populatedMatrices();
    evo.initPsi();

    evo.createOutDir();
    evo.initTimeInds();

    int j=0;
    DCE_Evolution::wvVec psi1=evo.oneStepEvolution(j,evo.Psi0);
    std::cout<<psi1<<std::endl;














}
