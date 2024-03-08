// MyProgram.cpp


#include "timeEvolution.hpp"

int main(int argc, char *argv[]){
    auto evo=DCE_Evolution(6,0);
    evo.populatedMatrices();
    evo.initPsi();

    evo.createOutDir();
    evo.initTimeInds();

    evo.evolution();













}
