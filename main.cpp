// MyProgram.cpp


#include "timeEvolution.hpp"

int main(int argc, char *argv[]){
    auto evo=DCE_Evolution(6,0);
    evo.populatedMatrices();
    evo.initPsi();
    evo.createOutDir();
    evo.initTimeInds();

    for(int fls=0;fls<evo.flushNum;fls++){
        evo.evolutionPerFlush(fls, DCE_Evolution::wvVec() );
    }










}
