//
// Created by polya on 2/22/24.
//

export module time_evolution;
#include <complex>
#include <string>


export class matricesHolder{};

export class DCE{

public:
    //parameters that need to be initialized
    int rowNum=-1;//row number of parameters to be read
    int group=-1;//computation group
    int j1H=-1;//initial photon number
    int j2H=-1;//initial phonon number
    double g0=0;//coupling constant
    double omegac=0;//frequency of photon
    double omegam=0;//frequency of phonon
    double omegap=0;//driving frequency
    double er=0;//magnification
    double lmd=0;//constant before quadratic terms of b^{2} and b^{\dagger 2}
public:

    ///
    /// @param csvFileName corresponding to a group
    /// @param r row number
    /// use parseCSV() to initialize the variables
    void parseCSV(const std::string &csvFileName,const int& r);

};

