//
// Created by polya on 2/26/24.
//

#ifndef DCE_CPP_TIMEEVOLUTION_HPP
#define DCE_CPP_TIMEEVOLUTION_HPP
#include <string>
#include <iostream>
#include <memory>
#include <cstdio>
#include <array>
#include <regex>



class DCE_Evolution {
public:
    DCE_Evolution(const int &group, const int &row) {
        this->groupNum = group;
        this->rowNum = row;

        this->parseCSV(group, row);
    }

public:
    int jH1 = -1;
    int jH2 = -1;
    double g0 = 0;
    double omegam = 0;
    double omegap=0;
    double omegac = 0;
    double er = 0;
    double thetaCoef = 0;
    int groupNum = -1;
    int rowNum = -1;

    ///
    /// @param group group number
    /// @param row row number
    ///parse csv with group number group
    void parseCSV(const int &group, const int &row);

    ///
    /// @param cmd python execution string
    /// @return signal from the python
    static std::string execPython(const char *cmd);

};


#endif //DCE_CPP_TIMEEVOLUTION_HPP
