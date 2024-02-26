//
// Created by polya on 2/26/24.
//

#ifndef DCE_CPP_TIMEEVOLUTION_HPP
#define DCE_CPP_TIMEEVOLUTION_HPP


class DCE_Evolution{
public:
    DCE_Evolution(const int& group,const int&row){
        this->groupNum=group;
        this->rowNum=row;

        this->parseCSV(group,row);
    }

public:
int jH1=-1;
int jH2=-1;
double g0=0;
double omegam=0;
double omegac=0;
double er=0;
double thetaCoef=0;
int groupNum=-1;
int rowNum=-1;


void parseCSV(const int &group, const int&row);


};


#endif //DCE_CPP_TIMEEVOLUTION_HPP
