#pragma once

#include <vector>
#include "MobMatrix.hpp"

class M_Ind_Mult{
public:
    double pC = 0.2;
    double pI = 0.1;
    double lambda = 0.01;
    double mu = 0.2;

    std::vector<double> rho;
    std::vector<double> Ieff, neff;
    std::vector<double> P, ProbInf;



    M_Ind_Mult(const MobMatrix& T);
    void Iteracion(const MobMatrix& T);
    void calculaIeffneff(const MobMatrix& T);
    double f(double neff, double area, double eps);
    double z(const MobMatrix& T);
};