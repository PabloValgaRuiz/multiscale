#include "M_Ind_Mult.hpp"
#include <stdio.h>
#include <iostream>
#include <math.h>


M_Ind_Mult::M_Ind_Mult(const MobMatrix& T){
    rho.resize(T.N);
    Ieff.resize(T.N);
    neff.resize(T.N);
    P.resize(T.N);
    ProbInf.resize(T.N);

    for(int i = 0; i < T.N; i++)
        rho[i] = 0.001;
}

void M_Ind_Mult::calculaIeffneff(const MobMatrix& T){
    double provneff, provIeff;
    for(int i = 0; i < T.N; i++){
        Ieff[i] = 0;
        neff[i] = 0;
        for(int j = 0; j < T.vecinos[i]; j++){
            provneff = T.population[T.Mvecinos[i][j]];
            if(T.cityPatch[i] == T.cityPatch[T.Mvecinos[i][j]]){
                provneff *= (1 - pC * T.mC[T.Mvecinos[i][j]]);
                if(i == T.Mvecinos[i][j]){
                    provneff *= (1 - pI * T.mI[T.Mvecinos[i][j]] );
                    
                }
                else{
                    provneff *= (pI * T.mI[T.Mvecinos[i][j]] * T.RinT[i][j]); //Rin porque es en la misma ciudad
                    
                }
            }
            else{
                provneff *= pC * T.mC[T.Mvecinos[i][j]] * T.RoutT[i][j];
                
            }
            //if(i == 299 && j == 18)
                //std::cout << T.Mvecinos[i][j] << "  " << T.RoutT[i][j] << std::endl;
            provIeff = provneff;
            provIeff *= rho[T.Mvecinos[i][j]];

            Ieff[i] += provIeff;
            neff[i] += provneff;
            

        }

        //printf("%lf\t%lf\n", Ieff[i], neff[i]);

    }
}

void M_Ind_Mult::Iteracion(const MobMatrix& T){
    calculaIeffneff(T);

    
    for(int i = 0; i < T.N; i++){
        P[i] = 1 - pow(1 - lambda * Ieff[i]/neff[i], z(T) * f(neff[i], T.area[i], 1.0));
    }

    double prov1, prov2;
    for(int i = 0; i < T.N; i++){

        prov1 = prov2 = 0;
        for(int j = 0; j < T.vecinos[i]; j++){
            if(T.cityPatch[i] == T.cityPatch[T.Mvecinos[i][j]] && i != T.Mvecinos[i][j])
                prov1 += T.Rin[i][j] * P[T.Mvecinos[i][j]];

            if(T.cityPatch[i] != T.cityPatch[T.Mvecinos[i][j]])
                prov2 += T.Rout[i][j] * P[T.Mvecinos[i][j]];
        }


        ProbInf[i] = (1 - pC * T.mC[i]) * (  (1 - pI * T.mI[i]) * P[i] + pI * T.mI[i] * prov1  ) + pC * T.mC[i] * prov2;
        //std::cout << ProbInf[i] << std::endl;
    }

    for(int i = 0; i < T.N; i++){
        rho[i] = (1 - mu) * rho[i] + (1 - rho[i]) * ProbInf[i];
    }

    

}

double M_Ind_Mult::f(double neff, double area, double eps){
    //std::cout << neff << std::endl;
    return neff/area;
}

double M_Ind_Mult::z(const MobMatrix& T){
    return 4;
}