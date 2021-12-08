#include "M_Ind_Mult.hpp"
#include <stdio.h>
#include <iostream>

int main(){
    
    std::string citPat = "citiesMult/Citypatch.txt";
    std::string mobNet = "citiesMult/mobnetwork.txt";
    std::string popAr = "citiesMult/Poparea.txt";
    
    const MobMatrix* M = new MobMatrix(citPat, mobNet, popAr);


    FILE* f = fopen("prueba.txt", "wt");
    
    
    double inf;
    M_Ind_Mult* A = new M_Ind_Mult(*M);


    for(int k = 0; k < 10; k++){
        A->Iteracion(*M);
        inf = 0;
        for(int i = 0; i < M->N; i++){
            //printf("%lf\t%d\n", A->rho[i], M->population[i]);
            inf += A->rho[i] * M->population[i];
        }
        printf("%d\t%lf\n", k, inf);
    }
    fclose(f);
    // for(int i = 0; i < M->N; i++){
    //     printf("%lf\n", A->rho[i]);
    // }


    return 0;
}