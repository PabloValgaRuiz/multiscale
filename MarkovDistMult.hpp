#pragma once

#include <vector>
#include "MobMatrix.hpp"
#include <cmath>

class MarkovDistMult{
public: //private
    double pI;
    double pC;
    double beta;
    double beta0;
    double mu = 0.2;
    double kD = 8, kN = 3; //media de contactos por agente
    double zD, zN; //normalizacion
    std::vector<std::vector<double>> rho, ProbInf;
    std::vector<double> IeffDia, IeffNoche, neff, PD, PN, fvector, sigma;
    
    void calculaIeffneff(const MobMatrix& T){
        neff.clear();
        IeffDia.clear();
        IeffNoche.clear();
        
        neff.resize(T.N, 0);
        IeffDia.resize(T.N, 0);
        IeffNoche.resize(T.N, 0);

        for(int i = 0; i < T.N; i++){
            for(int j = 0; j < T.vecinos[i]; j++){

                if(T.cityPatch[i] == T.cityPatch[T.Mvecinos[i][j]]){

                    neff[T.Mvecinos[i][j]] += pI * T.Mpesos[i][j];
                    IeffDia[T.Mvecinos[i][j]] += pI * T.Mpesos[i][j] * rho[i][j];

                    neff[i] += (1 - pI) * T.Mpesos[i][j];
                    IeffDia[i] += (1 - pI) * T.Mpesos[i][j] * rho[i][j];
                }
                else{
                    neff[T.Mvecinos[i][j]] += pC * T.Mpesos[i][j];
                    IeffDia[T.Mvecinos[i][j]] += pC * T.Mpesos[i][j] * rho[i][j];

                    neff[i] += (1 - pC) * T.Mpesos[i][j];
                    IeffDia[i] += (1 - pC) * T.Mpesos[i][j] * rho[i][j];
                }
            }
        }
        for(int i = 0; i < T.N; i++){
            for(int j = 0; j < T.vecinos[i]; j++){
                IeffNoche[i] += rho[i][j] * T.Mpesos[i][j];
            }
        }
    }

    void calculaP(const MobMatrix& T){
        for(int i = 0; i < T.N; i++){
            if(neff[i] != 0)
                PD[i] = 1 - pow(1 - beta * IeffDia[i]/neff[i], zD * fvector[i]);
            else PD[i] = 0;
            if(T.population[i] != 0)
                PN[i] = 1 - pow(1 - beta * IeffNoche[i]/T.population[i], zN * sigma[i]);
            else PN[i] = 0;
        }
    }

    void calculaProbInf(const MobMatrix& T){
        double prov1, prov2;
        for(int i = 0; i < T.N; i++){
            for(int j = 0; j < T.vecinos[i]; j++){
                if(T.cityPatch[i] == T.cityPatch[T.Mvecinos[i][j]]){
                    ProbInf[i][j] = (1 - pI) * (PD[i] + (1 - PD[i])*PN[i]) + pI * (PD[T.Mvecinos[i][j]] + (1 - PD[T.Mvecinos[i][j]]) * PN[i]);
                }
                else{
                    ProbInf[i][j] = (1 - pC) * (PD[i] + (1 - PD[i])*PN[i]) + pC * (PD[T.Mvecinos[i][j]] + (1 - PD[T.Mvecinos[i][j]]) * PN[i]);
                }
            }
            //std::cout << ProbInf[i] << std::endl;
        }
    }

    double f(double neff, double area){
    #ifdef topecontactos
        if(area != 0){
            if(neff/area > 53000)
                return 53000;
        }
    #endif
        double Kmax = 10;
        double psi = 1.0 / 1200; // 1/<d>

        if(area != 0)
            //return neff/area;
            return Kmax - (Kmax - 1) * exp(-psi * neff / area);
        else return 0;
    }

    void calcFZ(const MobMatrix& T){
        
        // z CONSTANTE CON P=0
        fvector.clear();
        fvector.resize(T.N);
        sigma.clear();
        sigma.resize(T.N);
        for(int i = 0; i < T.N; i++){
            fvector[i] = f(T.population[i], T.area[i]);
            sigma[i] = 1;
        }

        double temp = 0, tempN = 0;
        for(int i = 0; i < T.N; i++){
            temp += T.population[i] * fvector[i];
            tempN += T.population[i] * sigma[i];
        }
        zD = T.Pob * kD / temp;
        zN = T.Pob * kN / tempN;

        for(int i = 0; i < T.N; i++){
            //Calcular fvector para P distinto de cero finalmente
            fvector[i] = f(neff[i], T.area[i]);
        }

        //---------------------------------------------------------

        // z VARIANDO CON P
        // fvector.clear();
        // fvector.resize(T.N);
        // sigma.clear();
        // sigma.resize(T.N);
        // for(int i = 0; i < T.N; i++){
        //     fvector[i] = f(neff[i], T.area[i]);
        //     sigma[i] = 1;
        // }
        // double temp = 0;
        // double tempN = 0;
        // for(int i = 0; i < T.N; i++){
        //     temp += neff[i] * fvector[i];
        //     tempN += T.population[i] * sigma[i];
        // }
        // if(temp != 0) zD = T.Pob * kD / temp;
        // else zD = 0;
        // if(tempN != 0) zN = T.Pob * kN / tempN;
        // else zN = 0;
    }

    void calculaBeta0(const MobMatrix& T){
        double temp = 0, tempN = 0, F = 0, SIGMA = 0;

        for(int i = 0; i < T.N; i++){
            temp += T.population[i] * f(T.population[i], T.area[i]);
        }
        double ztemp = T.Pob * kD / temp;
        double zNtemp = zN;

        temp = tempN = 0;
        int I = 0;
        for(int i = 0; i < T.N; i++){
            temp = f(T.population[i], T.area[i]);
            tempN = sigma[i];
            if(ztemp*temp + zNtemp*tempN > ztemp*F + zNtemp*SIGMA){F = temp; SIGMA = tempN; I = i;}
        }
        beta0 = mu/(F * ztemp + SIGMA * zNtemp);
        //std::cout << I << " " << ztemp << " " << T.population[I] << "   " << F << " " << beta0 << std::endl;
    }

public:

    std::vector<double> infectados;
    double infectadosTotal;
    
    MarkovDistMult(double _pI, double _pC, const MobMatrix& T)
    :pI{_pI}, pC{_pC}{
        rho.resize(T.N);
        ProbInf.resize(T.N);
        for(int i = 0; i < T.N; i++){
            rho[i].resize(T.vecinos[i]);
            ProbInf[i].resize(T.vecinos[i]);
        }
        IeffDia.resize(T.N);
        IeffNoche.resize(T.N);
        neff.resize(T.N);
        PD.resize(T.N);
        PN.resize(T.N);
        calculaIeffneff(T);
        calcFZ(T);
        calculaBeta0(T);
    }

    void inicializar(double _rhoInicial){
        for(int i = 0; i < rho.size(); i++){
            for(int j = 0; j < rho[i].size(); j++){
                rho[i][j] = _rhoInicial;
            }
        }
    }

    void contarInfectados(const MobMatrix& T){
        infectados.clear();
        infectados.resize(T.N, 0);
        infectadosTotal = 0;
        for(int i = 0; i < T.N; i++){
            for(int j = 0; j < T.vecinos[i]; j++){
                infectados[i] += rho[i][j] * T.Mpesos[i][j];
            }
            infectadosTotal += infectados[i];
            if(T.population[i] != 0)
                infectados[i] /= T.population[i];
        }
        if(T.Pob != 0)
            infectadosTotal /= T.Pob;
    }

    void iteracion(const MobMatrix& T){
        calculaIeffneff(T);
        calculaP(T);
        calculaProbInf(T);
        
        for(int i = 0; i < T.N; i++){
            for(int j = 0; j < T.vecinos[i]; j++){
                rho[i][j] = (1 - mu) * rho[i][j] + (1 - rho[i][j]) * ProbInf[i][j];
            }
        }
    }

    //escoger la beta en funcion de beta0
    void setBeta(double _beta){
        beta = _beta * beta0;
        //std::cout << "Beta = " << beta << std::endl;
    }

    //Getters
    const double& get_Beta0() const {return beta0;}
    const double& get_mu() const {return mu;}
    const double& get_pI() const {return pI;}
    const double& get_pC() const {return pC;}
    const std::vector<double>& get_fvector() const {return fvector;}
    const std::vector<double>& get_sigma() const {return sigma;}
    const double& get_zD() const {return zD;}
    const double& get_zN() const {return zN;}
    const std::vector<double>& get_neff() const {return neff;}

};