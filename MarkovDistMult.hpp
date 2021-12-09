#pragma once

#include <vector>
#include "MobMatrix.hpp"
#include <cmath>
#include <stdexcept>

class MarkovDistMult{
public: //private
    double pI;
    double pC;
    double beta;
    double beta0;
    double mu = 0.2;
    double kW = 8, kH = 3, kN = 3; //media de contactos por agente
    double zW, zH, zN; //normalizacion
    std::vector<std::vector<double>> rho, ProbInf;
    std::vector<double> IeffW, IeffH, IeffNoche, neffW, neffH, PW, PH, PN, fWvector, fHvector, sigma;
    
    void calculaIeffneff(const MobMatrix& T){
        neffW.clear();
        neffH.clear();
        IeffW.clear();
        IeffH.clear();
        IeffNoche.clear();
        
        neffW.resize(T.N, 0);
        neffH.resize(T.N, 0);
        IeffW.resize(T.N, 0);
        IeffH.resize(T.N,0);
        IeffNoche.resize(T.N, 0);

        for(int i = 0; i < T.N; i++){
            for(int j = 0; j < T.vecinos[i]; j++){

                if(T.cityPatch[i] == T.cityPatch[T.Mvecinos[i][j]]){

                    neffW[T.Mvecinos[i][j]] += pI * T.Mpesos[i][j];
                    IeffW[T.Mvecinos[i][j]] += pI * T.Mpesos[i][j] * rho[i][j];

                    neffH[i] += (1 - pI) * T.Mpesos[i][j];
                    IeffH[i] += (1 - pI) * T.Mpesos[i][j] * rho[i][j];
                }
                else{
                    neffW[T.Mvecinos[i][j]] += pC * T.Mpesos[i][j];
                    IeffW[T.Mvecinos[i][j]] += pC * T.Mpesos[i][j] * rho[i][j];

                    neffH[i] += (1 - pC) * T.Mpesos[i][j];
                    IeffH[i] += (1 - pC) * T.Mpesos[i][j] * rho[i][j];
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
            if(neffW[i] != 0)
                PW[i] = 1 - pow(1 - beta * IeffW[i]/neffW[i], zW * fWvector[i]);
            else PW[i] = 0;

            if(neffH[i] != 0)
                PH[i] = 1 - pow(1 - beta * IeffH[i]/neffH[i], zH * fHvector[i]);
            else PH[i] = 0;

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
                    ProbInf[i][j] = (1 - pI) * (PH[i] + (1 - PH[i])*PN[i]) + pI * (PW[T.Mvecinos[i][j]] + (1 - PW[T.Mvecinos[i][j]]) * PN[i]);
                }
                else{
                    ProbInf[i][j] = (1 - pC) * (PH[i] + (1 - PH[i])*PN[i]) + pC * (PW[T.Mvecinos[i][j]] + (1 - PW[T.Mvecinos[i][j]]) * PN[i]);
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
        fWvector.clear();
        fHvector.clear();
        fWvector.resize(T.N);
        fHvector.resize(T.N);
        sigma.clear();
        sigma.resize(T.N);

        //p=0 for home contacts, p=1 for work contacts
        std::vector<double> neffWtemp, neffHtemp;
        neffWtemp.resize(T.N, 0);   neffHtemp.resize(T.N, 0);
        for(int i = 0; i < T.N; i++){
            for(int j = 0; j < T.vecinos[i]; j++){
                neffWtemp[T.Mvecinos[i][j]] += T.Mpesos[i][j];
                neffHtemp[i] += T.Mpesos[i][j];
            }
        }

        for(int i = 0; i < T.N; i++){
            fWvector[i] = f(neffWtemp[i], T.area[i]);
            fHvector[i] = f(neffHtemp[i], T.area[i]);
            sigma[i] = 1;
        }

        double tempW = 0, tempH = 0, tempN = 0;
        for(int i = 0; i < T.N; i++){
            tempW += neffWtemp[i] * fWvector[i];
            tempH += neffHtemp[i] * fHvector[i];
            tempN += T.population[i] * sigma[i];
        }

        if(tempW == 0 || tempH == 0 || tempN == 0)
            throw std::runtime_error("zW, zH or zN = infinity");
        zW = T.Pob * kW / tempW;
        zH = T.Pob * kH / tempH;
        zN = T.Pob * kN / tempN;

        for(int i = 0; i < T.N; i++){
            //Calcular fvector para p distinto de cero finalmente
            fWvector[i] = f(neffW[i], T.area[i]);
            fHvector[i] = f(neffH[i], T.area[i]);
            sigma[i] = 1;
        }

    }

    void calculaBeta0(const MobMatrix& T){

        //In the redefined mobility model, there is no people at work on p = 0, so it's always
        //the people at home that have contacts

        double tempH = 0, tempN = 0, FH = 0, SIGMA = 0;

        for(int i = 0; i < T.N; i++){
            tempH += T.population[i] * f(T.population[i], T.area[i]);
        }
        double zHtemp = 0;
        if(tempH != 0)
            zHtemp = T.Pob * kH / tempH;
        double zNtemp = zN;

        tempH = tempN = 0;
        int I = 0;
        for(int i = 0; i < T.N; i++){
            tempH = f(T.population[i], T.area[i]);
            tempN = sigma[i];
            if(zHtemp*tempH + zNtemp*tempN > zHtemp*FH + zNtemp*SIGMA){
                FH = tempH; SIGMA = tempN; I = i;
            }
        }
        if(FH * zHtemp + SIGMA * zNtemp != 0)
            beta0 = mu/(FH * zHtemp + SIGMA * zNtemp);
        else
            throw std::runtime_error("Beta_0 = infinity");
        //std::cout << I << "\t" << zHtemp << "\t" << FH << "\t" << zNtemp << "\t" << SIGMA << "\t" << zH << "\t" <<  beta0 << std::endl;
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
        IeffW.resize(T.N);
        IeffH.resize(T.N);
        IeffNoche.resize(T.N);
        neffW.resize(T.N);
        neffH.resize(T.N);
        PW.resize(T.N);
        PH.resize(T.N);
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
    const std::vector<double>& get_fWvector() const {return fWvector;}
    const std::vector<double>& get_fHvector() const {return fHvector;}
    const std::vector<double>& get_sigma() const {return sigma;}
    const double& get_zW() const {return zW;}
    const double& get_zH() const {return zH;}
    const double& get_zN() const {return zN;}
    const std::vector<double>& get_neffW() const {return neffW;}
    const std::vector<double>& get_neffH() const {return neffH;}
};