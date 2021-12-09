#pragma once

#include <vector>
#include "MobMatrix.hpp"
#include <cmath>
#include <random>
#include <stdexcept>

enum Compartment {S, I};
enum ContactPlace {H, W};
struct Individual{
    Compartment state{S};
    int Org{};
    int Des{};
    int Desplazamiento{};
    ContactPlace contactPlace{};

    //Constructor de individuos
    Individual(Compartment _state, size_t _Org, size_t _Des, size_t _Desplazamiento, ContactPlace _contactPlace)
    :state{_state}, Org{_Org}, Des{_Des}, Desplazamiento{_Desplazamiento}, contactPlace{_contactPlace}{}
};

class MCDistMult{
public: //private
    double pI;
    double pC;
    double beta;
    double beta0;
    double mu = 0.2;
    double kW = 8, kH = 3, kN = 3; //media de contactos por agente
    double zW, zH, zN; //normalizacion

    std::vector<double> fWvector, fHvector, sigma;
    std::vector<double> PW, PH, PN;
    std::vector <int> IeffW, IeffH, IeffNoche, neffW, neffH;
    std::vector<Individual> individuals;

    //Calcular la probabilidad de infectarse al contactar con un agente
    //en el lugar de destino, sin tener aun en cuenta la cantidad de
    //agentes con la que interactua (el exponente f*z)
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

    //Calcular la cantidad total de infectados tanto en el destino
    //como en el origen
    void calculaIeffneff(){
        for(int i = 0; i < IeffW.size(); i++){
            IeffW[i] = IeffH[i] = IeffNoche[i] = neffW[i] = neffH[i] = 0;
        }
        for(const Individual& k : individuals){
            if(k.state == I){
                if(k.contactPlace == W)
                    IeffW[k.Desplazamiento]++;
                else
                    IeffH[k.Desplazamiento]++;

                IeffNoche[k.Org]++;
            }
            if(k.contactPlace == W)
                neffW[k.Desplazamiento]++;
            else
                neffH[k.Desplazamiento]++;
        }
    }

    //Calcular de forma aleatoria si los individuos se desplazan
    //o se quedan en su lugar de origen. En caso de desplazarse
    //solamente tienen un posible lugar de destino
    void desplaza(const MobMatrix& T){
        std::random_device rd;
        std::mt19937 mt(rd());
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        for(Individual& k : individuals){
            if(T.cityPatch[k.Org] == T.cityPatch[k.Des]){
                if (dist(mt) < pI){		                //SI se desplaza
                    k.Desplazamiento = k.Des;
                    k.contactPlace = W;    
                }
                else{									//NO se desplaza
                    k.Desplazamiento = k.Org;
                    k.contactPlace = H;    
                }
            }
            else{
                if (dist(mt) < pC){		                //SI se desplaza
                    k.Desplazamiento = k.Des;
                    k.contactPlace = W;    
                }
                else{									//NO se desplaza
                    k.Desplazamiento = k.Org;
                    k.contactPlace = H;
                }
            }
        }
    }

    //Devuelve la funcion f para una poblacion y un area dados
    //Lento: mejor usar el vector fvector, precalculado para cada nodo
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

    //Calcular inicialmente fvector y las funciones de normalizacion zD y zN
    void calcFZ(const MobMatrix& T){
        // z CONSTANTE CON P=0
        fWvector.clear();
        fHvector.clear();
        fWvector.resize(T.N);
        fHvector.resize(T.N);
        sigma.clear();
        sigma.resize(T.N);

        //p=0 for home contacts, p=1 for work contacts
        std::vector<int> neffWtemp, neffHtemp;
        neffWtemp.resize(T.N, 0);   neffHtemp.resize(T.N, 0);
        for(Individual k : individuals){
            neffWtemp[k.Des]++;
            neffHtemp[k.Org]++;
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

        //Calcular fvector para p distinto de cero finalmente
        for(int i = 0; i < T.N; i++){
            fWvector[i] = f(neffW[i], T.area[i]);
            fHvector[i] = f(neffH[i], T.area[i]);
            sigma[i] = 1;
        }
    }

    //Calcular el umbral epidemico para pI=pC=0
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
    }


public:
    //Fraccion de infectados en cada nodo de origen
    //El algoritmo no utiliza esta variable, es solo para su visibilidad publica
    std::vector<double> infectados;

    //Fraccion total de infectados
    double infectadosTotal;

    //Métodos públicos para realizar iteraciones montecarlo. En general, utilizar en el orden que se muestran.
    //El método de iteración llama al de update de forma implícita, por lo que no hay que llamarlo cada vez
    //El método de contar infectados calcula Ieff y neff de nuevo, asi que no hace falta llamar update o
    //calculaIeffneff

    //Constructor del objeto montecarlo
    MCDistMult(double _pI, double _pC, const MobMatrix& T) :pI{_pI}, pC{_pC}{
        individuals.clear();
        for(int i = 0; i < T.N; i++){
            for(int j = 0; j < T.vecinos[i]; j++){
                for(int k = 0; k < T.Mpesos[i][j]; k++){
                    individuals.emplace_back(S, i, T.Mvecinos[i][j], i, H);
                }
            }
        }
        if(individuals.size() != T.Pob) throw std::logic_error{"MCDistMult: Population doesn't match."};

        IeffW.resize(T.N);
        IeffH.resize(T.N);
        IeffNoche.resize(T.N);
        neffW.resize(T.N);
        neffH.resize(T.N);
        PW.resize(T.N);
        PH.resize(T.N);
        PN.resize(T.N);
        calculaIeffneff();
        calcFZ(T);
        calculaBeta0(T);
    }

    //Establecer un numero inicial de infectados
    void inicializar(double _rhoInicial){
        std::random_device rd;
        std::mt19937 mt(rd());
        std::uniform_real_distribution<double> dist(0.0, 1.0);

        for(Individual& k : individuals){
            if(dist(mt) < _rhoInicial)
                k.state = I;
            else k.state = S;
        }
    }

    //escoger la beta en funcion de beta0
    void setBeta(double _beta){
        beta = _beta * beta0;
        //std::cout << "Beta = " << beta << std::endl;
    }

    //Avanzar un paso en el tiempo en la simulacion
    void iteracion(const MobMatrix& T){
        //Si zD * fi = 8.2, hay que hacer 8 contactos, y luego una probabilidad de 0.2 de tener un 9no contacto
        
        //Generador de numeros aleatorios
        std::random_device rd;
        std::mt19937 mt(rd());
        std::uniform_real_distribution<double> dist(0.0, 1.0);

        update(T);

        for(Individual& k : individuals){
            switch(k.state){
                case S:
                    // //Truncamiento de los contactos
                    // int contactosD = static_cast<int>(zD * fvector[k.Desplazamiento]);
                    // int contactosN = static_cast<int>(zN * sigma[k.Org]);
                    
                    // //Parte decimal de los contactos
                    // double probContactoD = zD * fvector[k.Desplazamiento] - static_cast<double>(contactosD);
                    // double probContactoN = zN * sigma[k.Org] - static_cast<double>(contactosN);

                    // //Probabilidad final de infectarse
                    // double probInfIndiv;
                    // if(dist(mt) < probContactoD)
                    //     probInfIndiv = 1 - pow(PD[k.Desplazamiento], contactosD + 1); //9 contactos
                    // else probInfIndiv = 1 - pow(PD[k.Desplazamiento], contactosD);     //8 contactos

                    // //Contagios por el dia
                    // if(dist(mt) < probInfIndiv)
                    //     k.state = I;
                    // else{   //Contagios por la noche
                    //     if(dist(mt) < probContactoN)
                    //         probInfIndiv = 1 - pow(PN[k.Org], contactosN + 1);
                    //     else probInfIndiv = 1 - pow(PN[k.Org], contactosN);
                    //     if(dist(mt) < probInfIndiv){
                    //         k.state = I;
                    //     }
                    // }

                    if(k.contactPlace == W){
                        if(dist(mt) < PW[k.Desplazamiento]){
                            k.state = I;
                        }
                        else{
                            if(dist(mt) < PN[k.Org]){
                                k.state = I;
                        }
                    }
                    }
                    else{
                        if(dist(mt) < PH[k.Desplazamiento]){
                            k.state = I;
                        }
                        else{
                            if(dist(mt) < PN[k.Org]){
                                k.state = I;
                        }
                    }
                    }

                    
                    
                    break;
                case I:
                    if(dist(mt) < mu)
                        k.state = S;
                    break;
            }
        }


    }

    //Contar infectados, recalcular f,sigma, calcular z's
    void update(const MobMatrix& T){
        desplaza(T);
        calculaIeffneff();
        calcFZ(T);
        calculaP(T);
    }

    //Contar la fraccion de infectados en cada nodo de origen y la total
    void calcInfectados(const MobMatrix& T){
        calculaIeffneff();

        infectados.clear();
        infectados.resize(T.N, 0);
        infectadosTotal = 0;
        for(int i = 0; i < T.N; i++){
            infectados[i] = IeffNoche[i];
            infectadosTotal += infectados[i];
            if(T.population[i] != 0)
                infectados[i] /= T.population[i];
        }
        if(T.Pob != 0)
            infectadosTotal /= T.Pob;
    }

    //GETTERS
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
    const std::vector<int>& get_neffW() const {return neffW;}
    const std::vector<int>& get_neffH() const {return neffH;}
};

