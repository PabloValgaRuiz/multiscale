#pragma once

#include <stdlib.h>
#include <vector>
#include <string>
#include <stdio.h>
#include <iostream>
#include <fstream>

class MobMatrix{
private:

    std::string city_patch;
    std::string mobility_network;
    std::string pop_area;

public:
    int N = 0, Nc = 0, Pob = 0;
    double Ratio = 0.2;
    long long Links = 0;
    std::vector<int> cityPatch;
    std::vector<int> population;
    std::vector<double> area;
    std::vector<int> vecinos;
    std::vector<int> vecinosT;
    std::vector<std::vector<int>> Mvecinos;
    std::vector<std::vector<int>> MvecinosT;
    std::vector<std::vector<double>> Mpesos;
    std::vector<std::vector<double>> MpesosT;
    std::vector<double> mC; //Matriz que cuenta la cantidad de gente de cada nodo que sale de su CIUDAD
    std::vector<double> mI; //Matriz que cuenta la cantidad de gente de cada nodo que no sale de su ciudad, pero sí de su NODO
    std::vector<std::vector<double>> Rin, RinT;
    std::vector<std::vector<double>> Rout, RoutT;

    MobMatrix(const std::string& _city_patch, const std::string& _mobility_network, const std::string& _pop_area){
        this->city_patch = _city_patch;
        this->mobility_network = _mobility_network;
        this->pop_area = _pop_area;
        readNCityPatch();
        readPopArea();
        leer_vecinos();
        leer_matrices();
        calculaTraspuesta();
    }

private:
    void readNCityPatch(){
        N = 0; Nc = 0;
        this->cityPatch.resize(0);

        std::ifstream inFile;
        inFile.open(this->city_patch);
        if (!inFile) {
            std::cout << "Unable to open file";
            std::exit(1); // terminate with error
        }

        int patch, city;
        while(inFile >> patch >> city){
            this->cityPatch.push_back(city);
            if(patch > N)
                N = patch;
            if(city > Nc)
                Nc = city;
        }
        N++; Nc++;
        std::cout << N << " " << Nc << std::endl;
        inFile.close();
    }

    void readPopArea(){
        Pob = 0;
        this->population.resize(N);
        this->area.resize(N);
        
        std::ifstream inFile;
        inFile.open(this->pop_area);
        if (!inFile) {
            std::cout << "Unable to open file";
            std::exit(1); // terminate with error
        }
        int i;
        int k = 0;
        while(inFile >> i >> population[i] >> area[i]){
            population[i] *= Ratio;
            Pob += population[i];
            //std::cout << population[i] << std::endl;
            if(i != k++) {
                std::cerr << "Fichero de áreas incompleto" << std::endl;
                }
        }
        inFile.close();
    }

    void leer_vecinos()
    {
        int i, j1, j2, trash1, trash2, trash;
        vecinos.resize(N);
        vecinosT.resize(N);
        for ( i = 0 ; i < N ; i++)
            vecinos[i] = 0;
        
        std::ifstream inFile;
        inFile.open(this->mobility_network);
        if (!inFile) {
            std::cout << "Unable to open file";
            std::exit(1); // terminate with error
        }
        Links = 0;
        while(inFile >> j1 >> j2 >> trash){
            vecinos.at(j1)++;
            vecinosT.at(j2)++;
            Links++;
        }
        inFile.close();
    }

    void leer_matrices()
    {
        Mvecinos.resize(N);
        Mpesos.resize(N);
        for(int i = 0; i < N; i++){
            Mvecinos[i].resize(vecinos[i]);
            Mpesos[i].resize(vecinos[i]);
        }
        std::ifstream inFile;
        inFile.open(this->mobility_network);
        if (!inFile){
            std::cout << "Unable to open file";
            std::exit(1); // terminate with error
        }
        int I,i,j,trash;
        int temp;

        //Recuento de la poblacion: los datos de pesos y poblaciones no dan una matriz Mpesos de numeros enteros
        //por lo que hay que truncarlos y recalcular las poblaciones
        Pob = 0;
        int population_temp;
        for(i = 0 ; i < N ; i++)
        {   
            temp = 0;
            for (j = 0 ; j < vecinos[i] ; j++)
            {
                inFile >> trash >> Mvecinos[i][j] >> Mpesos[i][j];
                temp += Mpesos[i][j];
            }
            population_temp = 0;
            for(j = 0; j < vecinos[i]; j++){
                Mpesos[i][j] = static_cast<int>(population[i] * Mpesos[i][j] / temp); //Mpesos/temp es equivalente a R_ij -> ahora Mpesos tiene informacion de poblacion total
                population_temp += Mpesos[i][j];
            }
            population[i] = population_temp;
            Pob += population_temp;
        }
        inFile.close();
    }

    void calculaTraspuesta()
    {  
        MvecinosT.resize(N);
        MpesosT.resize(N);
        for(int i = 0; i < N; i++){
            MvecinosT[i].resize(vecinosT[i]);
            MpesosT[i].resize(vecinosT[i]);
        }
        int k = 0;
        int i = 0, j = 0;
        int B;
        std::vector<int> AUX; AUX.resize(N);

        for (i = 0; i < N; i++){
            for (j = 0; j < vecinos[i]; j++){
                B = Mvecinos[i][j];
                MvecinosT[B][AUX[B]] = i;
                MpesosT[B][AUX[B]] = Mpesos[i][j];
                AUX[B] = AUX[B] + 1;
            }
        }
    }
};