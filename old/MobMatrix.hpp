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

    void readNCityPatch();    
    void readPopArea();
    void leer_vecinos(); //LOS FICHEROS NO PUEDEN ACABAR CON UNA LINEA EN BLANCO: DARÁ UN VECINO DE MÁS Y CRASHEARÁ
    void leer_vecinosT();
    void leer_matrices();
    void calculaTraspuesta();
    void setmC();
    void setmI();
    void setR();
    void calculaRTraspuesta();

public:
    int N = 0, Nc = 0, Links = 0, Pob = 0; Ratio = 1;
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

    MobMatrix(const std::string& _city_patch, const std::string& _mobility_network, const std::string& _pop_area);
};

MobMatrix::MobMatrix(const std::string& _city_patch, const std::string& _mobility_network, const std::string& _pop_area){
    this->city_patch = _city_patch;
    this->mobility_network = _mobility_network;
    this->pop_area = _pop_area;
    readNCityPatch();
    readPopArea();
    leer_vecinos();
    leer_matrices();
    calculaTraspuesta();
    setmC();
    setmI();
    setR();
    calculaRTraspuesta();
}

void MobMatrix::readNCityPatch(){
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

void MobMatrix::readPopArea(){
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
        if(i != k++) {
            std::cerr << "Fichero de áreas incompleto" << std::endl;
            }
    }
    inFile.close();
}

void MobMatrix::leer_vecinos()
{
    int i, j1, j2, trash1, trash2, trash;
    vecinos.resize(N);
    for ( i = 0 ; i < N ; i++)
        vecinos[i] = 0;
    
    std::ifstream inFile;
    inFile.open(this->mobility_network);
    if (!inFile) {
        std::cout << "Unable to open file";
        std::exit(1); // terminate with error
    }

    inFile >> trash >> trash;
    while(inFile >> j1 >> j2 >> trash){
        vecinos[j1]++;
        vecinosT[j2]++;
    }
    inFile.close();
}

void MobMatrix::leer_matrices()
{
    Mvecinos.resize(N);
    Mpesos.resize(N);
    for(int i = 0; i < N; i++){
        Mvecinos[i].resize(vecinos[i]);
        Mpesos[i].resize(vecinos[i]);
    }
    std::ifstream inFile;
    inFile.open(this->mobility_network);
    if (!inFile) {
        std::cout << "Unable to open file";
        std::exit(1); // terminate with error
    }
    int I,i,j,trash;
    int temp;
    for(i = 0 ; i < N ; i++)
    {   
        temp = 0;
        for ( j = 0 ; j < vecinos[i] ; j++)
        {
            inFile >> trash >> Mvecinos[i][j] >> Mpesos[i][j];
            temp += Mpesos[i][j];
        }
        for(j = 0; j < vecinos[i]; j++){
            Mpesos[i][j] = population[i] * Mpesos[i][j] / temp; //Mpesos/temp es equivalente a R_ij -> ahora Mpesos tiene informacion de poblacion total
        }
    }
    inFile.close();
}

void MobMatrix::calculaTraspuesta()
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
	int* AUX = (int*)malloc(N * sizeof(int));
	for (i = 0; i < N; i++)
		AUX[i] = 0;

	for (i = 0; i < N; i++){
		for (j = 0; j < vecinos[i]; j++){
            B = Mvecinos[i][j];
            MvecinosT[B][AUX[B]] = i;
		    MpesosT[B][AUX[B]] = Mpesos[i][j];
		    AUX[B] = AUX[B] + 1;
		}
	}
	free(AUX);
}

void MobMatrix::setmC(){
    mC.resize(N);
    double prov1, prov2;
    for(int i = 0; i < N; i++){
        prov1 = prov2 = 0;
        for(int j = 0; j < vecinos[i]; j++){
            prov2 += Mpesos[i][j];              //prov2 -> toda la movilidad del nodo, incluidos los autoloops (?)
            if(cityPatch[Mvecinos[i][j]] != cityPatch[i])
                prov1 += Mpesos[i][j];          //prov1 -> solo los que salen de su ciudad

        }
        if(prov2 != 0)
            mC[i] = prov1/prov2; 
        else
            mC[i] = 0;
        //printf("%lf\n", mC[i]);
    }
}

void MobMatrix::setmI(){
    mI.resize(N);
    double prov1, prov2;
    for(int i = 0; i < N; i++){
        prov1 = prov2 = 0;
        for(int j = 0; j < vecinos[i]; j++){
            if(cityPatch[Mvecinos[i][j]] == cityPatch[i]){
                prov2 += Mpesos[i][j];              //prov2 -> toda la movilidad del nodo dentro de su ciudad, incluidos los autoloops (?)
                if(i != Mvecinos[i][j])
                    prov1 += Mpesos[i][j];          //prov1 -> solo los que salen de su nodo
            }
        }
        if(prov2 != 0)
            mI[i] = prov1/prov2;
        else
            mI[i] = 0;
        //printf("%lf\n", mI[i]);
    }
}

void MobMatrix::setR(){
    Rout.resize(N);
    Rin.resize(N);
    for(int i = 0; i < N; i++){
        Rout[i].resize(vecinos[i]);
        Rin[i].resize(vecinos[i]);
    }

    double prov1;
    double prov2;
    for(int i = 0; i < N; i++){
        for(int j = 0; j < vecinos[i]; j++){
            //primero el Rout
            if(cityPatch[i] != cityPatch[Mvecinos[i][j]]){
                prov1 = prov2 = 0;
                prov1 = Mpesos[i][j];
                for(int k = 0; k < vecinos[i]; k++){
                    if(cityPatch[i] != cityPatch[Mvecinos[i][k]])
                        prov2 += Mpesos[i][k];
                }
                if(prov2 != 0)
                    Rout[i][j] = prov1/prov2;
                else
                    Rout[i][j] = 0;
                
                Rin[i][j] = 0; //El j pertenece a otra ciudad, por lo que no hay Rin
            }
            
            //despues el Rin
            if(cityPatch[i] == cityPatch[Mvecinos[i][j]]){
                prov1 = prov2 = 0;
                if(i == Mvecinos[i][j]){
                    Rin[i][j] = 0;
                }
                else{
                    prov1 = Mpesos[i][j];
                    for(int k = 0; k < vecinos[i]; k++){
                        if(cityPatch[i] == cityPatch[Mvecinos[i][k]] && i != Mvecinos[i][k])
                            prov2 += Mpesos[i][k];
                    }
                    if(prov2 != 0)
                        Rin[i][j] = prov1/prov2;
                    else
                        Rin[i][j] = 0;
                }
                Rout[i][j] = 0; //El j pertenece a la misma ciudad, por lo que no hay Rout
            }
            //printf("%lf\t%lf\n", Rin[i][j], Rout[i][j]);
        }
    }

}

void MobMatrix::calculaRTraspuesta(){
    RoutT.resize(N);
    RinT.resize(N);
    for(int i = 0; i < N; i++){
        RoutT[i].resize(vecinosT[i]);
        RinT[i].resize(vecinosT[i]);
    }
	int k = 0;
    int i = 0, j = 0;
    int B;
	int* AUX = (int*)malloc(N * sizeof(int));
	for (i = 0; i < N; i++)
		AUX[i] = 0;

	for (i = 0; i < N; i++){
		for (j = 0; j < vecinos[i]; j++){
            B = Mvecinos[i][j];
            MvecinosT[B][AUX[B]] = i;
            RoutT[B][AUX[B]] = Rout[i][j];
		    RinT[B][AUX[B]] = Rin[i][j];
		    AUX[B] = AUX[B] + 1;
		}
	}
	free(AUX);
}