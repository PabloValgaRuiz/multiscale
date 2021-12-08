//#define topecontactos


#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <future>
#include <unistd.h>
#include "MobMatrix.hpp"
#include "MarkovDistMult.hpp"
#include "MCDistMult.hpp"
#include "iteracion.hpp"
#include "mc.hpp"
#include <exception>
#include <numeric>
#include <ios>

static std::mutex s_Mutex;
static std::mutex threads_Mutex;
static int threads;
void markov(const MobMatrix& T, double beta, double p1, double p2, int tiempo);
void iteracionThreshold(const MobMatrix& T, double p1, double p2, std::vector<std::vector<std::pair<double, double>>>& heatMap, int i, int j);
int threshold_heatmap(const MobMatrix& T);
void iteracionHeatMap(const MobMatrix& T, double beta, double p1, double p2, int tiempo, std::vector<std::vector<double>>& heatMap, int i, int j);
int heatmap(const MobMatrix& T);
void iteracionMontecarlo(const MobMatrix& T, double beta, double p1, double p2);

static const std::string state = "ny";
int main(int argc, char* argv[]){

    std::string citPat = "citiesMult/"+ state +"/final/Citypatch.txt";
    std::string mobNet = "citiesMult/"+ state +"/final/mobnetwork.txt";
    std::string popAr = "citiesMult/"+ state +"/final/Poparea.txt";

    const MobMatrix T{citPat, mobNet, popAr};

    heatmap(T);

    // double pI, pC;
    // std::cin >> pI >> pC;
    // std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    // MarkovDistMult markov{pI,pC,T};
    // auto par = iteracion(T, markov);
    // std::cout << "done" << std::endl;
    // std::ofstream file{"out/" + state + "/eigenvectorDiaNoche.txt"};
    // std::ofstream file2{"out/tx/vectorentero.txt"};
    // double temp = 0;
    // for(int i = 0; i < T.N; i++){
    //     temp = 0;
    //     for(int j = 0; j < T.N; j++){
    //         file2 << i << "\t" << j << "\t" << par.first[i*T.N + j] << std::endl;
    //     }
    //     for(int j = 0; j < T.vecinos[i]; j++){
    //         temp += par.first(i*T.N + T.Mvecinos[i][j]) * T.Mpesos[i][j];
    //     }
    //     if(T.population[i] != 0)
    //         temp /= T.population[i];
    //     file << i << "\t" << temp << std::endl;
    // }
    // file.close();
    // file2.close();
    // std::cout << par.second << std::endl;
    return 0;
}

int threshold_heatmap(const MobMatrix& T){
    const double PI_MIN = 0;
    const double PI_MAX = 1;
    const double PC_MIN = 0;
    const double PC_MAX = 1;

    const int MAX_THREADS = 24;

    std::vector<std::vector<std::pair<double, double>>> heatMap;
    heatMap.resize(65 - 1); //Numero de pI's //2^7 + 1 =  129
    for(auto& i : heatMap){i.resize(65 - 1);} //Numero de pC's

    std::ofstream fileMk("out/"+ state +"/thresholds.txt");
    std::vector<std::future<void>> futures;
    double i = 0;
    threads = 0;
    for(double pI = 1.0/64; pI <= PI_MAX; pI += (PI_MAX - 0)/64.0){
        double j = 0;
        for(double pC = 1.0/64; pC <= PC_MAX; pC += (PC_MAX - 0)/64.0){
            std::cout << i  << "    " << j << std::endl;
            while(threads >= MAX_THREADS){usleep(1000);} //Esperar a que baje el numero de hilos (LINUX)

            threads_Mutex.lock();
            threads++;  //Sumar un hilo
            futures.push_back(std::async(std::launch::async, iteracionThreshold, std::ref(T), pI, pC, std::ref(heatMap), i, j));
            threads_Mutex.unlock();
            //La funcion al terminar resta el hilo
            j++;
        }
        futures.clear();
        j = 0;
        for(double pC = 1.0/64; pC <= PC_MAX; pC += (PC_MAX - PC_MIN)/64.0){
            fileMk << pI << "\t" << pC << "\t" << heatMap[i][j].first << "\t" << heatMap[i][j].second << "\n";
            j++;
        }
        fileMk << std::endl;
        i++;
    }
    return 0;
}

void iteracionThreshold(const MobMatrix& T, double p1, double p2, std::vector<std::vector<std::pair<double, double>>>& heatMap, int i, int j){
{
    MarkovDistMult markov{p1,p2,T};
    //Cada hilo cambia un par (i,j) distinto, no es necesario un mutex
    beginning:
    try{
    auto par = iteracion(T, markov);
    heatMap[i][j].first = 0;
    std::vector<double> vtemp; vtemp.resize(T.N);
    double temp = 0;
    for(int I = 0; I < T.N; I++){
        vtemp[I] = 0;
        for(int J = 0; J < T.vecinos[I]; J++){
            vtemp[I] += par.first[I*T.N + T.Mvecinos[I][J]] * T.Mpesos[I][J];
        }
        if(T.population[I] != 0){
            vtemp[I] = vtemp[I] / T.population[I];
        }
        temp += pow(vtemp[I], 2);
    }
    double norm = sqrt(temp);
    temp = 0;
    for(auto& v : vtemp){
        if(norm != 0)
            v /= norm;
        temp += pow(v, 4);
    }
    heatMap[i][j].first = sqrt(temp);
    
    heatMap[i][j].second = par.second;
    }
    catch(std::bad_alloc& except){
        std::cout << "Bad alloc exception thrown on " << i << ", " << j << ", reason: " << except.what() << std::endl;
        goto beginning;
    }
    catch(std::exception& except){
        std::cout << "Exception thrown on " << i << ", " << j << ", reason: " << except.what() << std::endl;
        goto beginning;
    }
    catch(...){
        std::cout << "Unknown exception thrown on " << i << ", " << j << std::endl;
        goto beginning;
    }
}   //Cerrar el scope antes de disminuir el valor de threads, para que se libere la memoria de las variables locales antes de crear el siguiente
    std::lock_guard<std::mutex> lock(threads_Mutex);
    threads--;
}

int heatmap(const MobMatrix& T){
    const double BETA_MIN = 0;
    const double BETA_MAX = 5;
    const double P_MIN = 0;
    const double P_MAX = 1;

    int tiempo = 400;
    const int MAX_THREADS = 24;

    std::vector<std::vector<double>> heatMap;
    heatMap.resize(129); //Numero de p's //2^7 + 1 = 129
    for(auto& i : heatMap){i.resize(129);} //Numero de lambdas
    
    std::ofstream fileMk("out/"+ state +"/heatmap_pIpC.txt");
    std::vector<std::future<void>> futures;
    double i = 0;
    threads = 0;

    for(double p = P_MIN; p <= P_MAX; p += (P_MAX - P_MIN)/(heatMap.size() - 1)){
        double pI = p;
        double pC = p;
        double j = 0;
        for(double beta = BETA_MIN; beta <= BETA_MAX; beta += (BETA_MAX - BETA_MIN)/(heatMap[i].size() - 1)){
            std::cout << i  << "    " << j << std::endl;
            while(threads >= MAX_THREADS){usleep(1000);} //Esperar a que baje el numero de hilos (LINUX)

            threads_Mutex.lock();
            threads++;  //Sumar un hilo
            futures.push_back(std::async(std::launch::async, iteracionHeatMap, std::ref(T), beta, pI, pC, tiempo, std::ref(heatMap), i, j));
            threads_Mutex.unlock();
            //La funcion al terminar resta el hilo
            j++;
        }
        futures.clear();

        j = 0;
        for(double beta = BETA_MIN; beta <= BETA_MAX; beta += (BETA_MAX - BETA_MIN)/(heatMap[i].size() - 1)){
            fileMk << p << "\t" << beta << "\t" << heatMap[i][j] << "\n";
            j++;
        }
        fileMk << "\n";
        fileMk.flush();
        i++;
    }
    return 0;

}

void iteracionHeatMap(const MobMatrix& T, double beta, double p1, double p2, int tiempo, std::vector<std::vector<double>>& heatMap, int i, int j){
{
    MarkovDistMult markov{p1, p2, T};
    markov.inicializar(0.002);
    markov.setBeta(beta);
    beginning:
    try{
        for(int t = 0; t < tiempo; t++){
            markov.iteracion(T);
        }
        markov.contarInfectados(T);
    }catch(std::exception except){
        std::cout << "Exception thrown on " << i << ", " << j << ", reason: " << except.what() << std::endl;
        goto beginning;
    }
    //Cada hilo cambia un par (i,j) distinto, no es necesario un mutex
    heatMap[i][j] = markov.infectadosTotal;
}//Cerrar el scope antes de disminuir el valor de threads, para que se libere la memoria de las variables locales antes de crear el siguiente
    std::lock_guard<std::mutex> lock(threads_Mutex);
    threads--;
}

void markov(const MobMatrix& T, double beta, double p1, double p2, int tiempo)
{
    MarkovDistMult markov{p1, p2, T};
    markov.inicializar(0.0000001);
    markov.setBeta(beta);

    try{
        for(int t = 0; t < tiempo; t++){
            markov.iteracion(T);
        }
        markov.contarInfectados(T);
    }catch(std::exception except){
        std::cout << "Exception thrown, reason: " << except.what() << std::endl;
        exit(1);
    }catch(...){
        std::cout << "Exception thrown" << std::endl;
        exit(1);
    }
    
    std::ofstream file("out/"+ state +"/markov.txt");
    for(int i = 0; i < T.N; i++){
        file << i << "\t" << markov.infectados[i] << std::endl;
    }
    file.close();
}

void iteracionMontecarlo(const MobMatrix& T, double beta, double p1, double p2){

    MCDistMult montecarlo(p1, p2, T);
    montecarlo.setBeta(beta);
    montecarlo.inicializar(0.0002);
    for(int t = 0; t < 300; t++){
        montecarlo.iteracion(T);
    }
    for(int t = 0; t < 500; t++){

    }



}