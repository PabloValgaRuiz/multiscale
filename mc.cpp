#include "mc.hpp"

static std::mutex s_Mutex;
static std::string state;
static double montecarlo0, markov0;
static int termalizacion = 150; //150

static void iteracionesMC(const MobMatrix& T, double &infectados, double beta, double pI, double pC, int tiempo, int i){
    MCDistMult montecarlo{pI, pC, T};
    montecarlo.setBeta(beta);
    montecarlo.inicializar(0.002);
    montecarlo.calcInfectados(T);
    for(int t = 0; t < tiempo && montecarlo.infectadosTotal > 0; ++t){
        montecarlo.iteracion(T);
        montecarlo.calcInfectados(T);
    }
    std::lock_guard<std::mutex> lock(s_Mutex);
    infectados += montecarlo.infectadosTotal;
}
int mc(const MobMatrix& T, const std::string state){   //NO VALIDO PARA SIS

    double pI = 0.125;
    double pC = 0.125;

    int tiempo = 400;
    int hilos = 24;

    //std::ofstream fileMC("out/MC_00_p" + name + ".txt");
    std::ofstream fileMC("out_mobdef/" + state + "/montecarlo_p02.txt");
    std::cout << "pI = " << pI << std::endl;
    std::cout << "pC = " << pC << std::endl;

    for(double beta = 0; beta <= 5; beta += 5.0/16){        
        std::cout << "beta = " << beta << std::endl;
        double infectados = 0;
        std::vector<std::future<void>> futures;
        for(int i = 0; i < hilos; i++)
            futures.push_back(std::async(std::launch::async, iteracionesMC, std::ref(T), std::ref(infectados), beta, pI, pC, tiempo, i));
        futures.clear();

        fileMC << beta << " " << infectados/hilos << std::endl;
    }
    fileMC.close();    

    return 0;
}