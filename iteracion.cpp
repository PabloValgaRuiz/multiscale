/*
Unidad de traduccion del algoritmo de obtencion de autovalores
Tarda bastante en compilar
*/

#include "iteracion.hpp"
#include <precompiled.hpp>

unsigned long long countNonZeros(const MobMatrix& T){
	unsigned long long result = 0;
	for(int i = 0; i < T.N; i++){
		for(int j = 0; j < T.vecinos[i]; j++){
			result += T.vecinos[i]
				+ T.vecinos[T.Mvecinos[i][j]]
				+ T.vecinosT[i]
				+ T.vecinosT[T.Mvecinos[i][j]]; 
		}
	}
	return result;
}

std::pair<Eigen::ArrayXd, double> iteracion_p0(const MobMatrix& T, const MarkovDistMult& markov){
	double temp = 0, tempN = 0;
    int I = 0;
	double F = 0, SIGMA = 0;
    for(int i = 0; i < T.N; i++){
        temp = markov.get_fvector()[i];
        tempN = markov.get_sigma()[i];
        if(markov.get_zD()*temp + markov.get_zN()*tempN > markov.get_zD()*F + markov.get_zN()*SIGMA){F = temp; SIGMA = tempN; I = i;}
    }
    
	Eigen::ArrayXd eigenvector{};
	eigenvector.setConstant(T.N * T.N, 0);
	double componente = 1.0 / sqrt(T.vecinos[I]); //para que tenga norma = 1
	for(int j = 0; j < T.vecinos[I]; j++){
		eigenvector[I*T.N + T.Mvecinos[I][j]] = componente;
	}
	return {eigenvector, 1.0};
}
std::pair<Eigen::ArrayXd, double> iteracion(const MobMatrix& T, const MarkovDistMult& markov)
{
	if((markov.get_pC() == markov.get_pI()) && (markov.get_pI() == 0)){
		return iteracion_p0(T, markov);
	}

    //==============================MATRIX CREATION====================================================//

	std::vector<Eigen::Triplet<double> > tripletList;
	unsigned long long sizeOfTripletList = countNonZeros(T);
	std::cout << "MegaBytes para almacenar:\t" << sizeOfTripletList * sizeof(Eigen::Triplet<double>)/1000000 << "\n";
	tripletList.reserve(sizeOfTripletList);
	Eigen::SparseMatrix<double, Eigen::ColMajor> MATRIX(T.N * T.N, T.N * T.N);	//SI SE UTILIZA ROWMAJOR CRASHEA EL GENEIGSSOLVER


	int N = T.N;
	//Version compacta: eliminando las columnas cuyas filas asociadas son cero
	for(int i = 0; i < T.N; i++){
		double tempNoche = 0;
		if(T.population[i] != 0)
			tempNoche = markov.get_zN() * markov.get_sigma()[i] / T.population[i];
		else tempNoche = 0;
		
		for(int j = 0; j < T.vecinos[i]; j++){
			double temp1P = 0, tempP = 0;

			if(T.cityPatch[i] == T.cityPatch[T.Mvecinos[i][j]]){
				if(markov.get_neff()[i] != 0)
					temp1P = (1 - markov.get_pI()) * markov.get_fvector()[i] * markov.get_zD() / markov.get_neff()[i];
				if(markov.get_neff()[T.Mvecinos[i][j]] != 0)
					tempP = markov.get_pI() * markov.get_fvector()[T.Mvecinos[i][j]] * markov.get_zD() / markov.get_neff()[T.Mvecinos[i][j]];
			}
			else {
				if(markov.get_neff()[i] != 0)
					temp1P = (1 - markov.get_pC()) * markov.get_fvector()[i] * markov.get_zD() / markov.get_neff()[i];
				if(markov.get_neff()[T.Mvecinos[i][j]] != 0)
					tempP = markov.get_pC() * markov.get_fvector()[T.Mvecinos[i][j]] * markov.get_zD() / markov.get_neff()[T.Mvecinos[i][j]];
			}

			for(int k = 0; k < T.vecinos[i]; k++){
				if(T.cityPatch[i] == T.cityPatch[T.Mvecinos[i][k]]){
					tripletList.emplace_back(i * N + T.Mvecinos[i][j], i * N + T.Mvecinos[i][k],
					(temp1P * (1 - markov.get_pI()) + tempNoche) * T.Mpesos[i][k] );
				}
				else tripletList.emplace_back(i * N + T.Mvecinos[i][j], i * N + T.Mvecinos[i][k],
					(temp1P * (1 - markov.get_pC()) + tempNoche) * T.Mpesos[i][k] );
			}
			for(int k = 0; k < T.vecinos[T.Mvecinos[i][j]]; k++){
				if(T.cityPatch[T.Mvecinos[i][j]] == T.cityPatch[T.Mvecinos[T.Mvecinos[i][j]][k]]){
					tripletList.emplace_back(i * N + T.Mvecinos[i][j], T.Mvecinos[i][j] * N + T.Mvecinos[T.Mvecinos[i][j]][k],
					tempP * (1 - markov.get_pI()) * T.Mpesos[T.Mvecinos[i][j]][k] );
				}
				else tripletList.emplace_back(i * N + T.Mvecinos[i][j], T.Mvecinos[i][j] * N + T.Mvecinos[T.Mvecinos[i][j]][k],
					tempP * (1 - markov.get_pC()) * T.Mpesos[T.Mvecinos[i][j]][k] );
			}
			for(int k = 0; k < T.vecinosT[i]; k++){
				if(T.cityPatch[i] == T.cityPatch[T.MvecinosT[i][k]]){
					tripletList.emplace_back(i * N + T.Mvecinos[i][j], T.MvecinosT[i][k] * N + i,
					temp1P * markov.get_pI() * T.MpesosT[i][k] );
				}
				else{
					tripletList.emplace_back(i * N + T.Mvecinos[i][j], T.MvecinosT[i][k] * N + i,
					temp1P * markov.get_pC() * T.MpesosT[i][k] );
				}
			}
			for(int k = 0; k < T.vecinosT[T.Mvecinos[i][j]]; k++){
				if(T.cityPatch[T.Mvecinos[i][j]] == T.cityPatch[T.MvecinosT[T.Mvecinos[i][j]][k]]){
					tripletList.emplace_back(i * N + T.Mvecinos[i][j], T.MvecinosT[T.Mvecinos[i][j]][k] * N + T.Mvecinos[i][j],
					tempP * markov.get_pI() * T.MpesosT[T.Mvecinos[i][j]][k] );
				}
				else {
					tripletList.emplace_back(i * N + T.Mvecinos[i][j], T.MvecinosT[T.Mvecinos[i][j]][k] * N + T.Mvecinos[i][j],
					tempP * markov.get_pC() * T.MpesosT[T.Mvecinos[i][j]][k] );
				}
			}
		}
	}

	/*
	for(int i = 0; i < T.N; i++){
		for(int j = 0; j < T.N; j++){

			double temp1P = 1, tempP = 1;
			if(T.cityPatch[i] == T.cityPatch[j]){
				temp1P = (1 - markov.get_pI()) * markov.get_fvector()[i];
				tempP = markov.get_pI() * markov.get_fvector()[j];
			}
			else {
				temp1P = (1 - markov.get_pC()) * markov.get_fvector()[i];
				tempP = markov.get_pC() * markov.get_fvector()[j];
			}


			for(int k = 0; k < T.vecinos[i]; k++)
			{
				//if(mPesosDense(i,j) != 0)//Eliminar las componentes sin poblacion
				if(T.cityPatch[i] == T.cityPatch[T.Mvecinos[i][k]]){
					tripletList.emplace_back(i * N + j, i * N + T.Mvecinos[i][k],
					temp1P * (1 - markov.get_pI()) * T.Mpesos[i][k] );
				}
				else tripletList.emplace_back(i * N + j, i * N + T.Mvecinos[i][k],
					temp1P * (1 - markov.get_pC()) * T.Mpesos[i][k] );
			}

			for(int k = 0; k < T.vecinos[j]; k++){
				//if(mPesosDense(j,i) != 0)
				if(T.cityPatch[j] == T.cityPatch[T.Mvecinos[j][k]]){
					tripletList.emplace_back(i * N + j, j * N + T.Mvecinos[j][k],
					tempP * (1 - markov.get_pI()) * T.Mpesos[j][k] );
				}
				else tripletList.emplace_back(i * N + j, j * N + T.Mvecinos[j][k],
					tempP * (1 - markov.get_pC()) * T.Mpesos[j][k] );
			}

			for(int k = 0; k < T.vecinosT[i]; k++)
			{
				if(T.cityPatch[i] == T.cityPatch[T.MvecinosT[i][k]]){
					tripletList.emplace_back(i * N + j, T.MvecinosT[i][k] * N + i,
					temp1P * markov.get_pI() * T.MpesosT[i][k] );
				}
				else{
					tripletList.emplace_back(i * N + j, T.MvecinosT[i][k] * N + i,
					temp1P * markov.get_pC() * T.MpesosT[i][k] );
				}
			}
			for(int k = 0; k < T.vecinosT[j]; k++){
				//if(mPesosDense(j,i) != 0)
				if(T.cityPatch[j] == T.cityPatch[T.MvecinosT[j][k]]){
					tripletList.emplace_back(i * N + j, T.MvecinosT[j][k] * N + j,
					tempP * markov.get_pI() * T.MpesosT[j][k] );
				}
				else {
					tripletList.emplace_back(i * N + j, T.MvecinosT[j][k] * N + j,
					tempP * markov.get_pC() * T.MpesosT[j][k] );
				}
			}
		}
	}
	*/
	
	MATRIX.setFromTriplets(tripletList.begin(), tripletList.end());
	tripletList.clear();

	//==============================EIGEN CALCULATION===================================================//

	// Construct matrix operation object using the wrapper class SparseGenMatProd
	Spectra::SparseGenMatProd<double> op(MATRIX);
	// Construct eigen solver object, requesting the largest eigenvalues
	Spectra::GenEigsSolver< double, Spectra::LARGEST_MAGN, Spectra::SparseGenMatProd<double> > eigs(&op, 1, 10);
	// Initialize and compute
	eigs.init();
	//std::cout << "Computing..." << std::endl;
	int nconv = eigs.compute();								//SI SE UTILIZA ROWMAJOR CRASHEA EL GENEIGSSOLVER

	double umbral = markov.get_mu()/(eigs.eigenvalues().real().coeffRef(0) * markov.get_Beta0());

	return {eigs.eigenvectors().real().array().abs(), umbral};
}