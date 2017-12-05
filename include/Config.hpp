#ifndef Config_hpp
#define Config_hpp

#include<fstream>
#include<iostream>

#include"Algo.hpp"
#include"Cell.hpp"
#include"Sample.hpp"
#include"Tenseur.hpp"

class Config{
	private:
		//Tableau de 4 parametres de controle (chaque direction) "v" -> vitesse "f"-> force
		//0 xx 1 xy 2 yx 3 yy
		 char BCU[4];
		 //Controle en vitesse:
		 Tensor2x2 LdUser_;
		 //Controle en force:
		 Tensor2x2 StressUser_;
		 double rho_;
	public:
		 Config();
		 ~Config();
		 Tensor2x2 returnLd() const { return LdUser_;}
		 Tensor2x2 returnStress() const { return StressUser_;}

		 char getBCxx() const { return BCU[0];}
		 char getBCxy() const { return BCU[1];}
		 char getBCyx() const { return BCU[2];}
		 char getBCyy() const { return BCU[3];}
		 char getBCU(int i) const { return BCU[i];}

		 //Interface d'initialistion
		 void init(std::ifstream&,Algo&,Cell&,Sample&);
};

#endif
