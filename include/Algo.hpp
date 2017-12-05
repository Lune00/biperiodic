#ifndef Algo_hpp
#define Algo_hpp

#include<vector>
#include<iostream>
#include<fstream>

class Cell;
class Particle;

class Algo{
	private:
		double dt_;
		int ns_;
		//Detection contacts
		int actualiseVerlet_;
		double distVerlet_;

		//Parametres DEM
		double Param_[10];

	public:
		Algo() { dt_ = 1. ; ns_ = 0 ;}
		~Algo(){};
		//void verletalgo(Cell&,std::vector<Particle>&);
		void verletalgo2(Cell&,std::vector<Particle>&);
		void init(std::ifstream&);
		bool initcheck();
		void run();
};

#endif
