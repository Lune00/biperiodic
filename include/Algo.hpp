#ifndef Algo_hpp
#define Algo_hpp
#include<vector>
#include<iostream>

class Cell;
class Particle;

class Algo{
	private:
		double dt_;
		//Detection contacts
		int actualiseVerlet_;
		double distVerlet_;

		//Parametres DEM
		double Param_[10];

	public:
		Algo() { dt_ = 1. ;}
		Algo(double dt): dt_(dt) {};
		~Algo(){};
		void verletalgo(Cell&,std::vector<Particle>&);
		void verletalgo2(Cell&,std::vector<Particle>&);
		void run();

};

#endif
