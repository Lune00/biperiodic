#ifndef Algo_hpp
#define Algo_hpp

#include<vector>
#include<iostream>
#include<fstream>

class Cell;
class Sample;
class Interactions;

class Algo{
	private:
		double dt_;
		double t_;
		int tic_;
		//Numerotation files
		int ticw_;
		int ns_;
		//Frequency recording sample+network
		int nrecord_;

		//Parametres DEM
		double Param_[10];

		//Plug:
		Interactions * Int_;
		Cell * cell_;
		Sample * spl_;

	public:
		Algo() { dt_ = 1. ; ns_ = 0 ; nrecord_ = 0; tic_= 0; t_=0.; ticw_ = 0;}
		~Algo(){};
		//void verletalgo(Cell&,std::vector<Particle>&);
		void verletalgo2();
		void init(std::ifstream&);
		bool initcheck();
		void run();
		void plug(Cell&,Sample&,Interactions&);
		void write();
};

#endif
