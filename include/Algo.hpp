#ifndef Algo_hpp
#define Algo_hpp

#include<vector>
#include<iostream>
#include<fstream>
#include<iomanip>

class Cell;
class Sample;
class Interactions;
class Analyse;

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
		//Frequency analyse during simulation
		int nana_;

		//Computed from Interactions/dt
		//Global parameters to print:
		//All writiin during computed during initcheck
		//calling checkSimulationParameters

		//Normal restitution coefficient
		double e_;
		//Max normal viscosity
		double gnmax_;
		//Max time step
		double dtmax_;


		//File where initial setup is written:
		std::string fsetup_;

		//Plug:
		Interactions * Int_;
		Cell * cell_;
		Sample * spl_;
		Analyse * ana_;

	public:
		Algo() { dt_ = 1. ; ns_ = 0 ; nrecord_ = 0; tic_= 0; t_=0.; ticw_ = 0; nana_ = 0 ; fsetup_ = "simusetup.txt";}
		~Algo(){};
		//void verletalgo(Cell&,std::vector<Particle>&);
		void verletalgo2();
		void init(std::ifstream&);
		bool initcheck();
		void run();
		void plug(Cell&,Sample&,Interactions&,Analyse&);
		void write();
		bool checkSimulationParameters();
		void computedtmax();
		void compute_gnmax_restitution();
		bool checktimestep()const;
		bool checkNormalViscosity() const;
		//Write simu setup in filesetup_;
		void writesetup() const;
};

#endif
