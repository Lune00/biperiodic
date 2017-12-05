#ifndef SAMPLE_HPP
#define SAMPLE_HPP

#include"Particle.hpp"
#include<vector>
#include<fstream>
#include<iostream>

class Sample{

	private:
		std::vector<Particle> spl_;
		std::string fichsample_;

		double xmin_;
		double xmax_;
		double ymin_;
		double ymax_;
		double rmin_;
		double rmax_;

		//Masse volumique particules:
		double rho_;

		//Check:
		bool sampleIsLoaded_;
		bool sampleIsFilled_;
		//Assure utilisateur defini rho
		bool rhodefined_;

	public:
		Sample();
		~Sample();
		void init(std::ifstream&);
		//Return 0 si ok, 1 sinon
		void loadSample();
		void write(std::ofstream&);
		void attributeMass();
		void setminmax();
		bool initcheck();
		bool isEmptySampleFile(std::ifstream&);

};


#endif
