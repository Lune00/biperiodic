#ifndef SAMPLE_HPP
#define SAMPLE_HPP

#include"Particle.hpp"
#include<vector>
#include<fstream>
#include<iostream>

//Lit un fichier préalablement genere
//Coordonnees reduites
//id r sx sy sdx sdy 

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

	public:
		Sample();
		~Sample();
		void init(std::ifstream&);
		void loadSample(std::ifstream&);
		void write(std::ofstream&);
		void attributeMass();
		void setminmax();

};


#endif
