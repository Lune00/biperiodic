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

	public:
		Sample();
		~Sample();
		void load(std::istream&);
		void write(std::ofstream&);
		void attributeMass(double density);

	private:
		std::vector<Particle> spl_;
};


#endif
