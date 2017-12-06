#ifndef SAMPLE_HPP
#define SAMPLE_HPP

#include"Particle.hpp"
#include<vector>
#include<fstream>
#include<iostream>


class Cell;

class Sample{

	private:
		std::vector<Particle> spl_;
		std::string fichsample_;
		std::string folder_;

		double xmin_;
		double xmax_;
		double ymin_;
		double ymax_;
		double rmin_;
		double rmax_;

		//Masse volumique particules:
		double rho_;
		//Sample total mass
		double M_;

		//Check:
		bool sampleIsLoaded_;
		bool sampleIsFilled_;
		//Assure utilisateur defini rho
		bool rhodefined_;

		Cell* cell_;

	public:
		Sample();
		~Sample();

		void init(std::ifstream&);
		void loadSample();
		void initReducedCoordinates(Cell&);
		void plugtoCell(Cell&);
		void write(std::ofstream&);
		void writeAbsolute(std::ofstream&);
		void attributeMass();
		void setminmax();
		void initfolder(std::string folder) { folder_ = folder;}

		bool initcheck();
		bool isEmptySampleFile(std::ifstream&);


		double getxmin() const { return xmin_;}
		double getxmax() const { return xmax_;}
		double getymin() const { return ymin_;}
		double getymax() const { return ymax_;}
		double getMass() const { return M_;}

		unsigned int getsize() const { return spl_.size();}
		//Reflechir a ne pas briser l'encapsulation plus tard...
		//Ca me plait moyen de donner mon vecteur de particules...
		std::vector<Particle>* getSample() { return &spl_;}

};


#endif
