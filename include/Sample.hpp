#ifndef SAMPLE_HPP
#define SAMPLE_HPP

#include<vector>
#include<fstream>
#include<iostream>
#include"Particle.hpp"
#include"globalfunctions.hpp"

class Cell;

class Sample{

	private:
		std::vector<Particle> spl_;
		//Input sample initial:
		std::string fsampleIni_;

		//outputs:
		std::string folder_;
		//Reduced coordinates
		std::string fsample_;
		//Absolute coordinates
		std::string fsampleA_;

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

		//If load or build from init sample
		bool load_sample_;
		//File to load if load
		unsigned int filetoload_;
		//Int to start writing if load
		unsigned int starting_;

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
		void initReducedCoordinates();
		void plugtoCell(Cell&);
		void attributeMass();
		void setminmax();
		bool initcheck();
		bool isEmptySampleFile(std::ifstream&);

		double getxmin() const { return xmin_;}
		double getxmax() const { return xmax_;}
		double getymin() const { return ymin_;}
		double getymax() const { return ymax_;}
		double getrmax() const { return rmax_;}
		double getrmin() const { return rmin_;}
		double getMass() const { return M_;}
		double getrho() const { return rho_;}

		unsigned int getsize() const { return spl_.size();}
		bool loaded() const { return load_sample_;}
		unsigned int filetoload() const { return filetoload_;}
		unsigned int startingTic() const { return starting_;}

		const std::vector<Particle>& inspectSample() const {return spl_;}
		//Only used to apply periodicity (used by cell)
		//Should be moved to sample
		std::vector<Particle>* getSample() {return &spl_;}

		void firstStepVerlet(const double);
		void secondStepVerlet(const double);

		Cell * getCell() const { return cell_;}

		//Writing outputs:
		void initfolder(std::string folder); 
		void write(int) const;
		void writeAbsolute(int) const;

		//TMP: debug writing
		void debug(int);
		void printSample() const;

		std::vector<Particle> getimages(double e) const;
		double getTKE() const;
		double getRKE() const;
		Vecteur returnrabs(const Particle&) const;
		Vecteur returnvabs(const Particle&) const;

		//Post processing:
		//Used to load network:
		//Not used
		Particle * getP(int id);
		void setrho(double rho) { rhodefined_ = true ; rho_ = rho;}
		void setfiletoload(const int i) { load_sample_ = true ; filetoload_ = i ;}


};


#endif
