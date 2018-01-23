#ifndef Cell_hpp
#define Cell_hpp

#include"Particle.hpp"
#include"Tenseur.hpp"

#include<iostream>
#include<fstream>
#include<cstdlib>

class Sample;


//Periodic Cell or "space" through which force/displacement are applied
//to the sample

class Cell{

	private:
		//Dimensions:
		double Lx_;
		double Ly_;
		//Coordonnees centre:
		double xc_;
		double yc_;
		//Cell mass:
		double mh_;

		//On impose soit une vitesse hd, soit une contrainte stress_ext
		//Il faut que la partie antisymetrique de Ldot soit nulle, equiv a h(hdot) symetrique impose
		char Control_[4];
		double Control_values_Init[4];

		//Metrics: collective degrees of freedom
		Tensor2x2 h_;
		//Garde en memoire la forme originelle cellule pour calcul d'engeniring strain (L(t)-L0/LO))
		Tensor2x2 h0_;
		Tensor2x2 hd_;
		Tensor2x2 hdd_;
		//Strain tensor: cumulative (time measure of the essai)
		Tensor2x2 s_;
		Tensor2x2 stress_ext;

		//Defined by initial sample
		bool L_auto_;
		bool mh_auto_;

		//Check:
		bool initCG_;
		bool initGeometry_;
		bool initMass_;

		//Folder where outpus are written:
		std::string folder_;
		std::string fcell_;
	public:
		Cell();
		//Init User control
		void init(std::ifstream&);
		void talkinit(Sample&);
		//INPUT/OUTPUTs:
		void load(const int);
		void write(const int) const;
		~Cell(){};

		void computeExternalStress(const Tensor2x2&);
		//Maj de h,hdd,Ld,s apres a la fin du pas de temps
		void updatehdd(const Tensor2x2 stress_int);
		void updatehd(const double dt);
		void firstStepVerlet(const double dt);

		//Met a jours la periodicite des particules en position
		void PeriodicBoundaries2(std::vector<Particle>*);
		void CalculStrainTensor();

		//Acces
		char getControlxx() const {return Control_[0];}
		char getControlxy() const {return Control_[1];}
		char getControlyx() const {return Control_[2];}
		char getControlyy() const {return Control_[3];}

		double getVolume() const ;
		double getMasse() const { return mh_;}

		double get_width() const;
		double get_height() const;

		double getxc() const ;
		double getyc() const ;

		Tensor2x2 geth0() const { return h0_;}
		Tensor2x2 geth() const { return h_;}
		Tensor2x2 gethd() const { return hd_;}
		Tensor2x2 gethdd() const { return hdd_;}
		const Tensor2x2& getStrainTensor() const { return s_;}

		Tensor2x2 getStressExt() const { return stress_ext;}

		void initfolder(std::string folder) { folder_ = folder;}
		void writeGeometry(const int) const;

		//Debug & track:
		bool initcheck();
		void debug(const int)const;

};

#endif
