#ifndef Cell_hpp
#define Cell_hpp

#include"Particle.hpp"
#include"Tenseur.hpp"

#include<iostream>
#include<fstream>
#include<cstdlib>

//Cellule periodique
class Sample;

class Cell{

	private:
		//Dimensions:
		double Lx_;
		double Ly_;
		//Coordonnees centre:
		double xc_;
		double yc_;
		//Cell mass:(a ajuster plus tard)
		double mh_;

		//Imposed boundary conditions: garde en memoire ce qui est controle et libre
		//On impose soit une vitesse hdd, soit une contrainte stress_ext
		//Il faut que la partie antisymetrique de Ldot soit nulle, equiv a h(hdot) symetrique impose
		char Control_[4];

		//Metrics: collective degrees of freedom
		Tensor2x2 h_;
		//Garde en memoire la forme originelle cellule pour calcul d'engeniring strain (L(t)-L0/LO))
		Tensor2x2 h0_;
		Tensor2x2 hd_;
		Tensor2x2 hdd_;
		//Strain tensor: cumulative (time measure of the essai)
		Tensor2x2 s_;
		//Velocity gradient tensor (a voir plus tard, mais BC plus generale possible)
		//Tensor2x2 Ld_;
		//Stress: int et ext
		Tensor2x2 stress_ext;
		//Maj par particules
		//Tensor2x2 stress_int;

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
		//Init Starting Geometry and mass from sample
		void initFromSample(Sample&);
		bool needSample();
		//A l'avenir constructeur prend en argument un truc qui init/garde en memoire  les DOF controles
		//Cell(double L,Config& config): L_(L),xc_(L_/2.),yc_(L_/2.){initCell(config);}
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
		double getLx() const { return Lx_;}
		double getLy() const { return Ly_;}

		double getxc() const ;
		double getyc() const ;

		Tensor2x2 geth() const { return h_;}
		Tensor2x2 gethd() const { return hd_;}
		Tensor2x2 gethdd() const { return hdd_;}

		Tensor2x2 getStressExt() const { return stress_ext;}

		//Writins outputs:
		void initfolder(std::string folder) { folder_ = folder;}
		void write(const int) const;
		void writeStrainTensor(std::ofstream&,double);

		//Debug & track:
		void affiche(){std::cout<<xc_<<" "<<yc_<<" "<<getVolume()<<std::endl;}
		bool initcheck();
		void debug(const int)const;

};

#endif
