#ifndef Cell_hpp
#define Cell_hpp

#include"Config.hpp"
#include"Particle.hpp"

#include<iostream>
#include<fstream>

//Cellule periodique
class Cell{

	private:
		//Dimension:
		double L_;
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
		//Strain tensor: cumulative (time measure of the essai)
		Tensor2x2 s_;
		//Velocity gradient tensor (a voir plus tard, mais BC plus generale possible)
		Tensor2x2 Ld_;
		//Stress: int et ext
		Tensor2x2 stress_ext;
		//Maj par particules
		Tensor2x2 stress_int;
	public:
		Cell(Config&);
		//A l'avenir constructeur prend en argument un truc qui init/garde en memoire  les DOF controles
		Cell(double L,Config& config): L_(L),xc_(L/2.),yc_(L/2.){initCell(config);}
		~Cell(){};

		//Init les parametres BC et de la grille par utilisateur
		void initCell(Config&);
		//Controle force/vitesse
		void ApplyBC();
		//Maj de h,hdd,Ld,s apres a la fin du pas de temps
		void update(Tensor2x2,Tensor2x2,double);

		//Met a jours la periodicite des particules
		void PeriodicBoundaries(Particle&);

		void write(std::ofstream&,std::ofstream&,double);

		//Acces
		double getVolume();
		double getMasse() const { return mh_;}
		double getL() const { return L_;}
		double getL2() const { return L_/2.;}
		Tensor2x2 geth() const { return h_;}
		Tensor2x2 gethd() const { return hd_;}
		Tensor2x2 getStressInt() const { return stress_int;}
		Tensor2x2 getStressExt() const { return stress_ext;}

		//Debug & track:
		void affiche(){std::cout<<xc_<<" "<<yc_<<" "<<getVolume()<<std::endl;}

};

#endif
