#ifndef ListeInteraction_hpp
#define ListeInteraction_hpp

#include<vector>
#include<iostream>
#include<string>
#include<utility>
#include<algorithm>
#include<iterator>

#include"Contact.hpp"
#include"Particle.hpp"
//Manage verlet list, contact list


//Car c'est different de l'algorithme (schema d'integration de Algo)



//liste interactions potentielles


//Toute interaction est considérée comme un contact
//Parmi les contacts, les contacts reels sont considérés comme actifs


class Sample;
class Cell;

class Interactions{

	private:

	  //Juset Verlet pour l'instant
		double dv_;

		unsigned int nv_;

		//User defined: Rmin, Rmax
		std::string scale_;

		//All interactions possible
		std::vector<Contact> pairs_; 

		//Verlet list
		//Points on pairs_
		std::vector<Contact*> vlist_;

		//Active contact list: id de vlist
		//Points on pairs_
		std::vector<Contact*> clist_;

		//plug:
		Sample * spl_;
		Cell * cell_;

		//Global check at initialisation:
		bool checkInteractions_;
		//Locals check at initialisation:
		bool initScale_;
		bool initdv_;

		bool initkn_;
		bool initkt_;
		bool initgn_;
		bool initgt_;
		bool initmus_;

		bool setgnmax_;

		std::string folder_;
		std::string fInteractions_;

		//Parametres DEM:

		//Mettre en place check sur les
		//bornes acceptables des valeurs des parametres
		//a comparer avec dt
		//Etudier un peu ca...

		//Normal and tangential stifness
		double kn_;
		double kt_;
		//Normal and tangential viscosities
		double gn_;
		double gt_;
		//Interparticle friction
		//Later will depend on particle Id TODO
		double mus_;


		//Interal stress tensor:
		//static part
		Tensor2x2 stress_s;
		//Kinetic part
		Tensor2x2 stress_c;
		//TMP: total stress
		Tensor2x2 stress_;

	public:
		Interactions();
		~Interactions();
		//Initialisations:
		void init(std::ifstream&);
		void initScale();
		bool initcheck();
		bool checkDEMparameters()const;
		void plug(Sample&,Cell&);
		void build();
		//Load file
		void talkinit();
		void load(const int);
		void updateverlet(const int);
		//Build verlet
		void updatevlist();
		//Build contact list
		void detectContacts();
		//Compute forces at contacts (end of 1st step verlet algo)
		void computeForces(const double);
		void computeInternalStress();

		//Writing outputs:
		void initfolder(std::string folder) { folder_ = folder;}
		void writeContacts(int) const;
		//Compute energy stored in contact deflection
		double getElasticEnergy() const;
		std::vector<double> getAverageMaxPenetration()const;
		void askNumberOfContacts() const;

		int getnc() const { return clist_.size();}
		int getnv() const { return nv_ ;}
		double getkn() const { return kn_;}
		double getgn() const { return gn_;}
		void setgn(double gn) { gn_ = gn;}
		bool setgnmax() const { return setgnmax_;}
		const Tensor2x2& stress() const { return stress_;}
		Tensor2x2 getStressInt() const { return stress_;}

		bool near(const Particle*,const Particle*,const double) const;
		//DEBUG
		void debug(const int)const;

		Vecteur getShortestBranch(const Particle*,const Particle*) const;

		const Contact& inspectContact(int k) const { return pairs_[k] ; } 

		void read(std::ifstream& is);
};



#endif

