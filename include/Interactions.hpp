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


//Peut etre algo de calcul des forces egalement
//Car c'est different de l'algorithme (schema d'integration de Algo)

//A besoin de communiquer avec Sample
//On va le plug a Sample


//liste interactions potentielles


//Toute interaction est considérée comme un contact
//Parmi les contacts, les contacts reels sont considérés comme actifs

//Liste de verlet et superListe de verlet (pour améliorer le temps)

class Sample;
class Cell;

class Interactions{

	private:

		//const Particle * i (svlist do not change particle)
		//but vlist need svlist and allow to change particle state
		struct particle_pair{
			Particle * i;
			Particle * j;
		};

		double dv_;
		double dsv_;

		unsigned int nv_;
		unsigned int nsv_;

		//factor for dv_ and dsv_ for list construction
		//User defined: Rmin, Rmax
		std::string scale_;

		//Verlet list
		std::vector<Contact> vlist_;
		//Super Verlet list
		std::vector<particle_pair> svlist_;
		//Active contact list: id de vlist
		std::vector<int> clist_;
		//std::vector<Contact> clist_;


		//Store tangential displacement for all possible contacts
		double * array_dt; 
		//Number of particles -> size of array_dt == size of spl
		int N_;

		//plug:
		Sample * spl_;
		Cell * cell_;


		//Global check at initialisation:
		bool checkInteractions_;
		//Locals check at initialisation:
		bool initScale_;
		bool initdv_;
		bool initdsv_;

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
		//Load file
		void talkinit();
		void load(const int);
		//Call updatevlist & updatsvlist
		void updateverlet(const int);
		//Build superverlet
		void updatevlist();
		//Build verlet
		void updatesvlist();
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
		int getnsv() const {return nsv_;}
		double getkn() const { return kn_;}
		double getgn() const { return gn_;}
		void setgn(double gn) { gn_ = gn;}
		bool setgnmax() const { return setgnmax_;}
		const Tensor2x2& stress() const { return stress_;}
		Tensor2x2 getStressInt() const { return stress_;}

		bool near(const Particle&,const Particle&,const Tensor2x2&,const double) const;
		//DEBUG
		void debug(const int)const;

		Vecteur getShortestBranch(const Particle&,const Particle&) const;

		void init_array_dt();
		//Used when reload file for simulation (load prev dt_)
		void fill_array_dt(std::ifstream&);
		//Return dt at contact, between particle i and j
		double get_dt(Contact&) const;
		void set_dt(Contact&) ;


		const Contact& inspectContact(int k) const { return vlist_[k] ; } 

};



#endif

