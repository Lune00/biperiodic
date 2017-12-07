#ifndef ListeInteraction_hpp
#define ListeInteraction_hpp

#include<vector>
#include<iostream>
#include<string>

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

		//plug:
		Sample * spl_;


		//Global check at initialisation:
		bool checkInteractions_;
		//Locals check at initialisation:
		bool initScale_;
		bool initdv_;
		bool initdsv_;


	public:
		Interactions();
		~Interactions();
		void updateverlet(const int);
		void updatevlist();
		void updatesvlist();
		int getnv() const { return nv_ ;}
		int getnsv() const {return nsv_;}
		void init(std::ifstream&);
		void initScale();
		bool initcheck();
		void plug(Sample&);
		bool near(const Particle&,const Particle&,const Tensor2x2&,const double) const;



};



#endif
