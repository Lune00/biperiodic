#ifndef ListeInteraction_hpp
#define ListeInteraction_hpp

#include<vector>
#include<iostream>

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

class Interactions{

	private:

		struct particle_pair{
			Particle * i;
			Particle * j;
		};

		double dv_;
		double dsv_;

		unsigned int nv_;
		unsigned int nsv_;

		double Rmax_;

		//Verlet list
		std::vector<Contact> vlist_;
		//Super Verlet list
		std::vector<Contact> svlist_;
		//Active contact list: id de vlist
		std::vector<int> clist_;

		//plug:
		Sample * spl_;


	public:
		Interactions();
		~Interactions();
		void updatevlist();
		int getnv() const { return nv_ ;}
		int getnsv() const {return nsv_;}
		void plug(Sample&);



};



#endif
