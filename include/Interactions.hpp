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

		double dv_;
		double dsv_;

		unsigned int nv_;
		unsigned int nsv_;

		struct particle_pair{
		  Particle * i ;
		  Particle * j ;
		};

		//User defined: Rmin, Rmax
		std::string scale_;

		//SVerlet list
		std::vector<particle_pair> svlist_;
		//Verlet list
		std::vector<Contact> vlist_;

		//Active contact list: id de vlist
		//Points on vlist_
		std::vector<int> clist_;

		//plug:
		Sample * spl_;
		Cell * cell_;

		//Store tangential deflection during integration
		double * dts_ ;
		//Number of particles (size of dts_)
		unsigned int N_ ;

		//Global check at initialisation:
		bool checkInteractions_;
		//Locals check at initialisation:
		bool initScale_;

		//Verlet et Sverlet init
		bool initdv_;
		bool initdsv_;

		bool initkn_;
		bool initkt_;
		bool initgn_;
		bool initgt_;
		bool initmus_;

		bool setgnmax_;
		bool setgtmax_;

		//Output:
		std::string folder_;
		std::string fInteractions_;

		//Parametres DEM:

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
		void build();
		void load(const int);
		void read_dt(std::ifstream&);

		//Post-processing:
		void loadnetwork(const int);
		void read_contact(std::ifstream&);

		//Verlet and SuperVerlet list managment:
		void updatevlist();
		void updatesvlist();
		void updateverlet(const int);

		//Build contact list
		void detectContacts();

		//Compute forces at contacts (end of 1st step verlet algo)
		void computeForces(const double,const double);
		void computeInternalStress();
		void addForce(Particle&,const double);

		//Writing outputs:
		void initfolder(std::string folder) { folder_ = folder;}
		void writeContacts(int) const;

		//Compute energy stored in contact deflection
		double getElasticEnergy() const;

		//Accesors/mutators:
		std::vector<double> getAverageMaxPenetration()const;
		void askNumberOfContacts() const;

		int getnc() const { return clist_.size();}
		int getnv() const { return nv_ ;}

		double getkn() const { return kn_;}
		double getkt() const { return kt_;}
		
		double getgn() const { return gn_;}
		double getgt() const { return gt_;}
		double getmu() const { return mus_;}

		void setgn(double gn) { gn_ = gn;}
		void setgt(double gt) { gt_ = gt;}

		bool setgnmax() const { return setgnmax_;}
		bool setgtmax() const { return setgtmax_;}
		double sign(double x){if(x<0.) return -1.;else return 1.;}

		const Tensor2x2& stress() const { return stress_;}
		Tensor2x2 getStressInt() const { return stress_;}
		Tensor2x2 getStressK() const { return stress_c;}
		Tensor2x2 getStressS() const { return stress_s;}

		const Contact* inspectContact(int k) const { return &vlist_[clist_[k]];}

		void read(std::ifstream& is);
		void print()const;

		//Store/manage tangential contact deflection:
		void init_array_dt();
		double get_dt(Contact&) const;
		void set_dt(Contact&);

		//Not used
		void reset_dt(Contact&);

		//DEBUG
		void debug(const int)const;

		//Post-processing
		void setparameters(double kn,double kt,double gn,double gt,double mu);

};



#endif

