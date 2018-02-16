#ifndef hpp_ANALYSE_hpp
#define hpp_ANALYSE_hpp

#include<vector>
#include<fstream>
#include"Probe.hpp"
#include"globalfunctions.hpp"


class Sample;
class Cell;
class Particle;
class Interactions;


//Ensuite faire une app aussi qui fait l'analyse en relecture
//des fichiers issus de la simulation. Possiblité de découpler
//les deux

class Analyse{

	private:

		Sample* spl_;
		Interactions* Int_;
		Cell* cell_;
		std::string folder_;
		//Thickness (in Rmax) of periodic band around
		//sample visualisation
		double e_;

		//Analyse types:
		bool printSample_;
		bool printh_;
		bool strain_;
		bool stress_;
		bool energy_;
		bool compacity_;
		bool fabric_;


		//Velocity Profile
		bool SP_;
		int nbinsSP_;

		//Interpenetration (track)
		bool interpenetration_;

		//Cooridnation
		bool coordination_;

	public:
		Analyse();
		~Analyse();
		//Init:
		void plug(Sample&,Cell&,Interactions&);
		void initfolder(std::string folder) {folder_ = folder;}
		void init(std::ifstream&);

		//Main: tic at wich analysis is carried
		void analyse(int,double,bool);
		//Init analyse types to false by default
		void allFalse();
		//Reset output files
		void cleanFiles();

		//Analyses:
		void printSample(int);
		void writePS(const std::string,const std::vector<Particle>&);
		//Energy in the system
		void computeEnergy(const double);
		//Cell strain tensor
		void strain(const double) const;
		//Internal and external stress
		void stress(const double) const;
		//Solid fraction
		void compacity(const double) const;

		//Profiles
		void ProfileVelocity(const double) const;

		//Interpenetration
		void Interpenetration(const double) const;
		//Fabric
		void fabric(const double) const;
		void Z(const double) const;
		void printh(const double) const;



};


#endif
