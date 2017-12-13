#ifndef hpp_ANALYSE_hpp
#define hpp_ANALYSE_hpp

#include<vector>
#include<fstream>


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

	public:
		Analyse();
		~Analyse();
		//Init:
		void plug(Sample&,Cell&,Interactions&);
		void initfolder(std::string folder) {folder_ = folder;}
		void init(std::ifstream&);

		//Main: tic at wich analysis is carried
		void analyse(int,double);

		//Analyses:
		void printSample(int);
		void writePS(const std::string,const std::vector<Particle>&);
		void computeEnergy(const int);

};


#endif
