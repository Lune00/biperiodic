#ifndef hpp_ANALYSE_hpp
#define hpp_ANALYSE_hpp

#include<vector>
#include<fstream>


class Sample;
class Cell;
class Particle;


//Ensuite faire une app aussi qui fait l'analyse en relecture
//des fichiers issus de la simulation. Possiblité de découpler
//les deux

class Analyse{

	private:

		Sample* spl_;
		Cell* cell_;
		std::string folder_;
		//Thickness (in Rmax) of periodic band around
		//sample visualisation
		double e_;

	public:
		Analyse();
		~Analyse();
		void printSample(int);
		void plug(Sample&,Cell&);
		void initfolder(std::string folder) {folder_ = folder;}
		void init(std::ifstream&);
		void analyse(int,double);
		void writePS(const std::string,const std::vector<Particle>&);

};


#endif
