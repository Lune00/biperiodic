#ifndef Algo_hpp
#define Algo_hpp

class Cell;
class Particle;

class Algo{
	private:
		double dt_;

	public:
		Algo(){};
		~Algo(){};
		void verletalgo(Cell&,std::vector<Particle>&,double);

};

#endif
