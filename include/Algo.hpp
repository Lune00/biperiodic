#ifndef Algo_hpp
#define Algo_hpp
#include<vector>
#include<iostream>


class Cell;
class Particle;

class Algo{
	private:
		double dt_;

	public:
		Algo(double dt): dt_(dt) {};
		~Algo(){};
		void verletalgo(Cell&,std::vector<Particle>&);

};

#endif
