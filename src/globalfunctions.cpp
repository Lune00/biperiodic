#include"globalfunctions.hpp"

using namespace std;

//Format output to current tic step k
string formatfile(const string folder, const string file, const int k){
	std::stringstream ss;
	ss << std::setw(wi)<< std::setfill('0') << k ;
	string ffile = folder + "/" + ss.str() + file;
	return ffile;
}

//A test:
Vecteur absVelocity(const Particle& p, const Tensor2x2& h, const Tensor2x2& hinv, const Tensor2x2 hd){

	return (h * p.getV() + hd * p.getR());

}

