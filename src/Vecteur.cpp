#include"Vecteur.hpp"

using namespace std;

void Vecteur::load(ifstream& is){
	is >> x_ >> y_ ;
}

void Vecteur::write(ofstream& os)const {
	os << x_ <<" "<<y_;
}

