#include"Probe.hpp"
#include"Vecteur.hpp"
#include"Tenseur.hpp"
#include"Particle.hpp"

using namespace std;

bool Probe::containCenter (const Particle& p ) const{
	double dy = fabs(p.getR().gety() - yc_);
	return ( dy < hh_ );
}

bool Probe::intersection (const Particle& p) const{
	return false;
}
