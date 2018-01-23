#include"Probe.hpp"
#include"Vecteur.hpp"
#include"Tenseur.hpp"
#include"Particle.hpp"

using namespace std;

bool Probe::containCenter (const Particle * p, const Tensor2x2& h) const{
	Vecteur r = h * p->getR() ;
	double dx = fabs(p->getR().getx() - xc_);
	double dy = fabs(p->getR().gety() - yc_);
	return ( dx < hl_ && dy < hh_ );
}
