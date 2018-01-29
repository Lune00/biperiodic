#ifndef H_PROBE_H
#define H_PROBE_H

class Particle;
class Tensor2x2;
class Vecteur;

//Rectangular probe along x direction
//On pourrait faire une probe qui map bien la cell deformee
//Chaque probe serait une petite cellule avec un tenseur qui decrit
//ses dimensions et sa geometrie

//Pour l'instant on reste a un proto de base avec des probes
//rectangulaires qui slicent l'empilement suivant x (toute la largeur)

class Probe{

	private:
		//Center coordinates (real)
		double yc_;
		//Half height
		double hh_;

	public:
		Probe(): yc_(0.), hh_(0.){}
		Probe(double y, double hh): yc_(y), hh_(hh){}
		~Probe(){};
		bool containCenter (const Particle&) const;
		bool intersection(const Particle&) const;
		//double area () const { return 4.* hh_ * hl_ ;} ;
		double gety() const { return yc_;}

};


#endif
