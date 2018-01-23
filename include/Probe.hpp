#ifndef H_PROBE_H
#define H_PROBE_H

class Particle;
class Tensor2x2;
class Vecteur;

//Rectangular probe
class Probe{

	private:
		//Center coordinates (real)
		double xc_;
		double yc_;
		//Half height, half length
		double hh_;
		double hl_;

	public:
		Probe(): xc_(0.), yc_(0.), hh_(0.), hl_(0.) {}
		Probe(double x, double y, double hh, double hl): xc_(x), yc_(y), hh_(hh), hl_(hl) {}
		~Probe(){};
		bool containCenter (const Particle*,const Tensor2x2&) const;
		double area () const { return 4.* hh_ * hl_ ;} ;

};


#endif
