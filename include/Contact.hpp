#ifndef hpp_CONTACT_hpp
#define hpp_CONTACT_hpp

#include<iostream>
#include"Vecteur.hpp"

class Particle;
class Tensor2x2;

class Contact{

	private:
		//Should be Particle * const i_/j_
		Particle * i_;
		Particle * j_;
		//Contact position
		Vecteur r_;
		//Relative velocity at contact
		Vecteur v_;
		//Normal to the contact
		Vecteur n_;
		//Tangent to the contact
		Vecteur t_;
		//Force at the contact
		Vecteur f_;

		//Interpenetration:
		double dn_;
		double dt_;

		bool isActif_;
		//Avoid truncature errors
		//Used by frame for threshold for interpenetration
		static const double tolerance_;

	public:
		Contact(){ i_ = NULL; j_ =NULL ; isActif_ = false; }
		Contact(Particle* i, Particle* j);
		~Contact(){};

		bool isActif() const { return isActif_;}
		void activate() { isActif_ = true;}
		void Frame(Tensor2x2&);
		void write(std::ofstream&) const;
		void updateRelativeVelocities();
		void computeForce();
		void updateAccelerations();
};



#endif
