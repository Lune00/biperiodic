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

		Vecteur r_;
		Vecteur n_;
		Vecteur t_;
		Vecteur f_;

		//Interpenetration:
		double dn_;
		double dt_;

		bool isActif_;

	public:
		Contact(){ i_ = NULL; j_ =NULL ; isActif_ = false;}
		Contact(Particle* i, Particle* j);
		~Contact(){};

		bool isActif() const { return isActif_;}
		void activate() { isActif_ = true;}
		void Frame(Tensor2x2&);
		void write(std::ofstream&) const;
};



#endif
