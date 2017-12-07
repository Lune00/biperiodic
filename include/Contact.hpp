#ifndef hpp_CONTACT_hpp
#define hpp_CONTACT_hpp

#include<iostream>

class Particle;

class Contact{

	private:
		//Should be Particle * const i_/j_
		Particle * i_;
		Particle * j_;

	public:
		Contact(){};
		Contact(Particle* i, Particle* j){ i_ = i ; j_ = j;}
		~Contact(){};

};



#endif
