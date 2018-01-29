#ifndef h_globalfunctions_h
#define h_globalfunctions_h

#include<string>
#include<sstream>
#include<iomanip>

#include"Tenseur.hpp"
#include"Vecteur.hpp"
#include"Particle.hpp"

static const int wi = 5 ;

//Concatenete a string and a const for file output formating
std::string formatfile(const std::string folder, const std::string file, const int k);

//Test

Vecteur absVelocity(const Particle& p, const Tensor2x2& h, const Tensor2x2& hinv, const Tensor2x2 hd);
#endif
