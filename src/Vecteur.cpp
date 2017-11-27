#include"Vecteur.hpp"
#include"Tenseur.hpp"

using namespace std;

Vecteur Vecteur::productTensor(Tensor2x2& T){
	Vecteur r;
	r.x_ = T.getxx() * this->x_ + T.getxy() * this->y_;
	r.y_ = T.getyx() * this->x_ + T.getyy() * this->y_;
	return r;
}
