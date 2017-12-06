#include<iomanip>
#include<cmath>
#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<string>
#include<vector>

using namespace std;

class Particule{

  private:
    double x_, y_,r_;
    double rot_, vx_, vy_, vrot_;
    string type_;
    unsigned int group_;
  public:
    Particule();
    ~Particule();
    double getr() const {return r_;};
    void setr(double r) { r_ = r ;};
    double getx() const {return x_;};
    void setx(double x) { x_ = x ;};
    double gety() const {return  y_;};
    void sety(double y) { y_ = y ;};
    double getrot() const {return rot_;};
    double getvrot() const {return vrot_;};
    double getvx()  const {return vx_;};
    double getvy() const {return vy_;};
    std::string gettype() const {return type_;};
    double getgroup()  const {return group_;};
};


Particule::Particule()
{
 type_ = "disk"; 
  group_ = 0 ; 
  x_ = 0. ;
  y_ = 0. ;
  r_ = 0. ;
  vx_ = 0. ;
  vy_ = 0. ;
  rot_ = 0. ;
  vrot_ = 0. ;
};

Particule::~Particule() {};



double rand(double min,double max){
  return ((max-min)*((double)rand()/(double)RAND_MAX)+min);
}

double powerlaw(double r1, double r2){
  double y = rand(0.,1.);
  double tmp = (1./r2-1./r1)*y + 1./r1;
  return 1./tmp;
}


int main(){


  //Paramètres echantillon et de la rugosite de la paroi

  //Echantillon:

  unsigned int nfree = 9 ;
  unsigned int ncouches = 3 ;
  double const r1 = 1. ;
  double const r2 = r1 ;
  double const rmean = 0.5 * ( r1 + r2 ) ;

  /* * * * * * * * * *  */

  unsigned int nslice = nfree / ncouches ;


  if(nfree == 0) {cout<<" Entrez un nombre de particules différent de 0."<<endl; return 0 ; }
  if(nslice ==0) {cout<<" Le paramètre nslice a une valeur non acceptable."<<endl; return 0; }
//  if(u < 0. || u * 2 * rparoi > 2 * r1) {cout<<"L'espacement entre particules de la paroi est trop grand (trous) ou trop faible (overlaping)."<<endl; return 0;}


  std::vector<Particule> sample(nfree);

  cout<<"--- Echantillon ---"<<endl;
  cout<<"Nombre de particules libres : "<<sample.size()<<endl;

  double xmin,xmax,ymin,ymax;

  xmin = 0. ;
  ymin = 0. ;
  // Initiatlisation des rayons:
  for(std::vector<Particule>::iterator it = sample.begin() ; it!= sample.end(); it++){
	  it->setr( powerlaw(r1,r2));
  }
  // Initialisation des positions sur la grille:
  unsigned int k = 0 ;
  unsigned int j = 0 ;
  for(std::vector<Particule>::iterator it = sample.begin() ; it!= sample.end(); it++){
	  double x, y;
	  double eps = 0. ; // rand(0.,1.) * r1 * 0.5   ;

	  if(it == sample.begin() || j == 0) 
	  { 
		  x = xmin ;
		  y = ymin + k * 2 * r2 + eps ;
	  }
	  else
	  {
		  x = (it - 1)->getx() + (it - 1)->getr() + it->getr() + eps ;
		  y = ymin + k * 2 * r2 ;
	  }

	  j++;

	  if( j == nslice ) { k++; j = 0 ; }

	  it->setx(x);
	  it->sety(y);
  }


  //Ecrture du fichier packing0.spl

  ofstream myFile ("packing0.spl",ios::out);
  int n = 0 ;
  for(std::vector<Particule>::iterator it = sample.begin() ; it!= sample.end(); it++){
	  myFile<<n<<" "<<it->getr()<<" "<<it->getx()<<" "<<it->gety()<<" "<<it->getvx()<<" "<<it->getvy()<<" "<<it->getrot()<<" "<<it->getvrot()<<endl;
	  n++;
  }
  myFile.close();

  ofstream plotsample("sample-gnuplot.txt");

  for(std::vector<Particule>::iterator it = sample.begin() ; it!= sample.end(); it++){
	  plotsample<<it->getr()<<" "<<it->getx()<<" "<<it->gety()<<" "<<it->getrot()<<" "<<it->getvx()<<" "<<it->getvy()<<" "<<it->getvrot()<<endl;
  }
  plotsample.close();
  return 0;
}
