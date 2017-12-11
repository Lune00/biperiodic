#include"Sample.hpp"
#include"Cell.hpp"

typedef std::vector<Particle>::iterator spit;

using namespace std;

Sample::Sample(){
	sampleIsLoaded_ = false ;
	sampleIsFilled_ = false ;
	rhodefined_ = false ;
	M_ = 0.;
	cell_ = NULL;
	fsampleIni_ = string();
	folder_ = string();

	//Un interet d'avoir absolu mis a part representation???
	//Avec reduced et h on peut recreer absolu quand on veut
	//A voir
	fsample_ = "reduced.txt";
	fsampleA_ = "absolute.txt";
}

Sample::~Sample(){
}


void Sample::plugtoCell(Cell& cell){
	cell_ = &cell;
}

void Sample::init(ifstream& is){

	string token;
	is >> token;
	while(is){
		if(token=="includeFile") is >> fsampleIni_;
		if(token=="rho") {
			rhodefined_= true;
			is >> rho_;
		}
		is >> token ;
	}

	loadSample();
	attributeMass();
	setminmax();
}

void Sample::loadSample(){

	ifstream is(fsampleIni_.c_str());
	if(!is){
		cerr<<"Sample::loadSample() : can not open file."<<endl;
		return;
	}
	else{
		if(isEmptySampleFile(is)) {
			cerr<<fsampleIni_<<" est vide!"<<endl;
			return ;
		}
		else
		{
			//Pb pour la derniere lecture...
			string token;
			while(is){
				if(is.eof()) break;
				Particle P(is);
				spl_.push_back(P);
			}
			//Solution provisoire:
			spl_.pop_back();

			cout<<"Nombre de particules: "<<spl_.size()<<endl;
			sampleIsLoaded_ = true ;
		}

	}
}


void Sample::initfolder(string folder){
	folder_ = folder ;
}

//Coordonnees reduites
void Sample::write(int k) const{

	string filename = formatfile( folder_, fsample_, k );
	ofstream file(filename.c_str());

	for(std::vector<Particle>::const_iterator it = spl_.begin(); it!= spl_.end(); it++){
		it->write(file);
	}
}

//Coordonnees absolues
void Sample::writeAbsolute(int k) const{
	Tensor2x2 h = cell_->geth();

	string filename = formatfile( folder_, fsampleA_, k );
	ofstream file(filename.c_str());

	for(std::vector<Particle>::const_iterator it = spl_.begin(); it!= spl_.end(); it++){
		it->write(file,h);
	}
}

//Coordonnees absolues: ecrit egalemement les particules periodiques dans une epaisseur e (distance en diametre max aux bords)a autour de la cellule
//Comment distinguer les particules periodiques des non periodiques sur l'image?
//on ecrit une autre fonction write Particule qui met a la fin un indice pour discriminer les deux
vector<Particle> Sample::getimages(double e) const{

	vector<Particle> images;

	Tensor2x2 h = cell_->geth();
	Tensor2x2 hinv = h.getInverse();

	e *= getrmax();
	//Surely a clever way to do it...
	//Tmp: build imaginary particules at the center of each boundary 
	Vecteur right(1.,0.5);
	Vecteur left(0.,0.5);
	Vecteur top(0.5,1.);
	Vecteur bottom(0.5,0.);

	//Turn into absolute coordinates:
	left = h * left;
	right = h * right;
	top = h * top;
	bottom = h * bottom;

	for(std::vector<Particle>::const_iterator it = spl_.begin(); it!= spl_.end(); it++){
		cout<<"Particule "<<it->getId()<<endl;
		Vecteur rabs = h * it->getR() ;
		//Distances to imaginary particules on the boundary
		//By definition should never be negative
		double dxleft = (rabs.getx() - left.getx() );
		double dxright = (right.getx() - rabs.getx() );
		double dytop = (top.gety() - rabs.gety() );
		double dybottom = (rabs.gety() -bottom.gety() );

		bool nearleftB =  dxleft < e;
		bool nearrightB =  dxright < e;
		bool neartopB =  dytop < e;
		bool nearbottomB = dybottom < e;

		cout<<"to right: "<<dxright<<" "<<e<<endl;
		cout<<"to left: "<<dxleft<<" "<<e<<endl;
		if( nearleftB){
			cout<<"Left"<<endl;
			Particle a = *(it);
			Vecteur dr(dxleft,0.); 
			dr = hinv * dr ;
			a.setrx(dr.getx()-1.);
			images.push_back(a);
		}
		if ( nearrightB){
			cout<<"Right"<<endl;
			Particle a = *(it);
			Vecteur dr(dxright,0.); 
			dr = hinv * dr ;
			a.setrx(1.-dr.getx());
			images.push_back(a);
		}
		if (neartopB){
			cout<<"Top"<<endl;
			Particle a = *(it);
			Vecteur dr(0.,dytop); 
			dr = hinv * dr ;
			a.setry(1. - dr.gety());
			images.push_back(a);
		}
		if (nearbottomB){
			cout<<"Bottom"<<endl;
			Particle a = *(it);
			Vecteur dr(0.,dybottom); 
			dr = hinv * dr ;
			a.setry(dr.gety()-1.);
			images.push_back(a);
		}
	}
	return images;

}

//Fill particles with mass and inertia
void Sample::attributeMass(){
	if(rhodefined_) {
		for(std::vector<Particle>::iterator it = spl_.begin(); it!= spl_.end(); it++){
			double s = it->getR() * it->getR() * M_PI;
			double m = s * rho_ ;
			it->setMasse(m);
			M_ += m;
		}
		sampleIsFilled_ = true ;
	}
}

//Calcul rmin,rmax,xmin,xmax,ymin,ymax de l'echantillon
//Sert a initialiser geometrie cellule au debut simulation
//Penser a recalculer ca au cours de la deformation?...
void Sample::setminmax(){

	if(!sampleIsLoaded_ || spl_.size() == 0){
		cerr<<"Sample::setminmax() Charger l'echantillon avant."<<endl;
		return;
	}

	xmin_ = spl_[0].getx() - spl_[0].getRadius();
	xmax_ = spl_[0].getx() + spl_[0].getRadius();
	ymin_ = spl_[0].gety() - spl_[0].getRadius();
	ymax_ = spl_[0].gety() + spl_[0].getRadius();

	rmin_ = spl_[0].getRadius();
	rmax_ = spl_[0].getRadius();

	for(spit it= spl_.begin(); it != spl_.end(); it++){
		rmax_ = max(rmax_, it->getRadius());
		rmin_ = min(rmin_, it->getRadius());
		xmin_ = min(xmin_, it->getx() - it->getRadius());
		xmax_ = max(xmax_, it->getx() + it->getRadius());
		ymin_ = min(ymin_, it->gety() - it->getRadius());
		ymax_ = max(ymax_, it->gety() + it->getRadius());
	}
	//cout<<"xmin = "<<xmin_<<endl;
	//cout<<"xmax = "<<xmax_<<endl;
	//cout<<"ymin = "<<ymin_<<endl;
	//cout<<"ymax = "<<ymax_<<endl;
}

bool Sample::initcheck() {
	if( sampleIsLoaded_ && sampleIsFilled_ ) return true;
	else return false;
}

bool Sample::isEmptySampleFile(ifstream& is){
	return is.peek() == std::ifstream::traits_type::eof();
}

//Shift and reduce
void Sample::initReducedCoordinates(Cell& cell){
	double Lx = cell.getLx();
	double Ly = cell.getLy();
	cout<<"Lx = "<<Lx<<endl;
	cout<<"Ly = "<<Ly<<endl;
	for(spit it= spl_.begin(); it != spl_.end(); it++){
		it->setr( (it->getx()-xmin_)/Lx, (it->gety()-ymin_)/Ly);
	}
}

