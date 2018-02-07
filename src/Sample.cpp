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
	load_sample_ = false;
	//Un interet d'avoir absolu mis a part representation???
	//Avec reduced et h on peut recreer absolu quand on veut
	//A voir
	fsample_ = "reduced.txt";
	fsampleA_ = "absolute.txt";
	filetoload_ = 0 ;
	starting_ = 0 ;
	ofstream os("track.txt");
	os.close();
}

Sample::~Sample(){
}


void Sample::plugtoCell(Cell& cell){
	cell_ = &cell;
}


//Build cell : build a new cell from sample
//Load sample/cell : load sample in reduced coord + associated cell
//Load from integer (number of the file) and look into
//sample/ and cell/ folder. If not found one of two abort.
//Cell is initialised anyway, and will be reset if needed.
//ex: build packing0.spl -> init cell from sample (in abs coordinates)
//    load 10 -> load sample/00010reduced.txt cell/00010cell.txt if found
void Sample::init(ifstream& is){

	string token;
	is >> token;
	while(is){
		if(token=="build") {
			is >> fsampleIni_;
			load_sample_ = false;
		}
		if(token=="load"){
			is >> filetoload_;
			is >> starting_ ;
			load_sample_ = true ;
		}
		if(token=="rho") {
			rhodefined_= true;
			is >> rho_;
		}
		if(token=="}") break;
		is >> token ;
	}

	loadSample();
}

//Load sample, set mass to each particle, and compute min/max
void Sample::loadSample(){

	spl_.clear();

	//Build fsampleIni_ from filetoload_
	if(load_sample_) {
		fsampleIni_ = formatfile( folder_, fsample_, filetoload_ );
	}
	//else it is already initialised

	ifstream is(fsampleIni_.c_str());

	cerr<<"Sample to load : "<<fsampleIni_<<endl;
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
			//Solution provisoire (!):
			spl_.pop_back();

			cout<<"Nombre de particules: "<<spl_.size()<<endl;
			sampleIsLoaded_ = true ;
		}

	}

	attributeMass();
	setminmax();
}


void Sample::initfolder(string folder){
	folder_ = folder ;
}

//Coordonnees reduites: allows reloading
void Sample::write(int k) const{

	string filename = formatfile( folder_, fsample_, k );
	ofstream file(filename.c_str());
	file.precision(12);

	for(std::vector<Particle>::const_iterator it = spl_.begin(); it!= spl_.end(); it++){
		it->write(file);
	}
}

//Only for analysis: 
//Coordonnees absolues
void Sample::writeAbsolute(int k) const{
	Tensor2x2 h = cell_->geth();
	Tensor2x2 hd = cell_->gethd();

	string filename = formatfile( folder_, fsampleA_, k );
	ofstream file(filename.c_str());

	for(std::vector<Particle>::const_iterator it = spl_.begin(); it!= spl_.end(); it++){
		it->write(file,h,hd);
	}
}


Vecteur Sample::returnrabs(const Particle& P) const{
	return cell_->geth() * P.getR() ;
}
Vecteur Sample::returnvabs(const Particle& P) const{
	return (cell_->geth() * P.getV() + cell_->gethd() * P.getR() );
}


//TODO: merge these two functions
//And add Ec and Ecrot as member variable of the class Sample

//Translational kinetic energy
double Sample::getTKE() const{

	double Ec = 0.;
	Tensor2x2 h = cell_->geth();
	Tensor2x2 hd = cell_->gethd();

	for(std::vector<Particle>::const_iterator it = spl_.begin(); it!= spl_.end(); it++){
		//Vecteur v = returnvabs(*it);
		Vecteur v = h *it->getV() + hd * it->getR();
		//Vecteur v =  hd * it->getR();
		Ec += 0.5 * it->getMasse() * v.getNorme2();
	}
	return Ec;
}


//Rotational kinetic energy
double Sample::getRKE() const {
	double Ecrot = 0.;
	for(std::vector<Particle>::const_iterator it = spl_.begin(); it!= spl_.end(); it++){
		Ecrot += 0.5 * it->getI() * it->getVrot() * it->getVrot();
	}
	return Ecrot;
}


void Sample::debug(int tic){
	ofstream os("track.txt",ios::app);

	Tensor2x2 h = cell_->geth();
	Tensor2x2 hd = cell_->gethd();
	Tensor2x2 hdd = cell_->gethdd();

	for(std::vector<Particle>::const_iterator it = spl_.begin(); it!= spl_.end(); it++){
		if(it->getId()==13){
			Vecteur sd = it->getV();
			Vecteur v = h * sd + hd * it->getR(); 
			Vecteur a = hdd * it->getR() + hd * sd * 2. + h * it->getA();
			Vecteur r = h * it->getR();
			//os<<tic<<" "<<it->getId()<<" "<<sd.getx()<<" "<<sd.gety()<<" "<<v.getx()<<" "<<v.gety()<<" "<<it->getVrot()<<" "<<a.getx()<<" "<<a.gety()<<" "<<it->getA().getx()<<" "<<it->getA().gety()<<endl;
			os<<tic<<" "<<r.getx()<<" "<<r.gety()<<" "<<it->getVrot()<<endl;
		}
	}
	os.close();
}

void Sample::printSample() const{

	for(std::vector<Particle>::const_iterator it = spl_.begin(); it!= spl_.end(); it++){
		it->print();
	}
}


//Return image particles within a const band width e around the cell in the lab frame. Band width is expressed in maximum radius of the sample
//Coordonnees absolues: ecrit egalemement les particules periodiques dans une epaisseur e (distance en diametre max aux bords)a autour de la cellule
vector<Particle> Sample::getimages(double e) const{

	vector<Particle> images;

	Tensor2x2 h = cell_->geth();
	Tensor2x2 hinv = h.getInverse();

	e *= getrmax();

	for(std::vector<Particle>::const_iterator it = spl_.begin(); it!= spl_.end(); it++){

		//Surely a clever way to do it...
		//Tmp: build imaginary particules at each boundary 
		Vecteur right(1.,it->gety());
		Vecteur left(0.,it->gety());
		Vecteur top(it->getx(),1.);
		Vecteur bottom(it->getx(),0.);

		//Turn into absolute coordinates:
		left = h * left;
		right = h * right;
		top = h * top;
		bottom = h * bottom;

		//cout<<"Particule "<<it->getId()<<endl;
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

		//cout<<"to right: "<<dxright<<" "<<e<<endl;
		//cout<<"to left: "<<dxleft<<" "<<e<<endl;
		if( nearleftB){
			Particle a = *(it);
			//Vecteur dr(dxleft,0.); 
			//dr = hinv * dr ;
			a.addrx(1.);
			images.push_back(a);
		}
		if ( nearrightB){
			Particle a = *(it);
			a.addrx(-1);
			images.push_back(a);
		}
		if (neartopB){
			Particle a = *(it);
			a.addry(-1.);
			images.push_back(a);
		}
		if (nearbottomB){
			Particle a = *(it);
			a.addry(1.);
			images.push_back(a);
		}
		//Diagonals for corners (more pretty)
		if(nearleftB && nearbottomB){
			Particle a = *(it);
			a.addry(1.);
			a.addrx(1.);
			images.push_back(a);
		}
		if(nearleftB && neartopB){
			Particle a = *(it);
			a.addry(-1.);
			a.addrx(1.);
			images.push_back(a);
		}
		if(nearrightB && neartopB){
			Particle a = *(it);
			a.addry(-1.);
			a.addrx(-1.);
			images.push_back(a);
		}
		if(nearrightB && nearbottomB){
			Particle a = *(it);
			a.addry(1.);
			a.addrx(-1.);
			images.push_back(a);
		}
	}
	return images;

}

//Fill particles with mass and inertia
void Sample::attributeMass(){
	if(rhodefined_) {
		cout<<"Set inertia to particles...";
		for(std::vector<Particle>::iterator it = spl_.begin(); it!= spl_.end(); it++){
			double s = it->getRadius() * it->getRadius() * M_PI;
			double m = s * rho_ ;
			it->setInertia(m);
			M_ += m;
		}
		sampleIsFilled_ = true ;
		cout<<"done."<<endl;
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
//Warnging: by default, the new build cell is assumed to be rectangular
void Sample::initReducedCoordinates(){

	double Lx = cell_->geth().getxx();
	double Ly = cell_->geth().getyy();

	for(spit it= spl_.begin(); it != spl_.end(); it++){
		it->setr( (it->getx()-xmin_)/Lx, (it->gety()-ymin_)/Ly);
	}
}

//First step Verlet algo, integration
void Sample::firstStepVerlet(const double dt_){

	double dt2_2 = 0.5 * dt_ * dt_ ;
	double dt_2 = 0.5 * dt_ ;

	for(spit it = spl_.begin(); it != spl_.end(); it++){

		//Positions
		it->updateR(dt_);
		it->updateRot(dt_);
		//First step vitesse: integrate over half step
		it->updateV(dt_2);
		it->updateVrot(dt_2);
		//Set accelerations to 0 for second step:
		it->resetA();
	}

	return;
}


//Update velocity and vrotation at the end of the time step
//Removes mean velocity to ensure that <sd>=0.
//The mean displacement is carried only by the cell (homogenous def)
void Sample::secondStepVerlet(const double dt_) {

	double dt_2 = 0.5 * dt_ ;

	Vecteur vmean;

	for(spit it =spl_.begin(); it != spl_.end(); it++){
		it->updateV(dt_2);
		it->updateVrot(dt_2);
		vmean = vmean + it->getV();
	}

		vmean = vmean / (double)spl_.size();
	
		for(spit it =spl_.begin(); it != spl_.end(); it++){
			it->removevmean(vmean);
		}
}

//Not used
Particle * Sample::getP(int i){

	for(spit it =spl_.begin(); it != spl_.end(); it++){
		if(it->getId()==i) {
			Particle * pp = &(*it) ;
			return pp;
		}
	}
	return NULL ;
}
