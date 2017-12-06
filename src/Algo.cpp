#include"Algo.hpp"
#include"Cell.hpp"
#include"Particle.hpp"
#include"Sample.hpp"

using namespace std;


void Algo::init(ifstream& is){
	string token;
	is >> token;
	while(is){
		if(token=="dt") is >> dt_;
		if(token=="ns") is >> ns_;
		if(token=="}") break;
		is >> token;
	}
}

//A adapter comme on souhaite en terme de tests
bool Algo::initcheck(){
	if (dt_ < 1. && ns_ != 0 ) return true;
	return false;
}


void Algo::run(Cell& cell, Sample& spl){

	double tfinal = ns_ * dt_ ;
	double t=0.;

	ofstream test("samplet.txt");
	ofstream testcell("cellt.txt");
	ofstream strain("strain.txt");
	cout<<"Simulation:"<<endl;
	cout<<"dt = "<<dt_<<endl;
	cout<<"ns = "<<ns_<<endl;
	cout<<"tfinal = "<<tfinal<<endl;
	while(t<tfinal){

		//Time step:
		verletalgo2(cell,spl);

		//Analyse, ecriture:
		spl.writeAbsolute(test);
		cell.write(testcell,t);
		cell.writeStrainTensor(strain,t);

		t+=dt_;
	}
	test.close();

}


//Pour l'instant on l'implemente de maniere naive
//et non optimale (ecriture condensee a l'aide tenseurs/vecteurs)
//On verra apres comment rendre ca plus compacte
void Algo::verletalgo2(Cell& cell,Sample& spl){

  double dt2_2 = 0.5 * dt_ * dt_ ;
  double dt_2 = 0.5 * dt_ ;

  vector<Particle>* ps = spl.getSample();

  for(std::vector<Particle>::iterator it = ps->begin(); it != ps->end(); it++){

	  Vecteur a = it->geta();
	  Vecteur r = it->getR();
	  Vecteur v = it->getV();
	  //Rotations, to add
	  r = r + v * dt_ + a * dt2_2;
	  v = v + a * dt_2;
	  //update position et acceleration debut pas temps
	  it->setRV(r,v);

  }

  //Periodicite en position des particules
  //Peut etre a bouger dans Sample plutot
  //Ca me parait plus etre un taff de sample de modifier les positions
  cell.PeriodicBoundaries2(ps);

  //Integration du mvt de la cellule:
  Tensor2x2 h = cell.geth();
  Tensor2x2 hd = cell.gethd();
  Tensor2x2 hdd = cell.gethdd();

  //h = h + hd * dt_ + hdd * dt2_2 ;
  //Controle en force ou controle en vitesse
  //Temporaire:
  double hxx, hxy , hyx , hyy ;
  double hdxx, hdxy , hdyx , hdyy ;

  if(cell.getControlxx() == 'f' ) {
	  hxx = h.getxx() + hd.getxx() * dt_ + hdd.getxx() * dt2_2 ;
	  hdxx = hd.getxx() + hdd.getxx() * dt_2 ;
  }
  else{
	  //Vitesse imposee reste la meme (hdd, definit par Ld au debut)
	  //hdd par definition nulle si control en vitesse sur hdd
	  hxx = h.getxx() + hd.getxx() * dt_ ; 
	  hdxx = hd.getxx();
  }

  if(cell.getControlxy() == 'f' ) {
	  hxy = h.getxy() + hd.getxy() * dt_  + hdd.getxy() * dt2_2 ;
	  hdxy = hd.getxy() + hdd.getxy() * dt_2 ;
  }
  else{
	  hxy = h.getxy() + hd.getxy() * dt_ ; 
	  hdxy = hd.getxy();
  }

  if(cell.getControlyx() == 'f' ) {
	  hyx = h.getyx() + hd.getyx() * dt_  + hdd.getyx() * dt2_2 ;
	  hdyx = hd.getyx() + hdd.getyx() * dt_2 ;
  }
  else{
	  hyx = h.getyx() + hd.getyx() * dt_ ; 
	  hdyx = hd.getyx();
  }

  if(cell.getControlyy() == 'f' ) {
	  hyy = h.getyy() + hd.getyy() * dt_  + hdd.getyy() * dt2_2 ;
	  hdyy = hd.getyy() + hdd.getyy() * dt_2 ;
  }
  else{
	  hyy = h.getyy() + hd.getyy() * dt_ ; 
	  hdyy = hd.getyy();
  }

  //Set new h
  h.set(hxx,hxy,hyx,hyy);
  hd.set(hdxx,hdxy,hdyx,hdyy);
  cell.update(h,hd);


  //Calcul des forces entre particules a la nouvelle position fin du pas de temps


  //Calcul du tenseur de contraintes internes: sigma_int

  //Calcul des vitesses a la fin du pas de temps:

  //Etape imcomprise: acceleration(i) = hinverse(t+dt) * acceleration(i)
  //Fin pas de temps vitesse
  //On retransforme la vitesse en coordonnee reduite
  Tensor2x2 hinv = h.getInverse();
  for(std::vector<Particle>::iterator it = ps->begin(); it != ps->end(); it++){

	  Vecteur a = it->geta();
	  Vecteur v = it->getV();
	  a = hinv * a ;
	  v = v + a * dt_2 ;
	  it->setV(v);
	  //if(it->getId() == 0) it->affiche();
  }

  //Apply stress_ext: si controle en force
  //stress ext est egal a celui impose
  //Sinon il est egal a -stressint et acc nulle
  cell.ApplyBC();
  double V = cell.getVolume();
  double mh = cell.getMasse();
  Tensor2x2 TotalStress = cell.getStressInt() + cell.getStressExt();
  hdd = hinv * (V/mh) * (TotalStress);
  //On a l'accleration en fin de pas, on peut integrer l'espace en vitesse a la fin du pas
  hd = hd + hdd * dt_2 ;
  cell.updatehd(hd);

  //Cell deformation
  cell.CalculStrainTensor();
}
