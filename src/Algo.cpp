#include"Algo.hpp"
#include"Cell.hpp"
#include"Particle.hpp"

using namespace std;


void Algo::run(){


}

//Le couplage passe par le tenseur de contrainte entre cellule et particules

//Update particules, calcul tenseur contrainte interne, update cellule

//On integre l'equation de la dynamique pour l'espace avec un algo de verlet
//Integre la dynamique de l'espace, de la cellule au temps t, a besoin tenseur contraintes
//Il y a couplage entre debut et fin du pas de temps
//Il faut integrer les deux en meme temps
//void Algo::verletalgo(Cell& cell,std::vector<Particle>& ps){
//
//	cout<<"Hello"<<endl;
//	//Attention: variable globale algo hinv, stress ext et int
//	//Init:
//	Tensor2x2 h = cell.geth();
//	Tensor2x2 hd = cell.gethd();
//	Tensor2x2 hinv = h.getInverse();
//
//	Tensor2x2 hdd0;
//	Tensor2x2 hdd1 = hdd0;
//
//	double mh=cell.getMasse();
//	double V0=cell.getVolume();
//
//	//Impose BC() en force: si une direction est controlee en vitesse, on impose hdd0 nulle dans cette dir
//	//Et on calcule le tenseur de contrainte internes a partir de l'equation ac terme de droite nul
//	//Cad que stressext=-stressint dans cette direction
//
//	//Si dir controle en vitesse, stressext dans cette dir = stress int
//	//On calcule les composantes de stressext dans cette direction
//	//Apres on fait le produit matriciel
//	//On recupere le stressint au debut du pas de temps (forces particules au debut du pas)
//	// computeInternalStress()
//	// Pour l'istant il est nul
//
//	//On applique les BCU pour recuperer les contraintes ext
//	cell.ApplyBC();
//	Tensor2x2 TotalStress = cell.getStressInt() +cell.getStressExt();
//
//
//	//Calcul de hdd au debut du pas, hdd0
//	//Eq (24)
//	//On a besoin de stressint au debut du pas de temps
//	hdd0=hinv*(V0/mh)*(TotalStress);
//
//	//Calcul position(h) avec hdd0 a la fin du pas de temps
//	h = h + hd*dt_ + hdd0*0.5*dt_*dt_;
//
//	//Calcul des positions particules fin du pas de temps
//
//	for(std::vector<Particle>::iterator it = ps.begin(); it != ps.end(); it++){
//
//		//Force pas de temps precedent:
//		Vecteur F = it->getForce(); 
//		Vecteur r = it->getR();
//		Vecteur v = it->getV();
//		double m = it->getMasse();
//
//		//Acceleration debut pas de temps
//		Vecteur tmp = hinv * r;
//		tmp = TotalStress * tmp; 
//		tmp = hinv * tmp;
//		tmp = tmp * (V0/mh);
//		//a0 a stocker
//		Vecteur a0 = F * (1./m) + tmp; 
//		r = r + v * dt_ + a0 * dt_ * dt_ * 0.5;
//		//update position et acceleration debut pas temps
//		it->setR(r,a0);
//	}
//
//	//On peut a present recalculer le tenseur des contraintes internes a partir des forces fin du pas de temps
//	//Calculer les forces a partir nouvelle position
//	//MD()
//	//computeInternalStress()
//	//Mettre a jout stressint de cell
//
//	//Volume fin du pas de temps
//	double V1=h.getDet();
//	//Calcul hdd fin du pas, hdd1
//	hinv=h.getInverse();
//
//	//Ensuite on calcule l'accleration de chaque particule fin du pas de temps en utilisant le nouveau h,V
//	//On a les vitesses des particules: position et vitesses particules fin pas de temps
//
//	cell.ApplyBC(); // permet de definir le stressext fin pas de temps
//	//TotalStress fin du pas de temps
//	TotalStress = cell.getStressInt() +cell.getStressExt(); 
//
//	TotalStress.affiche();
//	//Derniere boucle pour calcul des vitesses:
//	//Comment on applique les BCL en vitesse???
//	for(std::vector<Particle>::iterator it = ps.begin(); it != ps.end(); it++){
//
//		//Force recalculee fin du pas de temps:
//		Vecteur F = it->getForce(); 
//		//Position fin pas de temps
//		Vecteur r = it->getR();
//		//Acceleration debut pas de temps
//		Vecteur a0 = it->geta0();
//		double m = it->getMasse();
//
//		//Acceleration fin pas de temps
//		Vecteur tmp = hinv * r;
//		tmp = TotalStress * tmp; 
//		tmp = hinv * tmp;
//		tmp = tmp * (V1/mh);
//		Vecteur a1 = F * (1./m) + tmp; 
//		a1.print();
//		Vecteur v = it->getV() + (a1 + a0) * dt_ * 0.5 ;
//		v.print();
//		cout<<"vitesse:";
//		//update vitesse
//		it->setV(v);
//	}
//	//Fin integration particules
//
//
//	//Acceleration fin du pas: pareil que pour hdd0, on ompose hdd1 nul dans direction controle en vitesse
//	//Il faut prendre le nouveau tenseur contrainte interne et externe en appliqand les BCU fin pas de temps
//
//	hdd1=hinv*(V1/mh)*(TotalStress);
//
//	//Calcul de hd fin du pas de temps
//	//Alors la vitesse ds dir controle en vitesse est bien constante (hdd0+hdd1=0)
//	hd=hd+(hdd0+hdd1)*dt_*0.5;
//
//	//Update Cell:on fait remonter h,hd
//	cell.update(h,hd,dt_);
//}

//Pour l'instant on l'implemente de maniere naive
//et non optimale (ecriture condensee a l'aide tenseurs/vecteurs)
//On verra apres comment rendre ca plus compacte
void Algo::verletalgo2(Cell& cell,std::vector<Particle>& ps){

  double dt2_2 = 0.5 * dt_ * dt_ ;
  double dt_2 = 0.5 * dt_ ;
  //Integre mvt particules:
  for(std::vector<Particle>::iterator it = ps.begin(); it != ps.end(); it++){

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
  //Construire images si particules sortent de la cellule

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
  for(std::vector<Particle>::iterator it = ps.begin(); it != ps.end(); it++){

	  Vecteur a = it->geta();
	  Vecteur v = it->getV();
	  a = hinv * a ;
	  v = v + a * dt_2 ;
	  it->setV(v);
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
