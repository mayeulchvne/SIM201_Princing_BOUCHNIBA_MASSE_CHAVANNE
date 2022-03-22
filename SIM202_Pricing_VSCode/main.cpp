//
//  main.cpp
//  Maillage
//
//  Created by Mayeul CHAVANNE on 01/02/2022.
//

#include <vector>
#include <iostream>
#include <fstream>
#include "maillage.hpp"

int main( int argc, char * argv[] ) {
  //dimension du maillage
  double a(0) ;
  double b(1);
  double c(0);
  double d(1);

  //nombre de ligne et de colonnes calculs pour un pas h = max( (b-a)/m , (c-d)/n )
  // double h;
  int n(50);
  int m(50);

  //donnée du modèle de BS
  double K;
  vector<double> khi(4 , 0.);
  khi[0] =  0.04; khi[1] =  -0.024; khi[3] =  0.04; khi[2] =  -0.024; 
  double r(0.05);

  //nombre d'instants & pas de temps
  int M(100);
  double delta_t(1./M);

  //création du vecteur initiale Q0
  vector<double> Q0((n+1)*(m+1), 0.);

  //creation du maillage
  Maillage maillage(a, b , c, d, n, m, Q0);

  //remplissage du vecteur initiale Q0
  Point P(0.,0.);
  double tmp = 0.;
  for(int q =0 ; q < (n+1)*(m+1) ; q++) {
    P = maillage.sommets[q];
    tmp = P.x + P.y - K;
    if(tmp >= 0) Q0[q] = tmp;
    else Q0[q] = 0.;
  }

  //creation d'une video 
  Video video(M, a, b , c, d, n, m, Q0);

  //resolution du probleme temporel
  video.resolution(khi, r, delta_t);

  //écriture des resultats dans un fichier texte
  video.saveToFile("maillage1.txt");


  return 0;
}
