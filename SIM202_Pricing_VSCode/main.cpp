//
//  main.cpp
//  Maillage
//
//  Created by Mayeul CHAVANNE on 01/02/2022.
//

#include <iostream>
#include <fstream>
#include "maillage.hpp"

int main( int argc, char * argv[] )
{ //test de la classe Point
  Point O, P(0,1), Q(1,1), R(1);
  cout<<"O+P+2.*Q-R/2.="<<(O+P+2.*Q-R/2.) <<endl;
  cout<<"O!=P -> "<<(O!=P?"true":"false")<<" O==P -> "<<(O==P?"true":"false")<<endl;
  //test de la classe Numeros
  Numeros ns(1,2,3);
  cout<<"Numeros ns(1,2,3) -> "<<ns<<endl;
  cout<<"ns(2)="<<ns(2)<<endl;
  //test de la classe Maillage
  srand(time(NULL));
  cout<<"Maillage du carre [0,1]x[0,1]"<<endl;
  vector<double> Z(3*5, 0.);
  Maillage mc(4,2, Z);            //maillage carré unité
  mc.affiche();
  mc.saveToFile("maillage1.txt");
  cout<<"Maillage du rectangle [1,2]x[0,2]"<<endl;
  vector<double> W(5*5, 0.);
  Maillage mr(1,2,0,2,4,4, W); //maillage rectangle
  mr.affiche();
  mr.saveToFile("maillage2.txt");
  cout<<"Transformation affine du maillage du carre unite"<<endl;
  vector<double> A(4,0); A[0]=1;A[3]=3;
  vector<double> t(2,0); t[0]=-1;
  mc.tf_affine(A,t); //transformation affine de maillage
  mc.affiche();
//  cout<<"Concatenation de maillage"<<endl;
//  Maillage mf = Maillage(4,2) + mr; //concaténation de maillages
//  mf.affiche();
//  mf.saveToFile("mf.mail");
  return 0;
}
