//
//  maillage.hpp
//  Maillage
//
//  Created by Mayeul CHAVANNE on 01/02/2022.
//

#ifndef maillage_hpp
#define maillage_hpp

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <map>
#include <utility>
#include <cstdlib>
#include <ctime>

using namespace std;

void stop(const char * msg);                     // message d'arrêt


//---------------------------------------------------------------------------
//     classe Point  (point 2D)
//---------------------------------------------------------------------------
class Point
{public:
    double x;   //abcisse
    double y;   //ordonnée
    Point(double x0=0,double y0=0):x(x0),y(y0){} //constructeur commun
    Point& operator += (const Point & P) {x+=P.x;y+=P.y;return *this;}
    Point& operator -= (const Point & P) {x-=P.x;y-=P.y;return *this;}
    Point& operator *= (double a){x*=a;y*=a;return *this;}
    Point& operator /= (double a){x/=a;y/=a;return *this;}           //pas de test de division par 0!
    Point& tf_affine(const vector<double> &,const vector<double> &); //tf affine d'un point
};
//fonctions externes
Point operator + (const Point& ,const Point&);
Point operator - (const Point& ,const Point&);
Point operator * (const Point& ,double a);
Point operator * (double a,const Point & );
Point operator / (const Point& ,double a);
bool operator == (const Point& ,const Point& );
bool operator != (const Point& ,const Point& );
bool operator<(const Point&, const Point&);
ostream & operator <<(ostream &, const Point &);

//---------------------------------------------------------------------------
//     classe Numeros (triangles 2D)
//---------------------------------------------------------------------------

class Numeros : public vector<int>
{public :
    Numeros(int i1=0, int i2=0, int i3=0); //constructeur
    ~Numeros(){} //destructeur
    int operator() (const int i) const; //operateur d'acces en lecture
    int & operator() (const int i); //operateur d'acces en ecriture
};

ostream& operator <<(ostream& out, const Numeros & N);
//---------------------------------------------------------------------------
//     classe Maillage  (2D)
//---------------------------------------------------------------------------

class Maillage
{public:
    vector<Point> sommets;          //vecteur des sommets
    vector<double> valeurs;         //vecteur des valeurs prise par le maillage à chaque sommet
    list<Numeros> triangles;        //liste des numéros des sommets des triangles
    
    //fonction de base maillant le carré unité
    void maille_carre_unite(int n, int m);
    //fonction fixant la composante valeur du maillage
    void valeur(const vector<double> & );
    //constructeur de maillage carré appelant maille_carre_unite
    Maillage(int n, int m, const vector<double> & Q);
    //fonction effectuant la transformation affine Ax + t
    Maillage& tf_affine(const vector<double> & A, const vector<double> & t);
    //fonction permettant de mailler le rectangle [a,b]*[c,d] avec m*n cases
    void maille_rectangle(double a, double b, double c, double d, int n, int m);
    //constructeur de maillage rectangulaire [a,b]*[c,d] avec m*n cases appelant maille_rectangle
    Maillage(double a, double b, double c, double d, int n, int m, const vector<double> & Q) ;
    //calcul du maillage ulterieur à partir d'un maillage existant
    void resolution(vector<double> & khi, double r) ;
    //fonction affichant la liste des sommets, leur valeur et les numéros des triangles
    void affiche() const;
    //export du maillage dans un fichier
    void saveToFile(const char *fn) const;
};

//---------------------------------------------------------------------------
//     classe Video  (2D)
//---------------------------------------------------------------------------

class Video
{public:
    unsigned long M;                //nombre d'instants et d'images
    vector<Maillage> images;        //vecteur des maillages
    list<double> instants;          //liste des instants correspondants
    
    //constructeur par valeur pour M instant sur un rectangle [a,b]*[c,d] avec m*n cases
    Video(unsigned long M, double a, double b, double c, double d, int n, int m);
    //calcul du processus itératif pour un probleme donne
    void resolution(vector<double> & khi, double r);
};


#endif /* maillage_hpp */
