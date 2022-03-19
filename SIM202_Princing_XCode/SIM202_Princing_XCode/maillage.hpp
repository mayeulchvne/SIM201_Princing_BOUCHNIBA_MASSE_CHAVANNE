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
//     classe Numeros
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
    vector<Point> sommets;        //liste des sommets comptés une seule fois
    list<Numeros> triangles;        //liste des numéros des sommets des triangles
    
    void maille_carre_unite(int n, int m);    //fonction de base maillant le carré unité
    Maillage(int n, int m) //constructeur de maillage carré appelant maille_carre_unite
    {
        this->maille_carre_unite(n, m);
    }
    void affiche() const;    //fonction affichant la liste des sommets et des numéros des triangles
    Maillage& tf_affine(const vector<double> & A, const vector<double> & t);  //fonction effectuant la transformation affine Ax + t
    void maille_rectangle(double a, double b, double c, double d, int n, int m);  //fonction permettant de mailler le rectangle [a,b]*[c,d] avec m*n cases
    Maillage(double a, double b, double c, double d, int n, int m) //constructeur de maillage rectangulaire [a,b]*[c,d] avec m*n case appelant maille_rectangle
    {
        this->maille_rectangle(a, b, c, d, n, m);
    }
    void saveToFile(const char *fn) const;   //export du fichier dans un maillage
};

#endif /* maillage_hpp */
