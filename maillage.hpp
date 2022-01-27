#ifndef maillage.hpp
#define maillage.hpp

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
    Numeros(int i1, int i2, int i3); //constructeur
    int operator() (const int i) const; //operateur d'acces en lecture
    int operator() (const int i); //operateur d'acces en ecriture
};

ostream& operator <<(ostream& out, const Numeros & N);
//---------------------------------------------------------------------------
//     classe Maillage  (2D)
//---------------------------------------------------------------------------
class Maillage
{public:
    vector<Point> sommets;        //liste des sommets comptés une seule fois
    list<Numeros> numelts;        //liste des numéros des sommets des rectangles
    ...
};

#endif
