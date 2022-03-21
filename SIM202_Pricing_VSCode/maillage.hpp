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
    void operator = (const Numeros & ); //operateur d'affectation
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

    //constructeur par copie d'un maillage
    Maillage(Maillage & );
    //fonction de base maillant le carré unité avec m*n cases (n lignes, m colonnes)
    void maille_carre_unite(int n, int m);
    //fonction fixant la composante valeur du maillage
    void valeur(const vector<double> & );
    //constructeur de maillage carré appelant maille_carre_unite avec m*n cases (n lignes, m colonnes),  et fixant les valeurs initiales a Q
    Maillage(int n, int m, const vector<double> & Q);
    //fonction effectuant la transformation affine Ax + t
    Maillage& tf_affine(const vector<double> & A, const vector<double> & t);
    //fonction permettant de mailler le rectangle [a,b]*[c,d] avec m*n cases (n lignes, m colonnes)
    void maille_rectangle(double a, double b, double c, double d, int n, int m);
    //constructeur de maillage rectangulaire [a,b]*[c,d] avec m*n cases appelant maille_rectangle, et fixant les valeurs initiales a Q
    Maillage(double a, double b, double c, double d, int n, int m, const vector<double> & Q) ;

    //ASSEMBLAGE DE LA MATRICE B
    Matrice mat_B();
    //ASSEMBLAGE DE LA MATRICE MASSE
    Matrice mat_M();
    //ASSEMBLAGE DE LA MATRICE K
    Matrice mat_K();
    //MATRICE D
    Matrice mat_D();

    //INVERSION D'UN SYSTEME LINEAIRE PAR FACTORISATION LU
    vector<double>& Gauss() ;

    //calcul du maillage ulterieur à partir d'un maillage existant
    Maillage& resolution(vector<double> & khi, double r) ;

    //fonction affichant la liste des sommets, leur valeur et les numéros des triangles
    void affiche() const;
    //export du maillage dans un fichier
    void saveToFile(const char *fn) const;
};
//fonction externes

//CALCUL MATRICE SOURCE A
Matrice source_A(Point *P) ;
//CALCUL DE LA SOURCE V
vector<double> source_V(Point *P) ;
//CALCUL DE MATRICE ELEMENTAIRE M
Matrice mat_elem_M (Point P1,Point P2,Point P3) ;
//CALCUL MATRICE ELEMENTAIRE K
Matrice mat_elem_K (Point * P1,Point * P2,Point * P3) ;
//CALCUL MATRICE ELEMENTAIRE B
Matrice mat_elem_B(Point * P1,Point * P2,Point * P3)) ;



//---------------------------------------------------------------------------
//     classe Video  (2D)
//---------------------------------------------------------------------------

class Video
{public:
    Maillage maille;                 //maillage de base
    unsigned long M;                 //nombre d'images
    list<vector<double>> images;     //liste des vecteurs des valeurs du maillage a chaque instant
    
    //constructeur par valeur pour M instant sur un rectangle [a,b]*[c,d] avec m*n cases, et fixant les valeurs initiales a Q
    Video(unsigned long M, double a, double b, double c, double d, int n, int m, const vector<double> & Q);
    //calcul du processus itératif pour un probleme donne
    void resolution(vector<double> & khi, double r);
};


#endif /* maillage_hpp */
