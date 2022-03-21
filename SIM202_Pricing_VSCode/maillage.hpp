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
//     quelques fonctions utilitaires
//---------------------------------------------------------------------------

// message d'arrêt
void stop(const char * msg);                     

//transforme les element d'un vecteur de vecteur en vecteur
vector<double> vecteur(vector<vector<double> > vec, int i);

//===================================================================================================
//                                     Classe Matrice
//===================================================================================================

class Matrice {

protected:
    virtual double aff(int i, int j) const; 

public:
    int dim; //taille de la matrice carrée
    vector<int> Stack; //matrice des profils
    vector<double> Mat; //matrice des coefficients
    Matrice(){
        this->Stack.resize(0,0.);
        this->Mat.resize(0,0.);
        this->dim=0;
    }
    Matrice(const vector<int> S){
        this->dim=S.size();
        this->Stack=S;
        int h=this->Stack[dim-1];
        this->Mat.resize(h+1,0.);
    }
    Matrice(const Matrice& M){
        this->dim=M.dim;
        this->Stack=M.Stack;
        this->Mat=M.Mat;
    }
    virtual ~Matrice(){}; //Destructeur virtuel de la classe

    //accesseur
    int d() const{return this->dim;}
    //accesseurs virtuels purs
    virtual double operator()(int i, int j) const =0;
    virtual double& operator()(int i, int j) =0;
    //opérations globale
    

};

vector<int> profil(const vector<double>& V);



//===================================================================================================
//                                     Classe héritée Matrice_S
//===================================================================================================

class Matrice_S: public Matrice
{
public:
    Matrice_S()
        :Matrice(){}
    Matrice_S(const vector<int> S)
        :Matrice(S) {}
    Matrice_S(const Matrice& M)
        :Matrice(M) {}
    virtual ~Matrice_S(){};

    virtual double operator()(int i, int j) const ;
    virtual double& operator()(int i, int j) ;
    //opérateurs algébriques 
     Matrice_S& operator+=(const Matrice_S& M);
     Matrice_S& operator-=(const Matrice_S& M);
     Matrice_S& operator*=(const double& a);
     Matrice_S& operator/=(const double& a);
     Matrice_S& operator=(const Matrice_S& M);

};

Matrice_S operator+(const Matrice_S& M,const Matrice_S& H);
Matrice_S operator-(const Matrice_S& M,const Matrice_S& H);
Matrice_S operator*(const Matrice_S& M,const double& a);
Matrice_S operator/(const Matrice_S& M,const double& a);
vector<double> operator*(const Matrice_S& M,const vector<double>& X);
vector<double> operator*(const Matrice_S& M,const vector<int>& X);


//===================================================================================================
//                                     Classe héritée Matrice_PS
//===================================================================================================


class Matrice_PS: public Matrice
{
public:
    Matrice_PS()
        :Matrice(){}
    Matrice_PS(const vector<int> S)
        :Matrice(S) {this->Mat.resize(2*(this->Mat.size())-(this->Stack.size()),0.);}
    Matrice_PS(const Matrice_PS& M)
        :Matrice(M) {}
    Matrice_PS(const Matrice_S& M)
        :Matrice(M){
            auto it1=this->Mat.end();                           //création de l'espace de stockage supplémentaire pour la matrice non symétrique
            vector<double> u(this->Mat);                        //on insert Mat privée des coefficients diagonaux à la fin de Mat 
            auto it2=u.begin();
            auto it3=u.end();
            for(int i=0;i<this->d();i++){u.erase(it2+this->Stack[i]);}
            this->Mat.insert(it1,it2,it3);
        }
    virtual ~Matrice_PS(){};
    //accesseurs
    virtual double operator()(int i, int j) const;
    virtual double& operator()(int i, int j);
    //opérateurs algébriques
     Matrice_PS& operator+=(const Matrice_PS& M);
     Matrice_PS& operator-=(const Matrice_PS& M);
     Matrice_PS& operator*(const Matrice_PS& M);
     Matrice_PS& operator+=(const Matrice_S& M);
     Matrice_PS& operator-=(const Matrice_S& M);
     Matrice_PS& operator*=(const double& a);
     Matrice_PS& operator/=(const double& a);
     Matrice_PS& operator=(const Matrice_PS& M);
     Matrice_PS& operator=(const Matrice_S& M);
     //fonctions membres
    Matrice_PS* LU();
    //Matrice_PS test();
};


Matrice_PS operator+(const Matrice_PS& M,const Matrice_S& H);
Matrice_PS operator+(const Matrice_S& M,const Matrice_PS& H);
Matrice_PS operator+(const Matrice_S& M,const Matrice_PS& H);

Matrice_PS operator-(const Matrice_PS& M,const Matrice_PS& H);
Matrice_PS operator-(const Matrice_PS& M,const Matrice_S& H);
Matrice_PS operator-(const Matrice_S& M,const Matrice_PS& H);

Matrice_PS operator*(const Matrice_PS& M,const double& a);
Matrice_PS operator/(const Matrice_PS& M,const double& a);
vector<double> operator*(Matrice_PS& M,vector<double>& X);
vector<double> operator*(Matrice_PS& M,vector<int>& X);


//Flux de sortis

ostream& operator<<(ostream& os,const Matrice& M);

ostream& operator<<(ostream& os, const vector<double>& v);

ostream& operator<<(ostream& os, const vector<int>& v);


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

    //fonction generant le vecteur des indices des premiers coefficients non nuls de chaque 
    //lignes du profil des matrices associées au maillage m*n cases (n lignes, m colonnes)
    vector<int> indice_profil();
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
    list<vector<double> > images;     //liste des vecteurs des valeurs du maillage a chaque instant
    
    //constructeur par valeur pour M instant sur un rectangle [a,b]*[c,d] avec m*n cases, et fixant les valeurs initiales a Q
    Video(int M, double a, double b, double c, double d, int n, int m, const vector<double> & Q);
    //calcul du processus itératif pour un probleme donne
    void resolution(vector<double> & khi, double r);
};


#endif /* maillage_hpp */
