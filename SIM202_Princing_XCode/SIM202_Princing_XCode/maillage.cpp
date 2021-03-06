//
//  maillage.cpp
//  Maillage
//
//  Created by Mayeul CHAVANNE on 01/02/2022.
//

#include "maillage.hpp"
#include <cstdlib>
#include <iostream>

using namespace std;

// utilitaire de messages d'erreur
void stop(const char * msg) {
    cout << "ERREUR :" << msg << endl;
    exit(-1);
}

//---------------------------------------------------------------------------
//     classe POINT  (point 2D)
//---------------------------------------------------------------------------
//transformation affine Ax+t : A=(A11,A12,A21,A22) et t=(t1,t2)
Point& Point::tf_affine(const vector<double> & A, const vector<double> &t)
{ double xx=x,yy=y;
    x=A[0]*xx+A[1]*yy+t[0];
    y=A[2]*xx+A[3]*yy+t[1];
    return *this;
}
Point operator+(const Point & P,const Point &Q)
{Point R(P); return R+=Q;}
Point operator-(const Point & P,const Point &Q)
{Point R(P); return R-=Q;}
Point operator*(const Point & P,double a)
{Point R(P); return R*=a;}
Point operator* (double a,const Point & P)
{Point R(P); return R*=a;}
Point operator/ (const Point & P,double a)
{Point R(P); return R/=a;}
bool operator ==(const Point & P,const Point &Q)
{return (P.x==Q.x) && (P.y==Q.y);}
bool operator !=(const Point & P,const Point &Q)
{return !(P==Q);}
bool operator <(const Point & P,const Point &Q)
{if(P.x<Q.x) return true;
 if(P.x>Q.x) return false;
 if(P.y<Q.y) return true;
 return false;}
ostream & operator <<(ostream & os, const Point &P)
{os<<"("<<P.x<<","<<P.y<<")"; return os;}

//---------------------------------------------------------------------------
//     classe Numeros
//---------------------------------------------------------------------------

Numeros::Numeros(int i1, int i2, int i3){
    this->reserve(3);
    this->push_back(i1);
    this->push_back(i2);
    this->push_back(i3);
}

int Numeros::operator() (const int i) const{
    return this->operator[](i-1);
}

int & Numeros::operator() (const int i) {
    return this->operator[](i-1);
}

ostream& operator <<(ostream& out, const Numeros & N){
    out<<'('<< N(1)<<','<<N(2)<<','<<N(3)<<')'<<endl;
    return out;
}

//---------------------------------------------------------------------------
//     classe Maillage  (2D)
//---------------------------------------------------------------------------

//fonction de base maillant le carr?? unit??
void Maillage::maille_carre_unite(int m, int n)
{
    this->sommets.resize((n+1)*(m+1));
    int q;
    double x,y;
    //remplissage de la liste des sommets
    for(int j = 1; j<=n+1; j++)
    {
        for(int i = 1; i<=m+1; i++)
        {
            q = (j-1)*(m+1) + i - 1; //indice du sommet dans le vecteur this->sommets
            x = (i-1.0)/m;           //abscisse du sommet q
            y = (j-1.0)/n;           //ordonn??e du sommet q
            this->sommets[q] = Point(x,y);
        }
    }
    //remplissage de la liste des triangles
    srand(1000000000);
    int i,j;
    double r = 0.0;
    for(int k=1; k<=2*m*n; k++)
    {
        if((k%2)==1)  //selon la parit?? de l'indice du triangle
        {
            j = k/(2*m) + 1;        //calcul de l'unique couple (i,j) qui verifie
            i = ((k%(2*m)) + 1)/2;  // k = 2(j-1)m + 2i - 1, par division euclidienne
            q = (j-1)*(m+1) + i - 1;    //indice du sommet en bas a gauche de la case
            
            r = rand(); //choix al??atoire du d??coupage de la case en 2 triangles
            
            if (r>RAND_MAX/2)        //selon le choix du d??coupage de la case en 2 triangles
            {
                this->triangles.push_back(Numeros(q, q+1, q+m+2));
            }
            else
            {
                this->triangles.push_back(Numeros(q, q+1, q+m+1));
            }
        }
        else
        {
            j = (k-1)/(2*m) + 1;             //calcul de l'unique couple (i,j) qui verifie
            i = (((k-1)%(2*m)) + 1)/2;       // k = 2(j-1)m + 2i, par division euclidienne
            q = (j-1)*(m+1) + i - 1;             //indice du sommet en bas a gauche de la case
            if (r>RAND_MAX/2)           //selon le choix du d??coupage de la case en 2 triangles
            {
                this->triangles.push_back(Numeros(q, q+m+2, q+m+1));
            }
            else
            {
                this->triangles.push_back(Numeros(q+1, q+m+2, q+m+1));
            }
        }
    }
}

//fonction fixant la composante valeur du maillage
void Maillage::valeur(const vector<double> & Q)
{
    unsigned long n = Q.size();
    this->valeurs.resize(n);
    for(int i=0; i<n ; i++){
        this->valeurs[i] = Q[i];
    }
}

//constructeur de maillage carr?? appelant maille_carre_unite et valeurs
Maillage::Maillage(int n, int m, const vector<double> & Q) {
    if(Q.size() != (n+1)*(m+1)) {
        stop(" nombre de sommets et taille du maillage incompatibles");
    }
    this->maille_carre_unite(n, m);
    this->valeur(Q);
}



//fonction effectuant la transformation affine Ax + t
Maillage& Maillage::tf_affine(const vector<double> & A, const vector<double> & t)
{
    vector<Point>::iterator its=this->sommets.begin();
    for(; its!=this->sommets.end(); its++)
    {
        its->tf_affine(A, t);
    }
    return *this;
}

//fonction permettant de mailler le rectangle [a,b]*[c,d] avec m*n cases
void Maillage::maille_rectangle(double a, double b, double c, double d, int n, int m)
{
    this->maille_carre_unite(n, m);
    vector<double> A(4, 0.);
    A[0] = b-a;
    A[3] = d-c;
    vector<double> t(2, 0.);
    t[0] = a;
    t[1] = c;
    this->tf_affine(A, t);
 }

//constructeur de maillage rectangulaire [a,b]*[c,d] avec m*n cases appelant maille_rectangle
Maillage::Maillage( double a, double b, double c, double d, int n, int m, const vector<double> & Q) {
    if(Q.size() != (n+1)*(m+1)) {
        stop(" nombre de sommets et taille du maillage incompatibles");
    }
    this->maille_rectangle(a, b, c, d, n, m);
    this->valeur(Q);
}

//calcul du maillage ulterieur ?? partir d'un maillage existant
void Maillage::resolution(vector<double> & khi, double r) {
    
}

//fonction affichant la liste des sommets et des num??ros des triangles
void Maillage::affiche() const
{
    cout << "Liste des sommets-valeur (" << this->sommets.size() << " points)" << endl;
    vector<Point>::const_iterator its=this->sommets.begin();
    int i = 0;
    while(its!=this->sommets.end())
    {
        cout << "sommet " << i << " : " << *its << " ; valeur = " << this->valeurs[i] << endl;
        i++ ;
        its++ ;
    }
    cout << "Liste des triangles (" << this->triangles.size() << " triangles)" << endl;
    list<Numeros>::const_iterator itt=this->triangles.begin();
    i = 1;
    while(itt!=this->triangles.end())
    {
        cout << "triangle " << i << " : " << *itt << endl;
        i++ ;
        itt++ ;
    }
}

//export du maillage dans un fichier
void Maillage::saveToFile(const char *fn ) const
{
    
    string filename(fn);
    fstream os;
    
//    cout << " nom du fichier " << filename << endl;

    os.open(filename, std::ios_base::out);
    if (!os.is_open()) {
        stop(" echec de l'ouverture du fichier ");
    }
//    cout << "ouverture du fichier " << os.is_open() << endl;
    
    list<Numeros>::const_iterator itt=this->triangles.begin();
    for(; itt!=this->triangles.end() ; itt++)
    {
//        cout << "itt =" << endl;
        for(int i=0; i<3; i++)
        {
            const Point& p=this->sommets[(*itt)[i]];
            os << p.x  << " , " << p.y << endl;
        }
    }
    
    os.close();
    
}
