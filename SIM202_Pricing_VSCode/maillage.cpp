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

//---------------------------------------------------------------------------
//     quelques fonctions utilitaires
//---------------------------------------------------------------------------

// utilitaire de messages d'erreur
void stop(const char * msg) {
    cout << "ERREUR :" << msg << endl;
    exit(-1);
}

//transforme les element d'un vecteur de vecteur en vecteur
vector<double> vecteur(vector< vector<double> > vec, int i) {
    vector<double> vect;
    vector< vector<double> > vecteur_vecteur=vec;
    int n=vecteur_vecteur[i].size();
    for (int k=0; k<n; ++k) {
        vect[k]=vecteur_vecteur[i][k];
    }
return vect;
}

//===================================================================================================
//                                     Classe Matrice
//===================================================================================================


//opérateurs internes et fonctions membres

double Matrice::aff(int i,int j) const {return (*this)(i,j);}

//opérateurs externes
ostream& operator<<(ostream& os,const Matrice& M){
    int d=M.dim;
    os<<"\n";
    for(int i=1;i<=d;i++){
        os<<"|";
        for(int j=1;j<d;j++){os<<M(i,j)<<" ";}
        os<<M(i,d)<<"|"<<"\n";
    }
    cout<<endl;
    return os;
}



ostream& operator<<(ostream& os, const vector<double>& v) 
{
    os << "[";
    for (vector<double>::const_iterator ii = v.begin(); ii != v.end(); ++ii)
    {
        os << " " << *ii;
    }
    os << " ]";
    return os;
}
ostream& operator<<(ostream& os, const vector<int>& v) 
{
    os << "[";
    for (vector<int>::const_iterator ii = v.begin(); ii != v.end(); ++ii)
    {
        os << " " << *ii;
    }
    os << " ]";
    return os;
}




//===================================================================================================
//                                     Classe Matrice_S
//===================================================================================================



//opérateurs internes

double Matrice_S::operator()(int i, int j) const
{
        if(i>this->dim || j>this->dim){
            cout<<"Erreur, hors dimension, la dimension est"<<this->dim;
            exit(EXIT_FAILURE);
        }
        if(j>i){return (*this)(j,i);}
        if(i==j){return this->Mat[this->Stack[i-1]];}
        int Pm1=i-(this->Stack[i-1] - this->Stack[i-2]) + 1; //on détermine le profil de la ligne
        if(j<Pm1){return 0.;}
        return this->Mat[this->Stack[i-1]-(i-j)];
}

double& Matrice_S::operator()(int i, int j)
{
 if(i>this->dim || j>this->dim){
            cout<<"Erreur, hors dimension, la dimension est"<<this->dim;
            exit(EXIT_FAILURE);
        }
        if(j>i){return (*this)(j,i);}
        if(i==j){return this->Mat[this->Stack[i-1]];}
        int Pm2=i-(this->Stack[i-1] - this->Stack[i-2]) + 1; //on détermine le profil de la ligne
        if(j<Pm2){cout<<"Erreur, hors profil, j doit être sup à "<<Pm2<<" pour la ligne "<<i<<endl; exit(-1);}
        return this->Mat[this->Stack[i-1]-(i-j)];
}

//**opérateurs algébriques sur la même classe
Matrice_S&  Matrice_S::operator+=(const Matrice_S& M){
    if(this->dim!=M.d()){cout<<"les deux matrcies n'ont pas même dimensions"<<endl; exit(-1);}
    vector<int> S(M.d(),0);
    for(int i=1;i<M.d();i++){
        S[i]=S[i-1]+max(this->Stack[i]-this->Stack[i-1],M.Stack[i]-M.Stack[i-1]);    
    }
    Matrice_S U(S);
    U(1,1)=M(1,1)+(*this)(1,1);
    for(int i=1;i<U.d();i++){
        for(int j=i-(U.Stack[i]-U.Stack[i-1])+1;j<=i;j++){
            U(i+1,j+1)=M(i+1,j+1)+this->aff(i+1,j+1); // (*this)(i,j) fait expressement appel à la fonction non constante d'écriture et on ne peut pas utiliser const_cast   
        }
    }
    this->Stack=U.Stack;
    this->Mat=U.Mat;
    return *this;
}
Matrice_S& Matrice_S::operator-=(const Matrice_S& M){
    if(this->dim!=M.d()){cout<<"les deux matrcies n'ont pas même dimensions"<<endl; exit(-1);}  //On crée un profil adapté à la somme des deux matrices
    vector<int> S(M.d(),0);
    for(int i=1;i<M.d();i++){
        S[i]=S[i-1]+max(this->Stack[i]-this->Stack[i-1],M.Stack[i]-M.Stack[i-1]);    
    }
    Matrice_S U(S);                                                                 //On crée une Matrice temporaire à laquelle on attribue le profil crée
    U(1,1)=M(1,1)-(*this)(1,1);                                                     //On ne parcours que les indices intéressants par soucis d'efficacité
    for(int i=1;i<U.d();i++){
        for(int j=i-(U.Stack[i]-U.Stack[i-1])+1;j<=i;j++){
            U(i+1,j+1)=this->aff(i+1,j+1)-M(i+1,j+1); // (*this)(i,j) fait expressement appel à la fonction non constante d'écriture et on ne peut pas utiliser const_cast   
        }
    }
    this->Stack=U.Stack;                              //On attribue à l'élément courant les caractéristiques de U.
    this->Mat=U.Mat;
    return *this;
}
Matrice_S& Matrice_S::operator*=(const double& a)
{   for(int i=0;i<this->Mat.size();i++){this->Mat[i]*=a;}
    return *this;
}

Matrice_S& Matrice_S::operator/=(const double& a)
{
    if(a==0){cout<<"on ne divise pas par 0"<<endl; exit(-1);}
    return (*this)*=(1/a);
}

Matrice_S& Matrice_S::operator=(const Matrice_S& M)
{
    if(this->d()!=M.d()){cout<<"attention les matrices n'étaient pas de même taille, l'opération a quand même été réalisée"<<endl;}
    this->dim=M.d();
    this->Stack=M.Stack;
    this->Mat=M.Mat;
    return *this;
}



//fonctions externes
//**opérateurs algébriques
Matrice_S operator+(const Matrice_S& M,const Matrice_S& H){ Matrice_S A(M); return A+=H;}
Matrice_S operator-(const Matrice_S& M,const Matrice_S& H){Matrice_S A(M); return A-=H;}
Matrice_S operator*(const Matrice_S& M,const double& a){Matrice_S A(M); return A*=a;}
Matrice_S operator/(const Matrice_S& M,const double& a){Matrice_S A(M); return A/=a;}

vector<double> operator*(const Matrice_S& M,const vector<double>& X)
{   
    int d=X.size();
    if(M.d()!=d){cout<<"erreur dimensions pour le produit matrice vecteur "<<M.d()<<" ; "<<d<<endl; exit(-1);}
    vector<double> u(d,0.);
    for(int i=0;i<d;i++)
    {
        for(int j=0;j<d;j++){u[i]+=M(i+1,j+1)*X[j];}
    }
    return u;
}

vector<double> operator*(const Matrice_S& M,const vector<int>& X)
{   
    int d=X.size();
    if(M.d()!=d){cout<<"erreur dimensions pour le produit matrice vecteur "<<M.d()<<" ; "<<d<<endl; exit(-1);}
    vector<double> u(d,0.);
    for(int i=0;i<d;i++)
    {
        for(int j=0;j<d;j++){u[i]+=M(i+1,j+1)*X[j];}
    }
    return u;
}



//===================================================================================================
//                                     Classe Matrice_PS
//===================================================================================================


//opérateurs internes
double Matrice_PS::operator()(int i, int j) const {
    if(i>this->dim || j>this->dim){
        cout<<"Erreur, hors dimension, la dimension est"<<this->dim;
        exit(EXIT_FAILURE);
    }
    if(i==j){return this->Mat[this->Stack[i-1]];}
    if(j>i){
        int Pm1=j-(this->Stack[j-1] - this->Stack[j-2]) + 1;
        if(i<Pm1){return 0.;}
        return this->Mat[this->Stack[j-1]-(j-i) + (this->Stack[this->dim-1]+1)-(j-1)];
    }
    int Pm1=i-(this->Stack[i-1] - this->Stack[i-2]) + 1; //on détermine le profil de la ligne
    if(j<Pm1){return 0.;}
    return this->Mat[this->Stack[i-1]-(i-j)];
}

double& Matrice_PS::operator()(int i, int j) {
    if(i>this->dim || j>this->dim){
        cout<<"Erreur, hors dimension, la dimension est"<<this->dim;
        exit(EXIT_FAILURE);
    }
    if(i==j){return this->Mat[this->Stack[i-1]];}
    if(j>i){
        int Pm2=j-(this->Stack[j-1] - this->Stack[j-2]) + 1;
        if(i<Pm2){cout<<"Erreur, hors profil, i doit être sup à "<<Pm2; exit(-1);}
        return this->Mat[this->Stack[j-1]-(j-i) + (this->Stack[this->dim-1]+1)-(j-1)];
    }
    int Pm2=i-(this->Stack[i-1] - this->Stack[i-2]) + 1; //on détermine le profil de la ligne
    if(j<Pm2){cout<<"Erreur, hors profil, j doit être sup à "<<Pm2; exit(-1);}
    return this->Mat[this->Stack[i-1]-(i-j)];
}

//**opérateurs algébriques sur la même classe

Matrice_PS& Matrice_PS::operator+=(const Matrice_PS& M)
{
    if(this->dim!=M.d()){cout<<"les deux matrcies n'ont pas même dimensions"<<endl; exit(-1);}
    vector<int> S(M.d(),0);
    for(int i=1;i<M.d();i++){
        S[i]=S[i-1]+max(this->Stack[i]-this->Stack[i-1],M.Stack[i]-M.Stack[i-1]);    
    }
    Matrice_PS U(S);
    U(1,1)=M(1,1)+(*this)(1,1);
    for(int i=1;i<U.d();i++){
        U(i+1,i+1)=this->aff(i+1,i+1)+M(i+1,i+1);
        for(int j=i-(U.Stack[i]-U.Stack[i-1])+1;j<i;j++){
            U(i+1,j+1)=this->aff(i+1,j+1)+M(i+1,j+1); // (*this)(i,j) fait expressement appel à la fonction non constante d'écriture et on ne peut pas utiliser const_cast
            U(j+1,i+1)=this->aff(j+1,i+1)+M(j+1,i+1);
        }
    }
    this->Stack=U.Stack;
    this->Mat=U.Mat;
    return *this;

}

Matrice_PS& Matrice_PS::operator-=(const Matrice_PS& M)
{
    if(this->dim!=M.d()){cout<<"les deux matrcies n'ont pas même dimensions"<<endl; exit(-1);}
    vector<int> S(M.d(),0);
    for(int i=1;i<M.d();i++){
        S[i]=S[i-1]+max(this->Stack[i]-this->Stack[i-1],M.Stack[i]-M.Stack[i-1]);    
    }
    Matrice_PS U(S);
    U(1,1)=M(1,1)-(*this)(1,1);
    for(int i=1;i<U.d();i++){
        U(i+1,i+1)=this->aff(i+1,i+1)-M(i+1,i+1);
        for(int j=i-(U.Stack[i]-U.Stack[i-1])+1;j<i;j++){
            U(i+1,j+1)=this->aff(i+1,j+1)-M(i+1,j+1); // (*this)(i,j) fait expressement appel à la fonction non constante d'écriture et on ne peut pas utiliser const_cast
            U(j+1,i+1)=this->aff(j+1,i+1)-M(j+1,i+1);
        }
    }
    this->Stack=U.Stack;
    this->Mat=U.Mat;
    return *this;

}
Matrice_PS& Matrice_PS::operator*=(const double& a)
{   for(int i=0;i<this->Mat.size();i++){this->Mat[i]*=a;}
    return *this;
}

Matrice_PS& Matrice_PS::operator/=(const double& a)
{
    if(a==0){cout<<"on ne divise pas par 0"<<endl; exit(-1);}
    return (*this)*=(1/a);
}


//**opérateurs algébriques sur Matrice_S

Matrice_PS& Matrice_PS::operator+=(const Matrice_S& M)
{
    if(this->dim!=M.d()){cout<<"les deux matrcies n'ont pas même dimensions"<<endl; exit(-1);}
    vector<int> S(M.d(),0);
    for(int i=1;i<M.d();i++){
        S[i]=S[i-1]+max(this->Stack[i]-this->Stack[i-1],M.Stack[i]-M.Stack[i-1]);    
    }
    Matrice_PS U(S);
    U(1,1)=M(1,1)+(*this)(1,1);
    for(int i=1;i<U.d();i++){
        U(i+1,i+1)=this->aff(i+1,i+1)+M(i+1,i+1);
        for(int j=i-(U.Stack[i]-U.Stack[i-1])+1;j<i;j++){
            U(i+1,j+1)=this->aff(i+1,j+1)+M(i+1,j+1); // (*this)(i,j) fait expressement appel à la fonction non constante d'écriture et on ne peut pas utiliser const_cast
            U(j+1,i+1)=this->aff(j+1,i+1)+M(j+1,i+1);
        }
    }
    this->Stack=U.Stack;
    this->Mat=U.Mat;
    return *this;
}

Matrice_PS& Matrice_PS::operator-=(const Matrice_S& M)
{
    if(this->dim!=M.d()){cout<<"les deux matrcies n'ont pas même dimensions"<<endl; exit(-1);}
    vector<int> S(M.d(),0);
    for(int i=1;i<M.d();i++){
        S[i]=S[i-1]+max(this->Stack[i]-this->Stack[i-1],M.Stack[i]-M.Stack[i-1]);    
    }
    Matrice_PS U(S);
    U(1,1)=M(1,1)-(*this)(1,1);
    for(int i=1;i<U.d();i++){
        U(i+1,i+1)=this->aff(i+1,i+1)-M(i+1,i+1);
        for(int j=i-(U.Stack[i]-U.Stack[i-1])+1;j<i;j++){
            U(i+1,j+1)=this->aff(i+1,j+1)-M(i+1,j+1); // (*this)(i,j) fait expressement appel à la fonction non constante d'écriture et on ne peut pas utiliser const_cast
            U(j+1,i+1)=this->aff(j+1,i+1)-M(j+1,i+1);
        }
    }
    this->Stack=U.Stack;
    this->Mat=U.Mat;
    return *this;

}
Matrice_PS& Matrice_PS::operator=(const Matrice_PS& M)
{
    if(this->d()!=M.d()){cout<<"attention les matrices n'étaient pas de même taille, l'opération a quand même été réalisée"<<endl;}
    this->dim=M.d();
    this->Stack=M.Stack;
    this->Mat=M.Mat;
    return *this;
}

Matrice_PS& Matrice_PS::operator=(const Matrice_S& M)
{
    if(this->d()!=M.d()){cout<<"attention les matrices n'étaient pas de même taille, l'opération a quand même été réalisée"<<endl;}
    this->dim=M.d();
    this->Stack=M.Stack;
    this->Mat=M.Mat;
        vector<double>::iterator it1=this->Mat.end();
        vector<double> u(this->Mat);
        vector<double>::iterator it2=u.begin();
        vector<double>::iterator it3=u.end();
        for(int i=0;i<this->d();i++){u.erase(it2+this->Stack[i]);}
        this->Mat.insert(it1,it2,it3);
    return *this;
}



//Opérateurs externes
Matrice_PS operator+(const Matrice_PS& M,const Matrice_PS& H){Matrice_PS A(M); return A+=H;}
Matrice_PS operator+(const Matrice_PS& M,const Matrice_S& H){Matrice_PS A(M); return A+=H;}
Matrice_PS operator+(const Matrice_S& M,const Matrice_PS& H){Matrice_PS A(H); return A+=M;}

Matrice_PS operator-(const Matrice_PS& M,const Matrice_PS& H){Matrice_PS A(M); return A-=H;}
Matrice_PS operator-(const Matrice_PS& M,const Matrice_S& H){Matrice_PS A(M); return A-=H;}
Matrice_PS operator-(const Matrice_S& M,const Matrice_PS& H){Matrice_PS A(H); return A-=M;}

Matrice_PS operator*(const Matrice_PS& M,const double& a){Matrice_PS A(M); return A*=a;}
Matrice_PS operator/(const Matrice_PS& M,const double& a){Matrice_PS A(M); return A/=a;}

vector<double> operator*(Matrice_PS& M,vector<double>& X)
{   
    int d=X.size();
    if(M.d()!=d){cout<<"erreur dimensions pour le produit matrice vecteur "<<M.d()<<" ; "<<d<<endl; exit(-1);}
    vector<double> u(d,0.);
    for(int i=0;i<d;i++)
    {
        for(int j=0;j<d;j++){u[i]+=M(i+1,j+1)*X[j];}
    }
    return u;
}
vector<double> operator*(Matrice_PS& M,vector<int>& X)
{   
    int d=X.size();
    if(M.d()!=d){cout<<"erreur dimensions pour le produit matrice vecteur "<<M.d()<<" ; "<<d<<endl; exit(-1);}
    vector<double> u(d,0.);
    for(int i=0;i<d;i++)
    {
        for(int j=0;j<d;j++){u[i]+=M(i+1,j+1)*X[j];}
    }
    return u;
}


//factorisation LU
vector<Matrice_PS> Matrice_PS::LU()
{
    double d=this->d();
    Matrice_PS U((*this).Stack);
    Matrice_PS L((*this).Stack);
    double p=0;
    for(int i=1;i<=d;i++)
    {
        p=(*this)(i,i);
        L(i,i)=1;
        for(int j=i+1;j<=d;j++)
        {
            L(j,i)=(*this)(j,i)/p;
        }
        U(i,i)=p;
        for(int k=i+1;k<=d;k++){U(i,k)=(*this)(i,k);}
        for(int k=i+1;k<=d;k++)
        {
            for(int j=i+1;j<=d;j++){(*this)(k,j)=(*this)(k,j)-L(k,i)*U(i,j);}

        }
        
    }
    vector<Matrice_PS> T(2);
    T[0]=L;
    T[1]=U;
    return T;
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

//operateur d'affectation
void Numeros::operator = (const Numeros & N) {
    this->operator[](0) = N[0];
    this->operator[](1) = N[1];
    this->operator[](2) = N[2];
}

//---------------------------------------------------------------------------
//     classe Maillage  (2D)
//---------------------------------------------------------------------------

//constructeur par copie d'un maillage
Maillage::Maillage(Maillage & M) {
    unsigned long n = this->sommets.size();
    if(n != this->valeurs.size()) {
        stop(" nombre de sommets et taille du maillage incompatibles");
    }

    this->sommets.resize(n);
    this->valeurs.resize(n);
    this->triangles.resize(M.triangles.size());
    for(int i=0; i<n ; i++) {
        this->sommets[i] = M.sommets[i];
        this->valeurs[i] = M.valeurs[i];
    }
    list<Numeros>::const_iterator itt=M.triangles.begin();
    list<Numeros>::const_iterator itt=M.triangles.begin();

    while(itt!=this->triangles.end())
    {
        this->triangles.list<Numeros>::operator[](i) = *itt;
        itt++ ;
        i++
    }
}

//fonction de base maillant le carré unité avec m*n cases (n lignes, m colonnes)
void Maillage::maille_carre_unite(int m, int n) {
    this->sommets.resize((n+1)*(m+1));
    int q;
    double x,y;
    //remplissage de la liste des sommets ; j indice des lignes, i indices des colonnes
    for(int j = 1; j<=n+1; j++)
    {
        for(int i = 1; i<=m+1; i++)
        {
            q = (j-1)*(m+1) + i - 1; //indice du sommet dans le vecteur this->sommets
            x = (i-1.0)/m;           //abscisse du sommet q
            y = (j-1.0)/n;           //ordonnée du sommet q
            this->sommets[q] = Point(x,y);
        }
    }
    //remplissage de la liste des triangles
    srand(1000000000);
    int i,j;
    double r = 0.0;
    for(int k=1; k<=2*m*n; k++)
    {
        if((k%2)==1)  //selon la parité de l'indice du triangle
        {
            j = k/(2*m) + 1;        //calcul de l'unique couple (i,j) qui verifie
            i = ((k%(2*m)) + 1)/2;  // k = 2(j-1)m + 2i - 1, par division euclidienne
            q = (j-1)*(m+1) + i - 1;    //indice du sommet en bas a gauche de la case
            
            r = rand(); //choix aléatoire du découpage de la case en 2 triangles
            
            if (r>RAND_MAX/2)        //selon le choix du découpage de la case en 2 triangles
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
            if (r>RAND_MAX/2)           //selon le choix du découpage de la case en 2 triangles
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
void Maillage::valeur(const vector<double> & Q) {
    unsigned long n = Q.size();
    if(n != this->sommets.size()) {
        stop(" nombre de sommets et taille du maillage incompatibles");
    }
    this->valeurs.resize(n);
    for(int i=0; i<n ; i++){
        this->valeurs[i] = Q[i];
    }
}

//constructeur de maillage carré appelant maille_carre_unite avec m*n cases (n lignes, m colonnes), et fixant les valeurs initiales a Q
Maillage::Maillage(int n, int m, const vector<double> & Q) {
    if(Q.size() != (n+1)*(m+1)) {
        stop(" nombre de sommets et taille du maillage incompatibles");
    }
    this->maille_carre_unite(n, m);
    this->valeur(Q);
}

//fonction effectuant la transformation affine Ax + t
Maillage& Maillage::tf_affine(const vector<double> & A, const vector<double> & t) {
    vector<Point>::iterator its=this->sommets.begin();
    for(; its!=this->sommets.end(); its++)
    {
        its->tf_affine(A, t);
    }
    return *this;
}

//fonction permettant de mailler le rectangle [a,b]*[c,d] avec m*n cases
void Maillage::maille_rectangle(double a, double b, double c, double d, int n, int m) {
    this->maille_carre_unite(n, m);
    vector<double> A(4, 0.);
    A[0] = b-a;
    A[3] = d-c;
    vector<double> t(2, 0.);
    t[0] = a;
    t[1] = c;
    this->tf_affine(A, t);
}

//constructeur de maillage rectangulaire [a,b]*[c,d] avec m*n cases appelant maille_rectangle, et fixant les valeurs initiales a Q
Maillage::Maillage( double a, double b, double c, double d, int n, int m, const vector<double> & Q) {
    if(Q.size() != (n+1)*(m+1)) {
        stop(" nombre de sommets et taille du maillage incompatibles");
    }
    this->maille_rectangle(a, b, c, d, n, m);
    this->valeur(Q);
}

//ASSEMBLAGE DE LA MATRICE B
Matrice_PS& maillage::mat_B() { //A PARTIR DE LA LISTE DES TRIANGLES CONTENANT LES SOMMETS
    vector<int> vect_BB=this->indice_profil();
    Matrice_PS& BB(profil(vect_BB)); //matrice carree profilee

    vector<int> vect_B_el(3,1);
    Matrice_PS& Bel(profil(vect_B_el)); //matrice elementaire B
    for (int l=1; l<=this->triangles.size(); ++l) { //boucle sur le nombre de triangles
        Numeros Num(this->triangles[l]) //l-eme ligne de la matrice stockant les triangles
        int i1=Num[1]; int i2=Num[2]; int i3=Num[3]; //indices des sommets du l-eme triangle
        Bel=mat_elem_B(this->sommets[i1],this->sommets[i2],this->sommets[i3]); //matrice elementaire par acces aux coordonnees du (resp.) 1,2,3-eme sommet du l-eme triangle
        for (int i=1;i<=3;++i) {
            int I=this->Num[i];
            for (int j=1;j<=3;++j) {
                J=this->Num[j]; 
                BB(I,J)+=Bel(i,j);
            }
        }
    }
return BB;
}

//ASSEMBLAGE DE LA MATRICE MASSE
Matrice_S& Maillage::mat_M() { //A PARTIR DE LA LISTE DES TRIANGLES CONTENANT LES SOMMETS
    vector<int> vect_MM=this->indice_profil();
    Matrice_S& MM(profil(vect_MM)); //matrice carree taille nombre de sommets
    vector<int> vect_M_el(3,1);
    Matrice_S& Mel(profil(vect_M_el)); //matrice elementaire de masse
    for (int l=1; l<=this->triangles.size(); ++l) { //boucle sur le nombre de triangles
        Numeros Num(this->triangles[l]) //l-eme ligne de la matrice stockant les triangles
        int i1=Num[1]; int i2=Num[2]; int i3=Num[3]; //indices des sommets du l-eme triangle
        Mel=mat_elem_M(this->sommets[i1],this->sommets[i2],this->sommets[i3]); //matrice elementaire par acces aux coordonnees du (resp.) 1,2,3-eme sommet du l-eme triangle
        for (int i=1;i<=3;++i) {
            int I=this->Num[i];
            for (int j=1;j<=3;++j) {
                J=this->Num[j]; 
                MM(I,J)+=Mel(i,j);
            }
        }
    }
return MM;
}

//ASSEMBLAGE DE LA MATRICE K
Matrice Maillage::mat_K() { 
    vector<int> vect_KK=this->indice_profil();
    Matrice_S& KK(profil(vect_KK)); //matrice carree profilee
    vector<int> vect_K_el(3,1)
    Matrice_S& Kel(profil(vect_K_el)); //matrice elementaire de masse
    for (int l=1; l<=this->triangles.size(); ++l) { //boucle sur le nombre de triangles
        Numeros Num(this->triangles[l]) //l-eme ligne de la matrice stockant les triangles
        int i1=Num[1]; int i2=Num[2]; int i3=Num[3]; //indices des sommets du l-eme triangle
        Kel=mat_elem_K(this->sommets[i1],this->sommets[i2],this->sommets[i3]); //matrice elementaire par acces aux coordonnees du (resp.) 1,2,3-eme sommet du l-eme triangle
        for (int i=1;i<=3;++i) {
            int I=this->Num[i]; //index (entier) du sommet i du triangle l
            for (int j=1;j<=3;++j) {
                J=this->Num[j]; //index (entier) du sommet i du triangle l
                KK(I,J)+=Kel(i,j);
            }
        }
    }
return KK;
}

//MATRICE D
Matrice_PS& Maillage::mat_D() {
    vector<int> vect_DD=this->indice_profil();
    Matrice_PS& DD(profil(vect_DD));
    Matrice M(mat_M());
    Matrice B(mat_B());
    Matrice K(mat_M());
    DD=r*M+B+K;
return DD;
}

//CALCUL MATRICE SOURCE A
Matrice_S& source_A(Point *P) {
    vector<int> v(2,1);
    Matrice_S& A(profil(v));
    A(1,1)=khi[0]*P.x*P.x/2; A(2,1)=khi[1]*P.x*P.y/2; //calcul de chaque coefficient de la matrice A, khi1=khi2, A symétrique
    A(1,2)=khi[2]*P.x*P.y/2; A(2,2)=khi[3]*P.y*P.y/2;
return A;
}

//CALCUL DE LA SOURCE V
vector<double> source_V(Point *P) {
    vector<double> V; 
    V.push_back((khi[0]+khi[1]/2-r)*P.x); 
    V.push_back((khi[3]+khi[2]/2-r)*P.y);
return V;
}

//CALCUL DE MATRICE ELEMENTAIRE M
Matrice_S& mat_elem_M (Point * P1,Point * P2,Point * P3){
    double x1 = P1.x; double y1 = P1.y; //coordonnees des points du triangle choisi pour calcul sa matrice de masse (elementaire)
    double x2 = P2.x; double y2 = P2.y;
    double x3 = P3.x; double y3 = P3.y;

    double D = ((x2-x1)*(y3-y1) - (y2-y1)*(x3-x1)); //D est, au signe pres, deux fois l'aire du triangle


    vector<int> vect_M_elem(3,1);
    Matrice_S& M_elem(profil(vect_M_elem)); //creation matrice 3x3 elementaire, a voir selon syntaxe

    double Abs_Det_Bl = abs(D); //determinant jacobienne pour changement dans le triangle de base

    vector<int> vect_m_des_w(3,1);
    Matrice_S& M_des_w(profil(vect_m_des_w));
    M_des_w(1,1) = 1/3 ; M_des_w(1,2) = 1/6; M_des_w(1,3) = 1/6; //calcul des elements via les fonctions barycentriques(?)
    M_des_w(2,1) = 1/6; M_des_w(2,2) = 2/3; M_des_w(2,3) = 1/6;
    M_des_w(3,1) = 1/6; M_des_w(3,2) = 1/6; M_des_w(3,3) = 2/3;

    for (int i=1;i<=3;++i){
        for (int j=1;j<=3;++j){
            Terme_1 = (1/6)*M_des_w(i,1)*M_des_w(j,1);
            Terme_2 = (1/6)*M_des_w(i,2)*M_des_w(j,2);
            Terme_3 = (1/6)*M_des_w(i,3)*M_des_w(j,3);
		    M_elem(i,j) = (Terme_1 + Terme_2 + Terme_3)*Abs_Det_Bl; //calcul des coefficients de la matrice elementaire
        };
    };
return M_elem;
};

//CALCUL MATRICE ELEMENTAIRE K
Matrice_S& mat_elem_K (Point * P1,Point * P2,Point * P3) {
    double x1 = P1.x; double y1 = P1.y; //coordonnees des points du triangle choisi pour calcul sa matrice de masse (elementaire)
    double x2 = P2.x; double y2 = P2.y;
    double x3 = P3.x; double y3 = P3.y;

    double D = ((x2-x1)*(y3-y1) - (y2-y1)*(x3-x1)); //au signe pres, deux fois la taille du triangle choisi
    if (abs(D) <= eps) { //fixer un eps petit pour verifier que le triangle n'est pas plat
        exit('l aire d un triangle est nulle'); 
    }; 

    vector<int> vect_K_el(3,1);
    Matrice_S& Kel(profil(vect_K_el)); //matrice elementaire taille 3x3

    vector<int> vect_Bl(2,1);
    Matrice_PS& Bl(profil(vect_Bl)); //Matrice pour ponderer les changements de base des coordonnees des gradients
    Bl(1,1)= x2-x1; Bl(1,2)= x3-x1;
    Bl(2,1)= y2-y1; Bl(2,2)= y3-y1;

    double Abs_Det_Bl = abs(D);

    vector<int> vect_inv_Bl(2,1);
    Matrice_PS& inv_Bl(profil(vect_inv_Bl)); // transposee de BL^{-1}
    inv_Bl(1,1)= Bl(2,2)/D; inv_Bl(1,2)= -Bl(2,1)/D;
    inv_Bl(2,1)= -Bl(1,2)/D; inv_Bl(1,1)= Bl(1,1)/D;

    vector<int> vect_id(2,1);
    Matrice_S& Id(profil(vect_id));
    Id(1,1)= 1; Id(1,2)= 0;
    Id(2,1)= 0; Id(2,2)= 1;
    // utiliser Reftri qui est une liste qui sert de reference aux triangles
    vector<double> S;
    S.push_back(x1); S.push_back(y1);

    vector<double> S_chap_1;
    S_chap_1.push_back(1/6); S_chap_1.push_back(1/6);
    vector<double> S_chap_2;
    S_chap_2.push_back(2/3); S_chap_2.push_back(1/6);    
    vector<double> S_chap_3;
    S_chap_3.push_back(1/6); S_chap_3.push_back(2/3);

    vector<double> Fl_1;
    Fl_1=Bl*S_chap_1 + S;
    vector<double> Fl_2;
    Fl_2=Bl*S_chap_2 + S;
    vector<double> Fl_3;
    Fl_3=Bl*S_chap_3 + S;

    Point Fl__1(Fl_1[1],Fl_1[2]));
    Point Fl__2(Fl_2[1],Fl_2[2]);
    Point Fl__3(Fl_3[1],Fl_3[2]);
    
    vector<double> Mat_des_gradients_w_1;
    Mat_des_gradients_w_1.push_back(-1); Mat_des_gradients_w_1.push_back(-1);
    vector<double> Mat_des_gradients_w_2;
    Mat_des_gradients_w_2.push_back(1); Mat_des_gradients_w_2.push_back(0);
    vector<double> Mat_des_gradients_w_3;
    Mat_des_gradients_w_3.push_back(0); Mat_des_gradients_w_3.push_back(1);

    vector<vector<double> > Mat_des_gradients_w;
    Mat_des_gradients_w.push_back(Mat_des_gradients_w_1);
    Mat_des_gradients_w.push_back(Mat_des_gradients_w_2);
    Mat_des_gradients_w.push_back(Mat_des_gradients_w_3);

    for (int i=0;i<=2;++i) {
        for (int j=0;j<=2;++j) {
            double Terme_1 = (1/6)*Inv_Bl*Inv_Bl*(source_A(Fl__1)*vecteur(Mat_des_gradients_w[i]))*vecteur(Mat_des_gradients_w[j])*Abs_Det_Bl; // fonction qui permet d'isoler une colonne ou une ligne de matrice (:,i)
            double Terme_2 = (1/6)*Inv_Bl*Inv_Bl*(source_A(Fl__2)*Mat_des_gradients_w[i])*Mat_des_gradients_w[j]*Abs_Det_Bl;
            double Terme_3 = (1/6)*Inv_Bl*Inv_Bl*(source_A(Fl__3)*Mat_des_gradients_w[i])*Mat_des_gradients_w[j]*Abs_Det_Bl;
            Kel(i,j) = Terme_1 + Terme_2 + Terme_3;
        }
    }
return Kel;
}

//CALCUL MATRICE ELEMENTAIRE B
Matrice_PS& mat_elem_B(Point * P1,Point * P2,Point * P3)) {
    double x1 = P1.x; double y1 = P1.y; //coordonnees des points du triangle choisi pour calcul sa matrice de masse (elementaire)
    double x2 = P2.x; double y2 = P2.y;
    double x3 = P3.x; double y3 = P3.y;

    double D = ((x2-x1)*(y3-y1) - (y2-y1)*(x3-x1)); //au signe pres, deux fois la taille du triangle choisi

    vector<vector<double>> Bl;
    vector<double> Bl_1;
    Bl_1.push_back(y2-y3); Bl_1.push_back(x3-x2);
    vector<double> Bl_2;
    Bl_2.push_back(y3-y1); Bl_2.push_back(x1-x3);
    vector<double> Bl_3;
    Bl_3.push_back(y1-y2); Bl_3.push_back(x2-x1);
    Bl.push_back(Bl_1); Bl.push_back(Bl_2); Bl.push_back(Bl_3);

    vector<int> vect_M_des_w(3,1);
    Matrice_PS& M_des_w(profil(vect_M_des_w));
    M_des_w(1,1) = 1-x1-y1 ; M_des_w(1,2) = 1-x2-y2; M_des_w(1,3) = 1-x3-y3; //calcul des elements via les fonctions barycentriques(?)
    M_des_w(2,1) = x1;       M_des_w(2,2) = x2;      M_des_w(2,3) = x3;
    M_des_w(3,1) = y1;       M_des_w(3,2) = y2;      M_des_w(3,3) = y3;

    vector<vector<double> > S_chap;
    vector<double> S_chap_1;
    S_chap_1.push_back(1/6); S_chap_1.push_back(1/6);
    vector<double> S_chap_2;
    S_chap_2.push_back(2/3); S_chap_2.push_back(1/6);    
    vector<double> S_chap_3;
    S_chap_3.push_back(1/6); S_chap_3.push_back(2/3);
    S_chap.push_back(S_chap_1); S_chap.push_back(S_chap_2); S_chap.push_back(S_chap_3); 

    vector<int> vect_B_el(3,1);
    Matrice_PS& Bel(profil(vect_B_el));
    vector<int> vect_temp(3,1);
    Matrice temp(profil(vect_temp));
    for (int i=0; i<=2; ++i) {
        for (int j=0; j<=2; ++j) {
            temp(i+1,j+1)=vecteur(Bl[i])*vecteur(S_chap[j]);
        }
    }

    Bel=(1/6)*temp*M_des_w;

return Bel;
}

//INVERSION D'UN SYSTEME LINEAIRE PAR FACTORISATION LU
vector<double>& Maillage::Gauss(){
    vector<double> P_k = this->valeurs;
    int n=this->valeurs.size();


    Matrice_S& M((*this).mat_M());
    Matrice_PS& D( (*this).mat_D());


    Matrice_PS& E(M+(delta_t)*D); //à définir delta_t
    vector<double>& B;
    B=M*P_k;
    E.LU();

    Matrice_PS* G=E.LU();
    Matrice_PS& L(G[0]);
    Matrice_PS& U(G[1]);

    vector<double>& Y(n);
    Y[0]=B[0];
    for (int j=1; j<n; ++j) {
        for (int k=1; k<=i; ++k) {
            Y[j]=B[j]-L(j+1,k+1)*Y[k];
        }
    }

    vector<double>& X(n);
    X[0]=Y[0];
    for (int i=1; i<=n; ++i) {
        for (int k=i+1; k<n; ++i) {
            X[i]=Y[i]-U(i+1,k+1)*X[k];
        }
    }

    for (int p=0; p<n; ++p) {
        X[i]=X[i]/U(i,i);
    }
return X;
}

//calcul du maillage ulterieur à partir d'un maillage existant
Maillage& Maillage::resolution(vector<double> & khi, double r) {
    Maillage M(*this);
    M.valeurs = this->Gauss();
}

//fonction affichant la liste des sommets et des numéros des triangles
void Maillage::affiche() const{
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
void Maillage::saveToFile(const char *fn ) const {
    
    string filename(fn);
    fstream os;

    os.open(filename, std::ios_base::out);
    if (!os.is_open()) {
        stop(" echec de l'ouverture du fichier ");
    }
    
    list<Numeros>::const_iterator itt=this->triangles.begin();
    for(; itt!=this->triangles.end() ; itt++)
    {
        for(int i=0; i<3; i++)
        {
            const Point& p=this->sommets[(*itt)[i]];
            os << p.x  << " , " << p.y << endl;
        }
    }
    
    os.close();
    
}

//fonction generant le vecteur des indices des premiers coefficients non nuls de chaque 
//lignes du profil des matrices associées au maillage m*n cases (n lignes, m colonnes)
vector<int> indice_profil() {
    int N =this->sommets.size();
    vector<int> P(N, N);
    Numeros tmp(0, 0, 0);

    list<Numeros>::const_iterator itt=this->triangles.begin();
    while(itt!=this->triangles.end())
    {
        tmp = *itt;
        for(int i=0; i<3 ; i++) {
            for(int j=0; j<3; j++) {
                if(j<tmp[i]) { P[i] = j ;};
            }
        }
    }
    return P;
}

//---------------------------------------------------------------------------
//     classe Video  (2D)
//---------------------------------------------------------------------------

//constructeur par valeur pour M instant sur un rectangle [a,b]*[c,d] avec m*n cases, et fixant les valeurs initiales a Q
Video::Video(int M, double a, double b, double c, double d, int n, int m, const vector<double> & Q) {
    this->maille.maille_rectangle(a, b, c, d, n, m);
    this->maille.valeur(Q);
    this->images.resize(M);
    list<vector<double> >::const_iterator itim=this->images.begin();
    while(itim!=this->images.end()) {
        *itim vector<double>::operator= Q ; 
        itim++;
    }
}

//calcul du processus itératif pour un probleme donne
void Video::resolution(vector<double> & khi, double r) {
    list<vector<double> >::iterator itim=this->images.begin(); //itérateur sur les vecteurs image

    Maillage tmp(this->maille);                 //maillage temporaire initialisé à la valeur this->maille
    *itim = tmp.valeurs;
    itim++ ;
    while(itim!=this->images.end()) {
        tmp = tmp.resolution(khi, r);           //calcul du nouveau maillage temporaire
        *itim = tmp.valeurs;                    //stockage du vecteur image dans this->images
        itim++;
    }
}