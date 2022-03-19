#include Matrice.hpp
#include Maillage.hpp


double condition_initiale(Point * P) {
    return //fonction du point P;
};


double second_membre(Point * P) {
    return //fonction du point P;
};

Matrice khi(2)
//remplir la matrice de covariance

Matrice coeff_A(Point *P) {
    Matrice A(2);
    A(1,1)=khi(1,1)*P(1)*P(1)/2; // P(1)=x1, voir comment ca s'ecrit
    A(2,1)=khi(2,1)*P(1)*P(2)/2;
    A(1,2)=khi(1,2)*P(1)*P(2)/2;
    A(2,2)=khi(2,2)*P(2)*P(2)/2;
    return A;
}

Vecteur coeff_V(Point *P) {
    Vecteur V(2); // taille 2x1
    V(1,1)=(khi(1,1)+khi(2,1)/2-r)*P(1); // P(1)=x1, voir comment ca s'ecrit
    V(2,1)=(khi(2,2)+khi(1,2)/2-r)*P(2);
    return V;
}

//CALCUL DE MATRICES ELEMENTAIRES M
Matrice mat_elem_M (Point * P1,Point * P2,Point * P3){
    double x1 = P1(1); double y1 = P1(2); //REMPLACER SYNTAXE
    double x2 = P2(1); double y2 = P2(2);
    double x3 = P3(1); double y3 = P3(2);
    D = ((x2-x1)*(y3-y1) - (y2-y1)*(x3-x1));
    if (abs(D) <= eps) { //fixer un eps petit pour verifier que le triangle n'est pas plat
        exit('l aire d un triangle est nulle'); 
    };
    M_elem=Matrice(3); //creation matrice 3x3 elementaire

    double Abs_Det_Bl = abs(D);

    M_des_w=Matrice(3);
    M_des_w(1,1) = 1/3 ; M_des_w(1,2) = 1/6; M_des_w(1,3) = 1/6;
    M_des_w(2,1) = 1/6; M_des_w(2,2) = 2/3; M_des_w(2,3) = 1/6;
    M_des_w(3,1) = 1/6; M_des_w(3,2) = 1/6; M_des_w(3,3) = 2/3;

    for (int i=1;i<=3;++i){
        for (int j=1;j<=3;++j){
            Terme_1 = (1/6)*M_des_w(i,1)*M_des_w(j,1);
            Terme_2 = (1/6)*M_des_w(i,2)*M_des_w(j,2);
            Terme_3 = (1/6)*M_des_w(i,3)*M_des_w(j,3);
		    M_elem(i,j) = (Terme_1 + Terme_2 + Terme_3)*Abs_Det_Bl;
        };
    };
return M_elem;
};

Matrice mat_elem_K (Point * P1,Point * P2,Point * P3) {
    double x1 = P1(1); double y1 = P1(2); //REMPLACER SYNTAXE
    double x2 = P2(1); double y2 = P2(2);
    double x3 = P3(1); double y3 = P3(2);
    double D = ((x2-x1)*(y3-y1) - (y2-y1)*(x3-x1));
    if (abs(D) <= eps) { //fixer un eps petit pour verifier que le triangle n'est pas plat
        exit('l aire d un triangle est nulle'); 
    };
    Matrice Kel(3);

    Matrice Bl(2);
    Bl(1,1)= x2-x1; Bl(1,2)= x3-x1;
    Bl(2,1)= y2-y1; Bl(2,2)= y3-y1;

    double Abs_Det_Bl = abs(D);

    Matrice inv_Bl(2); // à transposer
    inv_Bl(1,1)= Bl(2,2)/D; inv_Bl(1,2)= -Bl(1,2)/D;
    inv_Bl(2,1)= -Bl(2,1)/D; inv_Bl(1,1)= Bl(1,1)/D;

    Matrice Id(2);
    Id(1,1)= 1; Id(1,2)= 0;
    Id(2,1)= 0; Id(2,2)= 1;
    // utiliser Reftri qui est une liste qui sert de reference aux triangles
    Vecteur S(2);
    S(1)=x1; S(2)=y1;

    Vecteur S_chap_1(2);
    S_chap_1(1) = 1/6; S_chap_1(2) = 1/6;
    Vecteur S_chap_2(2);
    S_chap_2(1) = 2/3; S_chap_2(2) = 1/6;    
    Vecteur S_chap_3(2);
    S_chap_3(1) = 1/6; S_chap_3(2) = 2/3;

    Matrice Fl_1 = Bl*S_chap_1 + S; // par copie
    Matrice Fl_2 = Bl*S_chap_2 + S;
    Matrice Fl_3 = Bl*S_chap_3 + S;

    Point Fl__1=(Fl_1(1),Fl_1(2));
    Point Fl__2=(Fl_2(1),Fl_2(2));
    Point Fl__3=(Fl_3(1),Fl_3(2));

    Matrice Mat_des_gradients_w(2,3);
    Mat_des_gradients_w(1,1)= -1; Mat_des_gradients_w(1,2)= -1;
    Mat_des_gradients_w(2,1)= 1; Mat_des_gradients_w(2,2)= 0;
    Mat_des_gradients_w(3,1)= 0; Mat_des_gradients_w(3,2)= 1;

    for (int i=1;i<=3;++i) {
        for (int j=1;j<=3;++j) {
            double Terme_1 = (1/6)*condition_initiale(Fl__1)*(Inv_Bl*Mat_des_gradients_w(:,i))*(Inv_Bl*Mat_des_gradients_w(:,j))*Abs_Det_Bl;
            double Terme_2 = (1/6)*sigma_1(Fl__2)*(Inv_Bl*Mat_des_gradients_w(:,i))*(Inv_Bl*Mat_des_gradients_w(:,j))*Abs_Det_Bl;
            double Terme_3 = (1/6)*sigma_1(Fl__3)*(Inv_Bl*Mat_des_gradients_w(:,i))*(Inv_Bl*Mat_des_gradients_w(:,j))*Abs_Det_Bl;
            Kel(i,j) = Terme_1 + Terme_2 + Terme_3;
        }
    }
}