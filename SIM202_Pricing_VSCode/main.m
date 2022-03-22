% lecture du maillage et affichage
% ---------------------------------
nom_maillage = "test.txt";
[video, Nbpt,Nbtri,Numtri,Coorneu]= lecture_cppp(nom_maillage);

% lecture du maillage et affichage
% ---------------------------------
affiche_video(video(:,1) , Numtri, Coorneu, "instant initial");



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%