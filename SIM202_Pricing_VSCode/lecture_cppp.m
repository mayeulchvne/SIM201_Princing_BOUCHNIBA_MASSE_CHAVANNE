%% Fonction de lecture du fichier texte edité par le script main.cpp
function [video,Nbpt,Nbtri,Numtri,Coorneu]=lecture_cppp(nomfile)
fid=fopen(nomfile,'r');
if fid <=0,
   msg=['Le fichier de maillage : ' nomfile ' n''a pas ete trouve'];
   error(msg);
end
Nbpt = str2int(fgetl(fid));
Coorneu = zeros(Nbpt,2);

tmp= str2double(fgetl(fid));
for i=1:Nbpt
    Coorneu(i,:) = tmp(end-1:end);
end


Nbtri = str2int(fgetl(fid));
Numtri = zeros(Nbtri,3);

tmp= str2int(fgetl(fid));
for i=1:Nbtri
    Numtri(i,:) = tmp(end-2:end);
end


temps=str2int(fgetl(fid));
video=zeros(Nbpt,temps);
valeurs=zeros(Nbpt,1);
tmp=str2double(fgetl(fid));

for j=1:temps
    for i=1:Nbvaleurs
        valeurs(i)=tmp(end);
        tmp=str2double(fgetl(fid));
    end
   
    video(:,j)=valeurs(:);
end

fclose(fid);