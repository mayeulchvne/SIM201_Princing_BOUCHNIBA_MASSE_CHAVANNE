%% affiche le maillage au deffirents instant
function affiche_video(UU,Numtri,Coorneu, titre)

% control on the input args
if (nargin<4), titre = ''; end

figure;
trisurf(Numtri,Coorneu(:,1),Coorneu(:,2),UU);
view(2);
shading interp
% shading faceted
% shading flat
colorbar;

% ajouter eventuellement un titre
title(titre);

