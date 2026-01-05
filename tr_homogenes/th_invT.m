 function kTi=th_invT(iTk)
    kOi=iTk(1:3,4);   % origine repere j exprimee dans repere i= iTk, lignes 1 a 3, colonne 4
    iRk=iTk(1:3,1:3); % iRk=vecteurs [xj,yj,zj] exprimes suivant [xi,yi,zi]=iTk lignes 1 a 3, colonnes 1 a 3
    kRi=iRk.';         % Rji=vecteurs [xi,yi,zi] exprimes suivant [xj,yj,zj] = transposee de iRk
    kOi=-kRi*kOi;     % Oji = origine repere i exprimee dans repere j = Rji . [-Oij]
    kTi=[[ kRi,kOi];[0,0,0,1 ] ];
  end
