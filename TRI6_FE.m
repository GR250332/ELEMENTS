function Fe = TRI6_FE(ma,ex,ey,es)
% function Fe = TRI6_FE(ma,ex,ey,es)
%--------------------------------------------------------------------------
% PURPOSE
%    Creation of an elementary internal force vector for TRI6
%--------------------------------------------------------------------------
% INPUT
%    ma       : material structutre
%    ex       : Elementary nodal X coordinates
%    ey       : Elementary nodal Y coordinates
%    es       : Stress at the integration point
%--------------------------------------------------------------------------
% OUTPUT
%    Fe       : Elementary internal force vector
%
%--------------------------------------------------------------------------
% REFERENCES
%     Romili PAREDES
%     25-07-2016
%--------------------------------------------------------------------------

%% Some information

% Thickness
t = ma.cara.epai;

%% Regles d'integration numerique
ngp = 4;

% Points d'intégration et poids
if ngp == 1
    
    g1  = 1/3;
    w1  = 1/2;
    
    gp  = [g1 g1];
    
    wp  = w1;
    
elseif ngp == 3
    
    g1 = 1/6;
    g2 = 2/3;
    
    w1 = 1/6;
    
    gp(:,1) = [g1; g2; g1];
    gp(:,2) = [g1; g1; g2];
    
    wp = [ w1; w1; w1];
    
elseif ngp == 4
    
    g1 = 1/5;
    g2 = 3/5;
    g3 = 1/3; 
    
    w1 = 25/(24*4);
    w2 = -27/(24*4);
    
    gp(:,1) = [g1; g2; g1; g3];
    gp(:,2) = [g1; g1; g2; g3];
    
    wp = [w1; w1; w1; w2];
    
elseif ngp == 7
    
    g1 = 1/3;
    A  = 0.470142064105115;
    B  = 0.101286507323456;
    
    w1 = 9/80;
    w2 = 0.066197076394253;
    w3 = 0.062969590272413;
    
    gp(:,1) = [g1; A; 1-2*A; A; B; 1-2*B; B];
    gp(:,2) = [g1; A; A; 1-2*A; B; B; 1-2*B];
    
    wp = [w1; w2; w2; w2; w3; w3; w3];
    
else
    
    error('Used number of integration points not implemented');
    
end;

%% Precalcule de quantités pour l'intégration

% Coordonnées locales dans le repère de référénce
xsi = gp(:,1);
eta = gp(:,2);

% Pour la stockage des valeurs des fonctions de formes et de leurs
% gradients
r2  = ngp * 2;

% Derivées des fonctions de formes
dNr(1:2:r2,1) = 4.*xsi + 4.*eta - 3;
dNr(1:2:r2,3) = 4.*xsi - 1;
dNr(1:2:r2,5) = 0;
dNr(1:2:r2,2) = 4 - 8.*xsi - 4.*eta;
dNr(1:2:r2,4) = 4.*eta;
dNr(1:2:r2,6) = -4.*eta;

dNr(2:2:r2+1,1) = 4.*xsi + 4.*eta - 3;
dNr(2:2:r2+1,3) = 0;
dNr(2:2:r2+1,5) = 4.*eta - 1;
dNr(2:2:r2+1,2) = -4.*xsi;
dNr(2:2:r2+1,4) = 4.*xsi;
dNr(2:2:r2+1,6) = 4 - 4.*xsi - 8.*eta;

%% On ne considère pas sigma_33 (soit nulle soit n'intervient pas dans Fe)
stress = es(:,[1 2 4]);

%% Boucle sur les points d'intégration

JT = dNr*[ex(1,:);ey(1,:)]';

Fe = zeros(12, 1);

for i=1:ngp
    
    indx=[ 2*i-1; 2*i ];
    
    detJ = det(JT(indx,:));
    
    if detJ < 10*eps
        
        error('Jacobian determinant equal or less than zero')
        
    end
    
    dNx=JT(indx,:)\dNr(indx,:);
    
    B(1,1:2:12-1) = dNx(1,:);
    B(2,2:2:12)   = dNx(2,:);
    B(3,1:2:12-1) = dNx(2,:);
    B(3,2:2:12)   = dNx(1,:);
    
    Fe = Fe + B'*stress(i,:)'*wp(i)*detJ*t;
    
end

