function Fe = QUA4_FE(ma,ex,ey,es)
% function Fe = QUA4_FE(ma,ex,ey,es)
%--------------------------------------------------------------------------
% PURPOSE
%    Creation of an elementary internal force vector for QUA4
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
%     Benjamin RICHARD
%     10-02-2016
%--------------------------------------------------------------------------

%% Some information

% Thickness
t = ma.cara.epai;

%% Regles d'integration numerique
ir  = 2;
ngp = ir * ir;

% Points d'intégration et poids
if ir==1
    
    g1 = 0.0;
    
    w1 = 2.0;
    
    gp = [ g1 g1 ];
    
    w  = [ w1 w1 ];
    
elseif ir==2
    
    g1 = 0.577350269189626;
    w1 = 1;
    
    gp(:,1) = [-g1; g1;-g1; g1];
    gp(:,2) = [-g1;-g1; g1; g1];
    
    w(:,1) = [ w1; w1; w1; w1];
    w(:,2) = [ w1; w1; w1; w1];
    
elseif ir==3
    
    g1 = 0.774596699241483;
    g2 = 0.;
    
    w1 = 0.555555555555555;
    w2 = 0.888888888888888;
    
    gp(:,1) = [-g1;-g2; g1;-g1; g2; g1;-g1; g2; g1];
    gp(:,2) = [-g1;-g1;-g1; g2; g2; g2; g1; g1; g1];
    
    w(:,1) = [ w1; w2; w1; w1; w2; w1; w1; w2; w1];
    w(:,2) = [ w1; w1; w1; w2; w2; w2; w1; w1; w1];
    
else
    
    error('Used number of integration points not implemented');
    
end

%% Precalcule de quantités pour l'intégration

% Matrice des poids
wp  = w(:,1).*w(:,2);

% Coordonnées locales dans le repère de référénce
xsi = gp(:,1);
eta = gp(:,2);

% Pour la stockage des valeurs des fonctions de formes et de leurs
% gradients
r2  = ngp * 2;

%% Derivées des fonctions de formes
dNr(1:2:r2,1) = -(1-eta)/4;
dNr(1:2:r2,2) =  (1-eta)/4;

dNr(1:2:r2,3) =  (1+eta)/4;
dNr(1:2:r2,4) = -(1+eta)/4;

dNr(2:2:r2+1,1) = -(1-xsi)/4;
dNr(2:2:r2+1,2) = -(1+xsi)/4;

dNr(2:2:r2+1,3) = (1+xsi)/4;
dNr(2:2:r2+1,4) = (1-xsi)/4;

%% On ne considère pas sigma_33 (soit nulle soit n'intervient pas dans Fe)
stress = es(:,[1 2 4]);

%% Boucle sur les points d'intégration

JT = dNr*[ex(1,:);ey(1,:)]';

Fe = zeros(8,1);

for i=1:ngp
    
    indx=[ 2*i-1; 2*i ];
    
    detJ = det(JT(indx,:));
    
    if detJ < 10*eps
        
        error('Jacobian determinant equal or less than zero')
        
    end
    
    dNx=JT(indx,:)\dNr(indx,:);
    
    B(1,1:2:8-1) = dNx(1,:);
    B(2,2:2:8)   = dNx(2,:);
    B(3,1:2:8-1) = dNx(2,:);
    B(3,2:2:8)   = dNx(1,:);
    
    Fe = Fe + B'*stress(i,:)'*wp(i)*detJ*t;
    
end

