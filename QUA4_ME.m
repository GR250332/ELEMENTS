function Me = QUA4_ME(ex,ey,ma,varargin)
% function Me = QUA4_ME(ex,ey,ma,type_module,varargin)
%--------------------------------------------------------------------------
% PURPOSE
%    Calculate the elementary mass matrix for a QUA4 element
%--------------------------------------------------------------------------
% INPUT
%    Ex          : Elementary nodal X coordinates
%    Ey          : Elementary nodal Y coordinates
%    ma          : Matarial structure
%    varargin    : element fields
%--------------------------------------------------------------------------
% OUTPUT
%    Ke          : local stiffness matrix
%--------------------------------------------------------------------------
% REFERENCES
%     Benjamin RICHARD
%     11-02-2016
%--------------------------------------------------------------------------

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

%% Fonctions de formes et gradients

% Fonction de forme
N(:,1)=(1-xsi).*(1-eta)/4;
N(:,2)=(1+xsi).*(1-eta)/4;
N(:,3)=(1+xsi).*(1+eta)/4;
N(:,4)=(1-xsi).*(1+eta)/4;

% Derivées des fonctions de formes
dNr(1:2:r2,1) = -(1-eta)/4;
dNr(1:2:r2,2) =  (1-eta)/4;

dNr(1:2:r2,3) =  (1+eta)/4;
dNr(1:2:r2,4) = -(1+eta)/4;

dNr(2:2:r2+1,1) = -(1-xsi)/4;
dNr(2:2:r2+1,2) = -(1+xsi)/4;

dNr(2:2:r2+1,3) = (1+xsi)/4;
dNr(2:2:r2+1,4) = (1-xsi)/4;

%% Initialisation avant integration numérique

% Matrice elementaire
Me = sparse(8,8);

% Jacobian
JT = dNr*[ex;ey]';

% Caracteristiques
t   = ma.cara.epai;
rho = ma.cara.rho;

%% Boucle sur les points de Gauss
for i=1:ngp
    
    indx=[2*i-1; 2*i ];
    
    detJ=det(JT(indx,:));
    
    if detJ < 10*eps
        
        error('Jacobian determinant equal or less than zero')
        
    end
    
    N2(1,1:2:8-1) = N(i,:);
    N2(2,2:2:8)   = N(i,:);
    
    Me = Me + rho*(N2')*detJ*N2*wp(i);
    
end

