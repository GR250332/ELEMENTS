function Me = BAR2_ME(ex,ey,ma,varargin)
% function Me = BAR2_ME(ex,ey,ma,varargin)
%--------------------------------------------------------------------------
% PURPOSE
%    Calculate the elementary mass matrix for a BARRE element
%--------------------------------------------------------------------------
% INPUT   
%    ex          : Elementary nodal X coordinates
%    ey          : Elementary nodal Y coordinates
%    ma          : Matarial structure
%    varargin    : element fields
%--------------------------------------------------------------------------
% OUTPUT
%    Me          : local mass matrix
%--------------------------------------------------------------------------
% REFERENCES
%     Benjamin RICHARD
%     26-02-2016
%--------------------------------------------------------------------------

%% Regles d'integration numerique
ngp = 2;

% Points d'intégration et poids
if ngp == 2
    
    g1 = 0.577350269189626; 
    w1 = 1;
    
    gp = [-g1 ; g1]; 
    
    w  = [w1 ; w1];
    
    wp = w(:,1);
        
else
    
    error('Used number of integration points not implemented');
    
end;

%% Coordonnées locales dans le repère de référénce
xsi = gp(:,1);

% Pour la stockage des valeurs des fonctions de formes et de leurs
% gradients
r2 = ngp * 1;

%% Fonctions de formes et gradients

% Fonction de forme
N(:,1) = 0.5 * (1 - xsi);
N(:,2) = 0.5 * (1 + xsi);

%% Initialisation avant integration numérique

% Matrice elementaire
Me = sparse(4,4);

%% Node coordinates
b = [ex(2)-ex(1) ; ey(2)-ey(1)];

%% Length
L = sqrt(b'*b);

% Some information
t   = ma.cara.epai;
rho = ma.cara.rho;
A   = ma.cara.sect;

% Boucle sur les points de Gauss
for i=1:ngp
        
    detJ = L/2;
    
    if detJ < 10*eps
        
        error('Jacobian determinant equal or less than zero')
        
    end
    
    N2(1,1:2:4-1) = N(i,:);
    N2(2,2:2:4) = N(i,:);
    
    Me=Me+rho*(N2')*N2*detJ*t*A*wp(i);
    
end