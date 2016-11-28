function Me = TET4_ME(ex,ey,ez,ma,varargin)
%  function Me = TET4_ME(ex,ey,ez,ma,varargin)
%--------------------------------------------------------------------------
% PURPOSE
%    Calculate the elementary mass matrix for a TET4 element
%--------------------------------------------------------------------------
% INPUT
%    Ex          : Elementary nodal X coordinates
%    Ey          : Elementary nodal Y coordinates
%    Ez          : Elementary nodal Z coordinates
%    ma          : Matarial structure
%    varargin    : element fields
%--------------------------------------------------------------------------
% OUTPUT
%    Me          : local mass matrix
%--------------------------------------------------------------------------
% REFERENCES
%     Benjamin RICHARD
%     23-02-2016
%--------------------------------------------------------------------------
%% Ouverture variables globales
global options

%% Regles d'integration numerique
ngp = 4;

% Points d'intégration et poids
if ngp == 1
    
    g1=1/4; 
    w1=1/6;
    
    gp=[g1 g1 g1];  
    w=[w1 w1 w1];
    
    wp = w(:,1);
    
elseif ngp == 4
    
    w1 = 1/24;
    
    a  = (5 -   sqrt(5)) / 20;
    b  = (5 + 3*sqrt(5)) / 20;
    
    gp(:,1) = [a; a; a; b];
    gp(:,2) = [a; a; b; a];
    gp(:,3) = [a; b; a; a];

    w(:,1) = [ w1; w1; w1; w1];
    
    wp=w(:,1)';

else
    
    error('Used number of integration points not implemented');
    
end;

% Coordonnées locales dans le repère de référénce
xsi=gp(:,1);
eta=gp(:,2);
zet=gp(:,3);

% Pour la stockage des valeurs des fonctions de formes et de leurs
% gradients
r2=ngp*3;

%% Fonctions de formes et gradients

% Fonction de forme
N(:,1) = 1 - xsi - eta - zet;
N(:,2) = xsi;
N(:,3) = eta;
N(:,4) = zet;

% Derivées des fonctions de formes
dNr(1:3:r2,1)=-1;
dNr(1:3:r2,2)= 1;
dNr(1:3:r2,3)= 0;
dNr(1:3:r2,4)= 0;

dNr(2:3:r2+1,1)=-1;
dNr(2:3:r2+1,2)= 0;
dNr(2:3:r2+1,3)= 1;
dNr(2:3:r2+1,4)= 0;

dNr(3:3:r2+2,1)=-1;
dNr(3:3:r2+2,2)= 0;
dNr(3:3:r2+2,3)= 0;
dNr(3:3:r2+2,4)= 1;

%% Initialisation avant integration numérique

% Matrice elementaire
Me = sparse(12,12);

% Jacobian
JT = dNr*[ex;ey;ez]';

% Masse volumique
rho = ma.cara.rho;

%% Mode TRID uniquement
if ~strcmp(options.mode,'TRID')
    
    error('Computation mode is not correct')
    
else
    
    % Boucle sur les points de Gauss
    for i=1:ngp
        
        indx=[ 3*i-2; 3*i-1; 3*i ];
        
        detJ=det(JT(indx,:));
        
        if detJ < 10*eps
            
            error('Jacobian determinant equal or less than zero')
            
        end
                
        N2(1,1:3:12-2) = N(i,:);
        N2(2,2:3:12-1) = N(i,:);
        N2(3,3:3:12)   = N(i,:);
        
        Me=Me+rho*(N2')*N2*detJ*wp(i);
        
    end
end
