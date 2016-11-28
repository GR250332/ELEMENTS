function Me = PRI6_ME(ex,ey,ez,ma,varargin)
%  function Me = PRI6_ME(ex,ey,ez,ma,type_module,varargin)
%--------------------------------------------------------------------------
% PURPOSE
%    Calculate the elementary mass matrix for a PRI6 element
%--------------------------------------------------------------------------
% INPUT
%    Ex          : Elementary nodal X coordinates
%    Ey          : Elementary nodal Y coordinates
%    Ez          : Elementary nodal Z coordinates
%    ma          : Matarial structure
%    varargin    : element fields
%--------------------------------------------------------------------------
% OUTPUT
%    Me          : local stiffness matrix
%--------------------------------------------------------------------------
% REFERENCES
%     Romili PAREDES
%     25-07-2016
%--------------------------------------------------------------------------
%% Ouverture variables globales
global options

%% Regles d'integration numerique
ngp = 6;

% Points d'int�gration et poids
if ngp == 6
    
    g1 = 1/sqrt(3);
    g2 = 1/2;
    g3 = 0;
    
    w1 = 1/6;
    
    gp(:,1)=[-g1; -g1; -g1; g1; g1; g1];
    gp(:,2)=[g2; g3; g2; g2; g3; g2];
    gp(:,3)=[g2; g2; g3; g2; g2; g3];
    
    wp = [w1; w1; w1; w1; w1; w1;];
    
elseif ngp == 8
    
    g1 = 0.577350269189626;
    g2 = 1/3;
    g3 = 0.6;
    g4 = 0.2;
    
    w1 = -27/96;
    w2 = 25/96;
    
    gp(:,1)=[-g1; -g1; -g1; -g1; g1; g1; g1; g1;];
    gp(:,2)=[g2; g3; g4; g4; g2; g3; g4; g4];
    gp(:,3)=[g2; g4; g3; g4; g2; g4; g3; g4];
    
    wp = [w1; w2; w2; w2; w1; w2; w2; w2];
    
else
    
    error('Used number of integration points not implemented');
    
end;

%% Precalcule de quantit�s pour l'int�gration

% Coordonn�es locales dans le rep�re de r�f�r�nce
xsi = gp(:,1);
eta = gp(:,2);
zet = gp(:,3);

% Pour la stockage des valeurs des fonctions de formes et de leurs
% gradients
r2=ngp*3;

%% Fonctions de formes et gradients

% Fonction de forme
N(:,2) = eta.*(1 - xsi);
N(:,3) = zet.*(1 - xsi);
N(:,1) = (1 - eta - zet).*(1 - xsi);
N(:,5) = eta.*(xsi + 1);
N(:,6) = zet.*(xsi + 1);
N(:,4) = (1 - eta - zet).*(xsi + 1);
N      = N/2; 

% Deriv�es des fonctions de formes
dNr(1:3:r2,2)= -eta/2;
dNr(1:3:r2,3)= -zet/2;
dNr(1:3:r2,1)= eta/2 + zet/2 - 1/2;
dNr(1:3:r2,5)= eta/2;
dNr(1:3:r2,6)= zet/2;
dNr(1:3:r2,4)= 1/2 - zet/2 - eta/2;
dNr(2:3:r2+1,2)= 1/2 - xsi/2;
dNr(2:3:r2+1,3)= 0;
dNr(2:3:r2+1,1)= xsi/2 - 1/2;
dNr(2:3:r2+1,5)= xsi/2 + 1/2;
dNr(2:3:r2+1,6)= 0;
dNr(2:3:r2+1,4)= -xsi/2 - 1/2;
dNr(3:3:r2+2,2)= 0;
dNr(3:3:r2+2,3)= 1/2 - xsi/2;
dNr(3:3:r2+2,1)= xsi/2 - 1/2;
dNr(3:3:r2+2,5)= 0;
dNr(3:3:r2+2,6)= xsi/2 + 1/2;
dNr(3:3:r2+2,4)= -xsi/2 - 1/2;

%% Initialisation avant integration num�rique

% Matrice elementaire
Me=sparse(18,18);

% Jacobian
JT=dNr*[ex;ey;ez]';

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
        
        if detJ<10*eps
            
            error('Jacobian determinant equal or less than zero')
            
        end
        
        dNx=JT(indx,:)\dNr(indx,:);
        
        B(1,1:3:18-2) = dNx(1,:);
        B(2,2:3:18-1) = dNx(2,:);
        B(3,3:3:18)   = dNx(3,:);
        B(4,1:3:18-2) = dNx(2,:);
        B(4,2:3:18-1) = dNx(1,:);
        B(5,1:3:18-2) = dNx(3,:);
        B(5,3:3:18)   = dNx(1,:);
        B(6,2:3:18-1) = dNx(3,:);
        B(6,3:3:18)   = dNx(2,:);
        
        N2(1,1:3:18-2) = N(i,:);
        N2(2,2:3:18-1) = N(i,:);
        N2(3,3:3:18)   = N(i,:);
        
        Me=Me+rho*(N2')*N2*detJ*wp(i);
        
    end
end
