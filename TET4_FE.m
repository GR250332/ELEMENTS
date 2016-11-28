function Fe = TET4_FE(ex,ey,ez,es)
% function Fe = TET4_FE(ex,ey,ez,es)
%--------------------------------------------------------------------------
% PURPOSE
%    Creation of an elementary internal force vector for TET4
%--------------------------------------------------------------------------
% INPUT
%    ex       : Elementary nodal X coordinates
%    ey       : Elementary nodal Y coordinates
%    es       : Stress at the integration points
%--------------------------------------------------------------------------
% OUTPUT
%    Fe       : Elementary internal force vector
%
%--------------------------------------------------------------------------
% REFERENCES
%     Benjamin RICHARD
%     24-02-2016
%--------------------------------------------------------------------------
%% Regles d'integration numerique
ngp = 1;

% Points d'intégration et poids
if ngp == 1
    
    g1 = 1/4; 
    w1 = 1/6;
    
    gp = [g1 g1 g1]; 
    
    w  = [w1 w1 w1];
    
    wp = w(:,1);
    
elseif ngp == 4
    
    w1 = 1/24;
    
    a  = (5 -   sqrt(5)) / 20;
    b  = (5 + 3*sqrt(5)) / 20;
    
    gp(:,1) = [a; a; a; b];
    gp(:,2) = [a; a; b; a];
    gp(:,3) = [a; b; a; a];

    w(:,1) = [ w1; w1; w1; w1];
    w(:,2) = [ w1; w1; w1; w1];
    w(:,3) = [ w1; w1; w1; w1];
    
    wp = w(:,1).*w(:,2).*w(:,3);
    
else
    
    error('Used number of integration points not implemented');
    
end;

%% Coordonnées locales dans le repère de référénce
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

%% Loop over the integration points

% Jacobian matrix
JT = dNr*[ex(1,:);ey(1,:);ez(1,:)]';

% Initialization
Fe = zeros(12,1);

ir = 0;

indx = [1:3]';

for i=1:ngp
    
    ir = ir + 1;
    
    detJ = det(JT(indx,:));
    
    if detJ < 10*eps
        
        error('Jacobian determinant equal or less than zero')
        
    end
    
    dNx   = JT(indx,:) \dNr(indx,:);
    
    B(1,1:3:12-2) = dNx(1,:);
    B(2,2:3:12-1) = dNx(2,:);
    B(3,3:3:12)   = dNx(3,:);
    
    B(4,1:3:12-2) = dNx(2,:);
    B(4,2:3:12-1) = dNx(1,:);
    B(5,1:3:12-2) = dNx(3,:);
    
    B(5,3:3:12)   = dNx(1,:);
    B(6,2:3:12-1) = dNx(3,:);
    B(6,3:3:12)   = dNx(2,:);
    
    Fe = Fe + B' *es(ir,:)' * wp(i) * detJ;
    
    indx = indx + 3;
    
end