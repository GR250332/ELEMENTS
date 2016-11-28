function Fe = PRI6_FE(ex,ey,ez,es)
% function Fe = PRI6_FE(ex,ey,ez,es)
%--------------------------------------------------------------------------
% PURPOSE
%    Creation of an elementary internal force vector for PRI6
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
%     Romili PAREDES
%     25-07-2016
%--------------------------------------------------------------------------
%% Numerical integration
ngp = 6;

% Points d'intégration et poids
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

%% Precalcule de quantités pour l'intégration

% Coordonnées locales dans le repère de référénce
xsi=gp(:,1);
eta=gp(:,2);
zet=gp(:,3);

% Pour la stockage des valeurs des fonctions de formes et de leurs
% gradients
r2=ngp*3;

%% Fonctions de formes et gradients

% Derivées des fonctions de formes
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

%% Set the stress
stress = es;

%% Loop over the integration points

% Jacobian matrix
JT = dNr*[ex(1,:);ey(1,:);ez(1,:)]';

% Initialization
Fe = zeros(18,1);

ir = 0;

indx = [1:3]';

for i=1:ngp
    
    ir = ir + 1;
    
    detJ = det(JT(indx,:));
    
    if detJ < 10*eps
        
        error('Jacobian determinant equal or less than zero')
        
    end
    
    dNx   = JT(indx,:) \dNr(indx,:);
    
    B(1,1:3:18-2) = dNx(1,:);
    B(2,2:3:18-1) = dNx(2,:);
    B(3,3:3:18)   = dNx(3,:);
    B(4,1:3:18-2) = dNx(2,:);
    B(4,2:3:18-1) = dNx(1,:);
    B(5,1:3:18-2) = dNx(3,:);
    B(5,3:3:18)   = dNx(1,:);
    B(6,2:3:18-1) = dNx(3,:);
    B(6,3:3:18)   = dNx(2,:);
    
    Fe = Fe + B' *es(ir,:)' * wp(i) * detJ;
    
    indx = indx + 3;
    
end