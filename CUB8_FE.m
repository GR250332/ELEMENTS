function Fe = CUB8_FE(ex,ey,ez,es)
% function Fe = CUB8_FE(ex,ey,ez,es)
%--------------------------------------------------------------------------
% PURPOSE
%    Creation of an elementary internal force vector for CUB8
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
%     10-02-2016
%--------------------------------------------------------------------------
%% Numerical integration
ir  = 2;
ngp = ir*ir*ir;

% Points d'intégration et poids
if ir==1
    
    g1=0.0; w1=2.0;
    gp=[ g1 g1 ];  w=[ w1 w1 ];
    
elseif ir==2
    
    g1=0.577350269189626;
    
    w1=1;
    
    gp(:,1)=[-1; 1; 1;-1;-1; 1; 1;-1]*g1;
    gp(:,2)=[-1;-1; 1; 1;-1;-1; 1; 1]*g1;
    gp(:,3)=[-1;-1;-1;-1; 1; 1; 1; 1]*g1;
    
    w(:,1)=[ 1; 1; 1; 1; 1; 1; 1; 1]*w1;
    w(:,2)=[ 1; 1; 1; 1; 1; 1; 1; 1]*w1;
    w(:,3)=[ 1; 1; 1; 1; 1; 1; 1; 1]*w1;
    
elseif ir==3
    
    g1=0.774596699241483;
    g2=0.;
    
    w1=0.555555555555555;
    w2=0.888888888888888;
    
    I1=[-1; 0; 1;-1; 0; 1;-1; 0; 1]';
    I2=[ 0;-1; 0; 0; 1; 0; 0; 1; 0]';
    
    gp(:,1)=[I1 I1 I1]'*g1;
    gp(:,1)=[I2 I2 I2]'*g2+gp(:,1);
    
    I1=abs(I1);
    I2=abs(I2);
    w(:,1)=[I1 I1 I1]'*w1;
    w(:,1)=[I2 I2 I2]'*w2+w(:,1);
    
    I1=[-1;-1;-1; 0; 0; 0; 1; 1; 1]';
    I2=[ 0; 0; 0; 1; 1; 1; 0; 0; 0]';
    gp(:,2)=[I1 I1 I1]'*g1;
    gp(:,2)=[I2 I2 I2]'*g2+gp(:,2);
    
    I1=abs(I1); I2=abs(I2);
    w(:,2)=[I1 I1 I1]'*w1;
    w(:,2)=[I2 I2 I2]'*w2+w(:,2);
    
    I1=[-1;-1;-1;-1;-1;-1;-1;-1;-1]';
    I2=[ 0; 0; 0; 0; 0; 0; 0; 0; 0]';
    I3=abs(I1);
    gp(:,3)=[I1 I2 I3]'*g1;
    gp(:,3)=[I2 I3 I2]'*g2+gp(:,3);
    
    w(:,3)=[I3 I2 I3]'*w1;
    w(:,3)=[I2 I3 I2]'*w2+w(:,3);
    
else
    
    error('Used number of integration points not implemented');
    
end;

%% Precalcule de quantités pour l'intégration
% Matrice des poids
wp=w(:,1).*w(:,2).*w(:,3);

% Coordonnées locales dans le repère de référénce
xsi=gp(:,1);
eta=gp(:,2);
zet=gp(:,3);

% Pour la stockage des valeurs des fonctions de formes et de leurs
% gradients
r2=ngp*3;

%% Fonctions de formes et gradients

% Fonction de forme
N(:,1)=(1-xsi).*(1-eta).*(1-zet)/8;
N(:,5)=(1-xsi).*(1-eta).*(1+zet)/8;
N(:,2)=(1+xsi).*(1-eta).*(1-zet)/8;
N(:,6)=(1+xsi).*(1-eta).*(1+zet)/8;
N(:,3)=(1+xsi).*(1+eta).*(1-zet)/8;
N(:,7)=(1+xsi).*(1+eta).*(1+zet)/8;
N(:,4)=(1-xsi).*(1+eta).*(1-zet)/8;
N(:,8)=(1-xsi).*(1+eta).*(1+zet)/8;

% Derivées des fonctions de formes

dNr(1:3:r2,1)=-(1-eta).*(1-zet);
dNr(1:3:r2,2)= (1-eta).*(1-zet);
dNr(1:3:r2,3)= (1+eta).*(1-zet);
dNr(1:3:r2,4)=-(1+eta).*(1-zet);
dNr(1:3:r2,5)=-(1-eta).*(1+zet);
dNr(1:3:r2,6)= (1-eta).*(1+zet);
dNr(1:3:r2,7)= (1+eta).*(1+zet);
dNr(1:3:r2,8)=-(1+eta).*(1+zet);
dNr(2:3:r2+1,1)=-(1-xsi).*(1-zet);
dNr(2:3:r2+1,2)=-(1+xsi).*(1-zet);
dNr(2:3:r2+1,3)= (1+xsi).*(1-zet);
dNr(2:3:r2+1,4)= (1-xsi).*(1-zet);
dNr(2:3:r2+1,5)=-(1-xsi).*(1+zet);
dNr(2:3:r2+1,6)=-(1+xsi).*(1+zet);
dNr(2:3:r2+1,7)= (1+xsi).*(1+zet);
dNr(2:3:r2+1,8)= (1-xsi).*(1+zet);
dNr(3:3:r2+2,1)=-(1-xsi).*(1-eta);
dNr(3:3:r2+2,2)=-(1+xsi).*(1-eta);
dNr(3:3:r2+2,3)=-(1+xsi).*(1+eta);
dNr(3:3:r2+2,4)=-(1-xsi).*(1+eta);
dNr(3:3:r2+2,5)= (1-xsi).*(1-eta);
dNr(3:3:r2+2,6)= (1+xsi).*(1-eta);
dNr(3:3:r2+2,7)= (1+xsi).*(1+eta);
dNr(3:3:r2+2,8)= (1-xsi).*(1+eta);
dNr=dNr/8.;

%% Set the stress
stress = es;

%% Loop over the integration points

% Jacobian matrix
JT = dNr*[ex(1,:);ey(1,:);ez(1,:)]';

% Initialization
Fe = zeros(24,1);

ir = 0;

indx = [1:3]';

for i=1:ngp
    
    ir = ir + 1;
    
    detJ = det(JT(indx,:));
    
    if detJ < 10*eps
        
        error('Jacobian determinant equal or less than zero')
        
    end
    
    dNx   = JT(indx,:) \dNr(indx,:);
    
    B(1,1:3:24-2) = dNx(1,:);
    B(2,2:3:24-1) = dNx(2,:);
    B(3,3:3:24)   = dNx(3,:);
    B(4,1:3:24-2) = dNx(2,:);
    B(4,2:3:24-1) = dNx(1,:);
    B(5,1:3:24-2) = dNx(3,:);
    B(5,3:3:24)   = dNx(1,:);
    B(6,2:3:24-1) = dNx(3,:);
    B(6,3:3:24)   = dNx(2,:);
    
    Fe = Fe + B' *es(ir,:)' * wp(i) * detJ;
    
    indx = indx + 3;
    
end