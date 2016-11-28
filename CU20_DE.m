function strain = CU20_DE(ex,ey,ez,ed)
% function strain = CU20_DE(ex,ey,ez,ed)
%--------------------------------------------------------------------------
% PURPOSE
%    Computation of the strain from a displacement vector ed for a CU20
%--------------------------------------------------------------------------
% INPUT
%    ex : X coordinates
%    ey : Y coordinates
%    ez : Z coordinates
%    ed : Nodal displacements
%--------------------------------------------------------------------------
% OUTPUT
%    strain
%--------------------------------------------------------------------------
% REFERENCES
%     Romili PAREDES
%     25-07-2016
%--------------------------------------------------------------------------

%% Position of the integration points
ngp = 27;

if ngp == 8
    
    g1 = 1/sqrt(3);
    
    gp(:,1)=[-g1; -g1; -g1; -g1; g1; g1; g1; g1];
    gp(:,2)=[-g1; -g1; g1; g1; -g1; -g1; g1; g1];
    gp(:,3)=[-g1; g1; -g1; g1; -g1; g1; -g1; g1];

elseif ngp == 27
    
    g1 = sqrt(3/5);
    g2 = 0;
    
    gp(:,1)=[-g1; -g1; -g1; -g1; -g1; -g1; -g1; -g1; -g1; ...
             g2; g2; g2; g2; g2; g2; g2; g2; g2; ...
             g1; g1; g1; g1; g1; g1; g1; g1; g1];
    gp(:,2)=[-g1; -g1; -g1; g2; g2; g2; g1; g1; g1; ...
             -g1; -g1; -g1; g2; g2; g2; g1; g1; g1; ...
             -g1; -g1; -g1; g2; g2; g2; g1; g1; g1];
    gp(:,3)=[-g1; g2; g1; -g1; g2; g1; -g1; g2; g1; ...
             -g1; g2; g1; -g1; g2; g1; -g1; g2; g1; ...
             -g1; g2; g1; -g1; g2; g1; -g1; g2; g1];
    
else
    
    error('Used number of integration points not implemented');
    
end;

%% Coordinates of the integration points in the reference elements
xsi = gp(:,1);
eta = gp(:,2);
zet = gp(:,3);
r2  = ngp*3;

%% Fonctions de formes et gradients

% Coordonnées dans l'espace de réference
ccr(:,1) = [-1; -1; -1; -1; -1; -1; -1; -1; 0; 0; 0; 0; 1; 1; 1; 1; 1; 1; 1; 1];
ccr(:,2) = [1; 0; -1; -1; -1; 0; 1; 1; 1; -1; -1; 1; 1; 0; -1; -1; -1; 0; 1; 1];
ccr(:,3) = [1; 1; 1; 0; -1; -1; -1; 0; 1; 1; -1; -1; 1; 1; 1; 0; -1; -1; -1; 0];

% Derivées des fonctions de formes

dNr = zeros(r2,20);

for i = [5 17 19 7 3 15 13 1]
	dNr(1:3:r2,i)   = (1/8).*ccr(i,1).*(1 + ccr(i,2).*eta).*(1 + ccr(i,3).*zet).*(2.*ccr(i,1).*xsi + ccr(i,2).*eta + ccr(i,3).*zet - 1);
    dNr(2:3:r2+1,i) = (1/8).*(1 + ccr(i,1).*xsi).*ccr(i,2).*(1 + ccr(i,3).*zet).*(ccr(i,1).*xsi + 2.*ccr(i,2).*eta + ccr(i,3).*zet - 1);
    dNr(3:3:r2+2,i) = (1/8).*(1 + ccr(i,1).*xsi).*(1 + ccr(i,2).*eta).*ccr(i,3).*(ccr(i,1).*xsi + ccr(i,2).*eta + 2.*ccr(i,3).*zet - 1);
end

for i = [11 12 10 9]
    dNr(1:3:r2,i)   = (1/4).*-2.*xsi.*(1 + ccr(i,2).*eta).*(1 + ccr(i,3).*zet);
    dNr(2:3:r2+1,i) = (1/4).*(1 - xsi.^2).*ccr(i,2).*(1 + ccr(i,3).*zet);
    dNr(3:3:r2+2,i) = (1/4).*(1 - xsi.^2).*(1 + ccr(i,2).*eta).*ccr(i,3);
end

for i = [18 6 14 2]
    dNr(1:3:r2,i)   = (1/4).*(1 - eta.^2).*ccr(i,1).*(1 + ccr(i,3).*zet);
    dNr(2:3:r2+1,i) = (1/4).*-2.*eta.*(1 + ccr(i,1).*xsi).*(1 + ccr(i,3).*zet);
    dNr(3:3:r2+2,i) = (1/4).*(1 - eta.^2).*(1 + ccr(i,1).*xsi).*ccr(i,3);
end

for i = [4 16 20 8]
    dNr(1:3:r2,i)   = (1/4).*(1 - zet.^2).*ccr(i,1).*(1 + ccr(i,2).*eta);
    dNr(2:3:r2+1,i) = (1/4).*(1 - zet.^2).*(1 + ccr(i,1).*xsi).*ccr(i,2);
    dNr(3:3:r2+2,i) = (1/4).*-2.*zet.*(1 + ccr(i,1).*xsi).*(1 + ccr(i,2).*eta);
end

%% Initialisation
et = [];

%% Jacobien
JT=dNr*[ex(1,:);ey(1,:);ez(1,:)]';

%% Loop over the integration points
for i=1:ngp
    
    indx=[ 3*i-2; 3*i-1; 3*i ];
    
    detJ = det(JT(indx,:));
    
    if detJ < 10*eps
        
        error('Jacobian determinant equal or less than zero')
        
    end
    
    dNx=JT(indx,:)\dNr(indx,:);
    
    B(1,1:3:60-2)=dNx(1,:);
    B(2,2:3:60-1)=dNx(2,:);
    B(3,3:3:60)  =dNx(3,:);
    B(4,1:3:60-2)=dNx(2,:);
    B(4,2:3:60-1)=dNx(1,:);
    B(5,1:3:60-2)=dNx(3,:);
    B(5,3:3:60)  =dNx(1,:);
    B(6,2:3:60-1)=dNx(3,:);
    B(6,3:3:60)  =dNx(2,:);
    
    ee=B*ed';
    
    et=[et; ee'];
    
end

%% Stockage et sortie
strain = et;