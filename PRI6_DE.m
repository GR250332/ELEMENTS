function strain = PRI6_DE(ex,ey,ez,ed)
% function strain = PRI6_DE(ex,ey,ez,ed)
%--------------------------------------------------------------------------
% PURPOSE
%    Computation of the strain from a displacement vector ed for a PRI6
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
ngp = 6;

if ngp == 6
    
    g1 = 1/sqrt(3);
    g2 = 1/2;
    g3 = 0;
    
    gp(:,1)=[-g1; -g1; -g1; g1; g1; g1];
    gp(:,2)=[g2; g3; g2; g2; g3; g2];
    gp(:,3)=[g2; g2; g3; g2; g2; g3];

elseif ngp == 8
    
    g1 = 0.577350269189626;
    g2 = 1/3;
    g3 = 0.6;
    g4 = 0.2;
    
    gp(:,1)=[-g1; -g1; -g1; -g1; g1; g1; g1; g1;];
    gp(:,2)=[g2; g3; g4; g4; g2; g3; g4; g4];
    gp(:,3)=[g2; g4; g3; g4; g2; g4; g3; g4];
    
else
    
    error('Used number of integration points not implemented');
    
end;

%% Coordinates of the integration points in the reference elements
xsi = gp(:,1);
eta = gp(:,2);
zet = gp(:,3);
r2  = ngp*3;

%% Gradient of the shape functions
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
    
    B(1,1:3:18-2)=dNx(1,:);
    B(2,2:3:18-1)=dNx(2,:);
    B(3,3:3:18)  =dNx(3,:);
    B(4,1:3:18-2)=dNx(2,:);
    B(4,2:3:18-1)=dNx(1,:);
    B(5,1:3:18-2)=dNx(3,:);
    B(5,3:3:18)  =dNx(1,:);
    B(6,2:3:18-1)=dNx(3,:);
    B(6,3:3:18)  =dNx(2,:);
    
    ee=B*ed';
    
    et=[et; ee'];
    
end

%% Stockage et sortie
strain = et;