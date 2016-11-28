function strain = TRI6_DE(ex,ey,ed)
% function strain = TRI6_DE(ex,ey,ed)
%--------------------------------------------------------------------------
% PURPOSE
%    Computation of the strain from a displacement vector ed for a TRI6
%--------------------------------------------------------------------------
% INPUT
%    ex : X coordinates
%    ey : Y coordinates
%    ed : Nodal displacements
%--------------------------------------------------------------------------
% OUTPUT
%    strain
%--------------------------------------------------------------------------
% REFERENCES
%     Romili PAREDES
%     25-07-2016
%--------------------------------------------------------------------------

%% Regles d'integration numerique
ngp = 4;

% Points d'intégration et poids
if ngp == 1
    
    g1  = 1/3;
    
    gp  = [g1 g1];
    
elseif ngp == 3
    
    g1 = 1/6;
    g2 = 2/3;
    
    gp(:,1) = [g1; g2; g1];
    gp(:,2) = [g1; g1; g2];
    
elseif ngp == 4
    
    g1 = 1/5;
    g2 = 3/5;
    g3 = 1/3; 
    
    gp(:,1) = [g1; g2; g1; g3];
    gp(:,2) = [g1; g1; g2; g3];
    
elseif ngp == 7
    
    g1 = 1/3;
    A  = 0.470142064105115;
    B  = 0.101286507323456;
    
    gp(:,1) = [g1; A; 1-2*A; A; B; 1-2*B; B];
    gp(:,2) = [g1; A; A; 1-2*A; B; B; 1-2*B];
      
else
    
    error('Used number of integration points not implemented');
    
end;

%% Coordinates of the integration points in the reference elements
xsi = gp(:,1);
eta = gp(:,2);
r2  = ngp*2;

% Derivées des fonctions de formes
dNr(1:2:r2,1) = 4.*xsi + 4.*eta - 3;
dNr(1:2:r2,3) = 4.*xsi - 1;
dNr(1:2:r2,5) = 0;
dNr(1:2:r2,2) = 4 - 8.*xsi - 4.*eta;
dNr(1:2:r2,4) = 4.*eta;
dNr(1:2:r2,6) = -4.*eta;

dNr(2:2:r2+1,1) = 4.*xsi + 4.*eta - 3;
dNr(2:2:r2+1,3) = 0;
dNr(2:2:r2+1,5) = 4.*eta - 1;
dNr(2:2:r2+1,2) = -4.*xsi;
dNr(2:2:r2+1,4) = 4.*xsi;
dNr(2:2:r2+1,6) = 4 - 4.*xsi - 8.*eta;

%% Initialisation
et = [];

%% Loop over the elements (if needed)

% Jacobien
JT = dNr*[ex(1,:);ey(1,:)]';

for i=1:ngp
    
    indx = [2*i-1 ; 2*i ];
    
    detJ = det(JT(indx,:));
    
    if detJ < 10*eps
        
        error('Jacobian determinant equal or less than zero')
        
    end
    
    dNx = JT(indx,:)\dNr(indx,:);
    
    B(1,1:2:12-1) = dNx(1,:);
    B(2,2:2:12)   = dNx(2,:);
    B(3,1:2:12-1) = dNx(2,:);
    B(3,2:2:12)   = dNx(1,:);
    
    ee = zeros(4,1);
    
    ee([1 2 4])=B*ed';
    
    et=[et; ee'];
    
end

%% Stockage et sortie
strain = et;