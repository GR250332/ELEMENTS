function strain = QUA4_DE(ex,ey,ed)
% function strain = QUA4_DE(ex,ey,ed)
%--------------------------------------------------------------------------
% PURPOSE
%    Computation of the strain from a displacement vector ed for a QUA4
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
%     Benjamin RICHARD
%     11-02-2016
%--------------------------------------------------------------------------

%% Position of the integration points
ir  = 2;
ngp = ir*ir;

if ir==1
    
    g1 = 0.0;
    gp = [ g1 g1 ];
    
elseif ir==2
    
    g1 = 0.577350269189626;
    
    gp(:,1) = [-g1; g1;-g1; g1];
    gp(:,2) = [-g1;-g1; g1; g1];
    
elseif ir==3
    
    g1 = 0.774596699241483;
    g2 = 0.;
    
    gp(:,1) = [-g1;-g2; g1;-g1; g2; g1;-g1; g2; g1];
    gp(:,2) = [-g1;-g1;-g1; g2; g2; g2; g1; g1; g1];
    
else
    
    error('Case not implemented')
    
end

%% Coordinates of the integration points in the reference elements
xsi = gp(:,1);
eta = gp(:,2);
r2  = ngp*2;

%% Gradient of the shape functions
dNr(1:2:r2,1)   = -(1-eta)/4;
dNr(1:2:r2,2)   =  (1-eta)/4;
dNr(1:2:r2,3)   =  (1+eta)/4;
dNr(1:2:r2,4)   = -(1+eta)/4;
dNr(2:2:r2+1,1) = -(1-xsi)/4;
dNr(2:2:r2+1,2) = -(1+xsi)/4;
dNr(2:2:r2+1,3) =  (1+xsi)/4;
dNr(2:2:r2+1,4) =  (1-xsi)/4;

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
    
    B(1,1:2:8-1) = dNx(1,:);
    B(2,2:2:8)   = dNx(2,:);
    B(3,1:2:8-1) = dNx(2,:);
    B(3,2:2:8)   = dNx(1,:);
    
    ee = zeros(4,1);
    
    ee([1 2 4])=B*ed';
    
    et=[et; ee'];
    
end

%% Stockage et sortie
strain = et;