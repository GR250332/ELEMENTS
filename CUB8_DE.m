function strain = CUB8_DE(ex,ey,ez,ed)
% function strain = CUB8_DE(ex,ey,ez,ed)
%--------------------------------------------------------------------------
% PURPOSE
%    Computation of the strain from a displacement vector ed for a CUB8
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
%     Benjamin RICHARD
%     11-02-2016
%--------------------------------------------------------------------------

%% Position of the integration points
ir  = 2;
ngp = ir*ir*ir;

if ir==1
    
    error('Case not implemented')
    
elseif ir==2
    
    g1=0.577350269189626;
    
    gp(:,1)=[-1; 1; 1;-1;-1; 1; 1;-1]*g1;
    gp(:,2)=[-1;-1; 1; 1;-1;-1; 1; 1]*g1;
    gp(:,3)=[-1;-1;-1;-1; 1; 1; 1; 1]*g1;
    
elseif ir==3
    
    g1=0.774596699241483;
    g2=0.;
    
    I1=[-1; 0; 1;-1; 0; 1;-1; 0; 1]';
    I2=[ 0;-1; 0; 0; 1; 0; 0; 1; 0]';
    
    gp(:,1)=[I1 I1 I1]'*g1;
    gp(:,1)=[I2 I2 I2]'*g2+gp(:,1);
    
    I1=[-1;-1;-1; 0; 0; 0; 1; 1; 1]';
    I2=[ 0; 0; 0; 1; 1; 1; 0; 0; 0]';
    
    gp(:,2)=[I1 I1 I1]'*g1;
    gp(:,2)=[I2 I2 I2]'*g2+gp(:,2);
    
    I1=[-1;-1;-1;-1;-1;-1;-1;-1;-1]';
    I2=[ 0; 0; 0; 0; 0; 0; 0; 0; 0]';
    I3=abs(I1);
    
    gp(:,3)=[I1 I2 I3]'*g1;
    gp(:,3)=[I2 I3 I2]'*g2+gp(:,3);
    
else
    
    error('Case not implemented')
    
end

%% Coordinates of the integration points in the reference elements
xsi = gp(:,1);
eta = gp(:,2);
zet = gp(:,3);
r2  = ngp*3;

%% Gradient of the shape functions
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
    
    B(1,1:3:24-2)=dNx(1,:);
    B(2,2:3:24-1)=dNx(2,:);
    B(3,3:3:24)  =dNx(3,:);
    B(4,1:3:24-2)=dNx(2,:);
    B(4,2:3:24-1)=dNx(1,:);
    B(5,1:3:24-2)=dNx(3,:);
    B(5,3:3:24)  =dNx(1,:);
    B(6,2:3:24-1)=dNx(3,:);
    B(6,3:3:24)  =dNx(2,:);
    
    ee=B*ed';
    
    et=[et; ee'];
    
end

%% Stockage et sortie
strain = et;