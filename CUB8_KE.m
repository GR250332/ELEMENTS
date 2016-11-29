function Ke = CUB8_KE(ex,ey,ez,ma,type_module,varargin)
%  function Ke = CUB8_KE(ex,ey,ez,ma,type_module,varargin)
%--------------------------------------------------------------------------
% PURPOSE
%    Calculate the elementary stiffness matrix for a CUB8 element
%--------------------------------------------------------------------------
% INPUT
%    Ex          : Elementary nodal X coordinates
%    Ey          : Elementary nodal Y coordinates
%    Ez          : Elementary nodal Z coordinates
%    ma          : Matarial structure
%    type_module : keyword according to the type of mudulus
%    varargin    : element fields
%--------------------------------------------------------------------------
% OUTPUT
%    Ke          : local stiffness matrix
%--------------------------------------------------------------------------
% REFERENCES
%     Benjamin RICHARD
%     11-02-2016
%--------------------------------------------------------------------------
%% Ouverture variables globales
global options

%% Integration point number
ngp = 8;

%% Calcul du module selon type_module
switch type_module
    
    case 'INITIAL'
        
        for ipint = 1:ngp
            
            D{ipint} = ghooke(ma);
        
        end
        
    case 'SECANT'
        
        % Looking for the fields
        field = varargin{1};
        
        % Computation of the elastic stiffness
        E = ghooke(ma);
        
        for ipint = 1:ngp
            
            % Computation of the secant modulus
            D{ipint} = (1 - field.varf(ipint,1)) * E;
            
        end
        
    case 'TANGENT'
        
        error('Case not implemented')
        % Il faudra donner un champ par element de variables internes a
        % minima et eventuellement contraintes + deformations
        % D = ghooke_tangent(ma,varargin);
        
    otherwise
        
        error('Case not implemented')
        
end


%% Regles d'integration numerique

% Points d'intégration et poids
if ngp==1
    
    g1=0.0; w1=2.0;
    gp=[ g1 g1 ];  w=[ w1 w1 ];
    
elseif ngp==8
    
    g1=0.577350269189626;
    
    w1=1;
    
    gp(:,1)=[-1; 1; 1;-1;-1; 1; 1;-1]*g1;
    gp(:,2)=[-1;-1; 1; 1;-1;-1; 1; 1]*g1;
    gp(:,3)=[-1;-1;-1;-1; 1; 1; 1; 1]*g1;
    
    w(:,1)=[ 1; 1; 1; 1; 1; 1; 1; 1]*w1;
    w(:,2)=[ 1; 1; 1; 1; 1; 1; 1; 1]*w1;
    w(:,3)=[ 1; 1; 1; 1; 1; 1; 1; 1]*w1;
    
elseif ngp==27
    
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
%% Initialisation avant integration numérique

% Matrice elementaire
Ke=sparse(24,24);

% Jacobian
JT=dNr*[ex;ey;ez]';

%% Mode TRID uniquement
if ~strcmp(options.mode,'TRID')
    
    error('Computation mode is not correct')
    
else
    
    % Boucle sur les points de Gauss
    for i=1:ngp
        
        indx=[ 3*i-2; 3*i-1; 3*i ];
        
        detJ=abs(det(JT(indx,:)));
        
        dNx=JT(indx,:)\dNr(indx,:);
        
        B(1,1:3:24-2) = dNx(1,:);
        B(2,2:3:24-1) = dNx(2,:);
        B(3,3:3:24)   = dNx(3,:);
        B(4,1:3:24-2) = dNx(2,:);
        B(4,2:3:24-1) = dNx(1,:);
        B(5,1:3:24-2) = dNx(3,:);
        B(5,3:3:24)   = dNx(1,:);
        B(6,2:3:24-1) = dNx(3,:);
        B(6,3:3:24)   = dNx(2,:);
        
        N2(1,1:3:24-2) = N(i,:);
        N2(2,2:3:24-1) = N(i,:);
        N2(3,3:3:24)   = N(i,:);
        
        Ke=Ke+B'*D{i}*B*detJ*wp(i);
        
    end
end
