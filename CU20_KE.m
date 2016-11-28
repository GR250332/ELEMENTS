function Ke = CU20_KE(ex,ey,ez,ma,type_module,varargin)
%  function Ke = CU20_KE(ex,ey,ez,ma,type_module,varargin)
%--------------------------------------------------------------------------
% PURPOSE
%    Calculate the elementary stiffness matrix for a CU20 element
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
%     Romili PAREDES
%     25-07-2016
%--------------------------------------------------------------------------
%% Ouverture variables globales
global options

%% Integration point number
ngp = 27;

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
if ngp == 8
    
    g1 = 1/sqrt(3);
    
    w1 = 1;
    
    gp(:,1)=[-g1; -g1; -g1; -g1; g1; g1; g1; g1];
    gp(:,2)=[-g1; -g1; g1; g1; -g1; -g1; g1; g1];
    gp(:,3)=[-g1; g1; -g1; g1; -g1; g1; -g1; g1];
    
    wp = [w1; w1; w1; w1; w1; w1; w1; w1];
    
elseif ngp == 27
    
    g1 = sqrt(3/5);
    g2 = 0;
    
    w1 = 5/9;
    w2 = 8/9;
    
    gp(:,1)=[-g1; -g1; -g1; -g1; -g1; -g1; -g1; -g1; -g1; ...
             g2; g2; g2; g2; g2; g2; g2; g2; g2; ...
             g1; g1; g1; g1; g1; g1; g1; g1; g1];
    gp(:,2)=[-g1; -g1; -g1; g2; g2; g2; g1; g1; g1; ...
             -g1; -g1; -g1; g2; g2; g2; g1; g1; g1; ...
             -g1; -g1; -g1; g2; g2; g2; g1; g1; g1];
    gp(:,3)=[-g1; g2; g1; -g1; g2; g1; -g1; g2; g1; ...
             -g1; g2; g1; -g1; g2; g1; -g1; g2; g1; ...
             -g1; g2; g1; -g1; g2; g1; -g1; g2; g1];
    
    wp = [w1^3; w1^2*w2; w1^3; w1^2*w2; w2^2*w1; w1^2*w2; w1^3; w1^2*w2; w1^3; ...
          w1^2*w2; w2^2*w1; w1^2*w2; w2^2*w1; w2^3; w2^2*w1; w1^2*w2; w2^2*w1; w1^2*w2; ...
          w1^3; w1^2*w2; w1^3; w1^2*w2; w2^2*w1; w1^2*w2; w1^3; w1^2*w2; w1^3];
    
else
    
    error('Used number of integration points not implemented');
    
end;

%% Precalcule de quantités pour l'intégration

% Coordonnées locales dans le repère de référénce
xsi = gp(:,1);
eta = gp(:,2);
zet = gp(:,3);

% Pour la stockage des valeurs des fonctions de formes et de leurs
% gradients
r2 = ngp*3;

%% Fonctions de formes et gradients

% Coordonnées dans l'espace de réference

ccr(:,1) = [-1; -1; -1; -1; -1; -1; -1; -1; 0; 0; 0; 0; 1; 1; 1; 1; 1; 1; 1; 1];
ccr(:,2) = [1; 0; -1; -1; -1; 0; 1; 1; 1; -1; -1; 1; 1; 0; -1; -1; -1; 0; 1; 1];
ccr(:,3) = [1; 1; 1; 0; -1; -1; -1; 0; 1; 1; -1; -1; 1; 1; 1; 0; -1; -1; -1; 0];

% Fonction de forme

for i = [5 17 19 7 3 15 13 1]
	N(:,i) = (1/8).*(1 + ccr(i,1).*xsi).*(1 + ccr(i,2).*eta).*(1 + ccr(i,3).*zet).*(ccr(i,1).*xsi + ccr(i,2).*eta + ccr(i,3).*zet - 2);
end

for i = [11 12 10 9]
    N(:,i) = (1/4).*(1 - xsi.^2).*(1 + ccr(i,2).*eta).*(1 + ccr(i,3).*zet);
end

for i = [18 6 14 2]
    N(:,i) = (1/4).*(1 - eta.^2).*(1 + ccr(i,1).*xsi).*(1 + ccr(i,3).*zet);
end

for i = [4 16 20 8]
    N(:,i) = (1/4).*(1 - zet.^2).*(1 + ccr(i,1).*xsi).*(1 + ccr(i,2).*eta);
end

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

%% Initialisation avant integration numérique

% Matrice elementaire
Ke=sparse(60,60);

% Jacobian
JT=dNr*[ex;ey;ez]';

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
        
        B(1,1:3:60-2) = dNx(1,:);
        B(2,2:3:60-1) = dNx(2,:);
        B(3,3:3:60)   = dNx(3,:);
        B(4,1:3:60-2) = dNx(2,:);
        B(4,2:3:60-1) = dNx(1,:);
        B(5,1:3:60-2) = dNx(3,:);
        B(5,3:3:60)   = dNx(1,:);
        B(6,2:3:60-1) = dNx(3,:);
        B(6,3:3:60)   = dNx(2,:);
        
        N2(1,1:3:60-2) = N(i,:);
        N2(2,2:3:60-1) = N(i,:);
        N2(3,3:3:60)   = N(i,:);
        
        Ke=Ke+B'*D{i}*B*detJ*wp(i);
        
    end
end
