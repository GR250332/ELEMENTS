function Ke = QUA4_KE(ex,ey,ma,type_module,varargin)
% function Ke = QUA4_KE(ex,ey,ma,type_module,varargin)
%--------------------------------------------------------------------------
% PURPOSE
%    Calculate the elementary stiffness matrix for a QUA4 element
%--------------------------------------------------------------------------
% INPUT   
%    Ex          : Elementary nodal X coordinates
%    Ey          : Elementary nodal Y coordinates
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

%% Number of integration points
ngp = 4;

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
        
        % Looking for the fields
        field = varargin{1};
                
        D = ghooke_tangent(ma,field);
        
    otherwise
        
        error('Case not implemented')
        
end

%% Regles d'integration numerique

% Points d'intégration et poids
if ngp==1
    
    g1 = 0.0;
    
    w1 = 2.0;
    
    gp = [ g1 g1 ];
    
    w  = [ w1 w1 ];
    
elseif ngp==4
    
    g1 = 0.577350269189626;
    w1 = 1;
    
    gp(:,1) = [-g1; g1;-g1; g1];
    gp(:,2) = [-g1;-g1; g1; g1];
    
    w(:,1) = [ w1; w1; w1; w1];
    w(:,2) = [ w1; w1; w1; w1];
    
elseif ngp==9
    
    g1 = 0.774596699241483;
    g2 = 0.;
    
    w1 = 0.555555555555555;
    w2 = 0.888888888888888;
    
    gp(:,1) = [-g1;-g2; g1;-g1; g2; g1;-g1; g2; g1];
    gp(:,2) = [-g1;-g1;-g1; g2; g2; g2; g1; g1; g1];
    
    w(:,1) = [ w1; w2; w1; w1; w2; w1; w1; w2; w1];
    w(:,2) = [ w1; w1; w1; w2; w2; w2; w1; w1; w1];
    
else
    
    error('Used number of integration points not implemented');
    
end

%% Precalcule de quantités pour l'intégration

% Matrice des poids
wp  = w(:,1).*w(:,2);

% Coordonnées locales dans le repère de référénce
xsi = gp(:,1);
eta = gp(:,2);

% Pour la stockage des valeurs des fonctions de formes et de leurs
% gradients
r2  = ngp * 2;

%% Fonctions de formes et gradients

% Fonction de forme
N(:,1)=(1-xsi).*(1-eta)/4;
N(:,2)=(1+xsi).*(1-eta)/4;
N(:,3)=(1+xsi).*(1+eta)/4;
N(:,4)=(1-xsi).*(1+eta)/4;

% Derivées des fonctions de formes
dNr(1:2:r2,1) = -(1-eta)/4;
dNr(1:2:r2,2) =  (1-eta)/4;
dNr(1:2:r2,3) =  (1+eta)/4;
dNr(1:2:r2,4) = -(1+eta)/4;

dNr(2:2:r2+1,1) = -(1-xsi)/4;
dNr(2:2:r2+1,2) = -(1+xsi)/4;
dNr(2:2:r2+1,3) = (1+xsi)/4;
dNr(2:2:r2+1,4) = (1-xsi)/4;

%% Initialisation avant integration numérique

% Matrice elementaire
Ke = sparse(8,8);

% Jacobian
JT = dNr*[ex;ey]';

% Epaisseur
t  = ma.cara.epai;

%% Switch selon le cas
switch options.mode
    
    case 'PLAN_CONT'
            
                    
        % Boucle sur les points de Gauss
        for i=1:ngp
            
            % On prend la matrice réduite
            Dm = D{i};

            indx=[2*i-1; 2*i ];
            
            detJ=det(JT(indx,:));
            
            if detJ < 10*eps
                
                error('Jacobian determinant equal or less than zero')
                
            end
                        
            dNx   = JT(indx,:)\dNr(indx,:);
            
            B(1,1:2:8-1)  = dNx(1,:);
            B(2,2:2:8)    = dNx(2,:);
            B(3,1:2:8-1)  = dNx(2,:);
            B(3,2:2:8)    = dNx(1,:);
            
            Ke = Ke+B'*Dm*B*detJ*wp(i)*t;
            
        end
        
    case 'PLAN_DEFO'
                
        % Boucle sur les points d'intégration
        for i=1:ngp
            
            % On ne prend que la matrice réduite
            Dm = D{i};
            Dm = Dm([1 2 4],[1 2 4]);
            
            indx=[ 2*i-1; 2*i ];
            
            detJ=det(JT(indx,:));
            
            if detJ<10*eps
                
                error('Jacobian determinant equal or less than zero')
                
            end
            
            dNx=JT(indx,:)\dNr(indx,:);
            
            B(1,1:2:8-1)=dNx(1,:);
            B(2,2:2:8)  =dNx(2,:);
            B(3,1:2:8-1)=dNx(2,:);
            B(3,2:2:8)  =dNx(1,:);
            
            Ke=Ke+B'*Dm*B*detJ*wp(i)*t;
            
        end
        
    otherwise
        
        error('Case not implemented')
        
end