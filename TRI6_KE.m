function Ke = TRI6_KE(ex,ey,ma,type_module,varargin)
%  function Ke = TRI6_KE(ex,ey,ez,ma,type_module,varargin)
%--------------------------------------------------------------------------
% PURPOSE
%    Calculate the elementary stiffness matrix for a TET4 element
%--------------------------------------------------------------------------
% INPUT
%    Ex          : Elementary nodal X coordinates
%    Ey          : Elementary nodal Y coordinates
%    Ez          : Elementary nodal Z coordinates
%    ma          : Material structure
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
if ngp == 1
    
    g1  = 1/3;
    w1  = 1/2;
    
    gp  = [g1 g1];
    
    wp  = w1;
    
elseif ngp == 3
    
    g1 = 1/6;
    g2 = 2/3;
    
    w1 = 1/6;
    
    gp(:,1) = [g1; g2; g1];
    gp(:,2) = [g1; g1; g2];
    
    wp = [ w1; w1; w1];
    
elseif ngp == 4
    
    g1 = 1/5;
    g2 = 3/5;
    g3 = 1/3; 
    
    w1 = 25/(24*4);
    w2 = -27/(24*4);
    
    gp(:,1) = [g1; g2; g1; g3];
    gp(:,2) = [g1; g1; g2; g3];
    
    wp = [w1; w1; w1; w2];
    
elseif ngp == 7
    
    g1 = 1/3;
    A  = 0.470142064105115;
    B  = 0.101286507323456;
    
    w1 = 9/80;
    w2 = 0.066197076394253;
    w3 = 0.062969590272413;
    
    gp(:,1) = [g1; A; 1-2*A; A; B; 1-2*B; B];
    gp(:,2) = [g1; A; A; 1-2*A; B; B; 1-2*B];
    
    wp = [w1; w2; w2; w2; w3; w3; w3];
    
else
    
    error('Used number of integration points not implemented');
    
end;

%% Precalcule de quantités pour l'intégration

% Coordonnées locales dans le repère de référénce
xsi = gp(:,1);
eta = gp(:,2);

% Pour la stockage des valeurs des fonctions de formes et de leurs
% gradients
r2  = ngp * 2;

%% Fonctions de formes et gradients

% Fonction de forme
N(:,1) = -(1 - xsi - eta).*(1 - 2.*(1 - xsi - eta));
N(:,3) = -xsi.*(1 - 2.*xsi);
N(:,5) = -eta.*(1 - 2.*eta);
N(:,2) = 4.*xsi.*(1 - xsi - eta);
N(:,4) = 4.*xsi.*eta;
N(:,6) = 4.*eta.*(1 - xsi - eta);

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

%% Initialisation avant integration numérique

% Matrice elementaire
Ke = sparse(12, 12);

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

            indx=[2*i - 1; 2*i];
            
            detJ=det(JT(indx,:));
            
            if detJ < 10*eps
                
                error('Jacobian determinant equal or less than zero')
                
            end
                        
            dNx   = JT(indx,:)\dNr(indx,:);
            
            B(1,1:2:12-1)  = dNx(1,:);
            B(2,2:2:12)    = dNx(2,:);
            B(3,1:2:12-1)  = dNx(2,:);
            B(3,2:2:12)    = dNx(1,:);
            
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
            
            B(1,1:2:12-1)=dNx(1,:);
            B(2,2:2:12)  =dNx(2,:);
            B(3,1:2:12-1)=dNx(2,:);
            B(3,2:2:12)  =dNx(1,:);
            
            Ke=Ke+B'*Dm*B*detJ*wp(i)*t;
            
        end
        
    otherwise
        
        error('Case not implemented')
        
end



