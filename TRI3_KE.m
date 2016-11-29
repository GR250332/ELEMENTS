function Ke = TRI3_KE(ex,ey,ma,type_module,varargin)
% function Ke = TRI3_KE(ex,ey,ma,type_module,varargin)
%--------------------------------------------------------------------------
% PURPOSE
%    Calculate the elementary stiffness matrix for a TRI3 element
%--------------------------------------------------------------------------
% INPUT   
%    ex          : Elementary nodal X coordinates
%    ey          : Elementary nodal Y coordinates
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
% COMMENTS 
%     - 01/03/2016 - BR : Replacement of field.varf by field.varf(:,1) due
%     to the fact 2 internal variables are considered (D and eps_eq) for
%     Mazars' model
%
%--------------------------------------------------------------------------

%% Ouverture variables globales
global options

%% Calcul du module selon type_module
switch type_module
    
    case 'INITIAL'
        
        D = ghooke(ma);
        
    case 'SECANT'
              
        % Looking for the fields
        field = varargin{1};

        
        D = ghooke_secant(ma,field);
        
        D = cell2mat(D);
        
    case 'TANGENT'
        
        % Looking for the fields
        field = varargin{1};
                
        D = ghooke_tangent(ma,field);
        
        D = cell2mat(D);
        
    otherwise
        
        error('Case not implemented')
        
end

%% Quelques informations récupérées
t = ma.cara.epai;

%% Integration analytique
C = [ 1  ex(1) ey(1)   0     0       0
    0    0     0     1   ex(1)   ey(1)
    1  ex(2) ey(2)   0     0       0
    0    0     0     1   ex(2)   ey(2)
    1  ex(3) ey(3)   0     0       0
    0    0     0     1   ex(3)   ey(3)];

A = 1/2*abs(det([ones(3,1) ex' ey']));

%% Switch selon le mode de calcul
switch options.mode
    
    case 'PLAN_CONT'
        
       B=[0 1 0 0 0 0
          0 0 0 0 0 1
          0 0 1 0 1 0]/C;
              
       Dm = D;

       Ke = sparse(B'*Dm*B*A*t);

    case 'PLAN_DEFO'
        
       B=[0 1 0 0 0 0
           0 0 0 0 0 1
           0 0 1 0 1 0]/C;
       
       Dm = D([1 2 4],[1 2 4]);
       
       Ke = sparse(B'*Dm*B*A*t);
       
    otherwise
        
        error('Case not implemented')
   
end