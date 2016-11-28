function Ke = BAR2_KE(ex,ey,ma,type_module,varargin)
% function Ke = BAR2_KE(ex,ey,ma,type_module,varargin)
%--------------------------------------------------------------------------
% PURPOSE
%    Calculate the elementary stiffness matrix for a BAR2 element
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
%     25-02-2016
%--------------------------------------------------------------------------

%% Calcul du module selon type_module
switch type_module
    
    case 'INITIAL'
        
        E = ma.cara.youn;
        A = ma.cara.sect;
        
    case 'SECANT'
        
        % Looking for the fields
        field = varargin{1};
        
        % Safety check
        if size(field.var0,2) > 1
            
            % Looking for the fields
            field = varargin{1};
            
            % Computation of the elastic stiffness
            E = ma.cara.youn;
            
            % Computation of the secant modulus
            E = (1 - field.varf(:,1)) * E;
            
        else
            
            E = ma.cara.youn;
            
        end
        
        A = ma.cara.sect;
        
    case 'TANGENT'
        
        error('Case not implemented')
        % Il faudra donner un champ par element de variables internes a
        % minima et eventuellement contraintes + deformations
        % D = ghooke_tangent(ma,varargin{1});
        
    otherwise
        
        error('Case not implemented')
        
end

%% Node coordinates
b = [ ex(2)-ex(1); ey(2)-ey(1) ];

%% Length
L = sqrt(b'*b);

%% Local matrix
Kle = E*A/L*[ 1 -1;
    -1  1];

n = b'/L;

%% Expansion
G = [   n      zeros(size(n));
    zeros(size(n))     n   ];

Ke = G' * Kle * G;
