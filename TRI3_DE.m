function strain = TRI3_DE(ex,ey,ed)
% function strain = TRI3_DE(ex,ey,ed)
%--------------------------------------------------------------------------
% PURPOSE
%    Computation of the strain from a displacement vector ed for a TRI3
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

%% Ouverture de variables
global options

%% Loop over the elements (if needed)
C = [ 1  ex(1,1) ey(1,1)   0          0          0
    0         0        0   1   ex(1,1)   ey(1,1)
    1  ex(1,2) ey(1,2)   0          0          0
    0         0        0   1   ex(1,2)   ey(1,2)
    1  ex(1,3) ey(1,3)   0          0          0
    0         0        0   1   ex(1,3)   ey(1,3)];

B = [0 1 0 0 0 0;
    0 0 0 0 0 1;
    0 0 1 0 1 0]/C;

strain_tmp = (B * ed')';

%% Switch selon le cas
switch options.mode
    
    case 'PLAN_CONT'
        
        strain(1,1) = strain_tmp(1);
        strain(1,2) = strain_tmp(2);
        strain(1,3) = 0;
        strain(1,4) = strain_tmp(3);
        
    case 'PLAN_DEFO'
        
        strain(1,1) = strain_tmp(1);
        strain(1,2) = strain_tmp(2);
        strain(1,3) = 0;
        strain(1,4) = strain_tmp(3);
        
    otherwise
        
        error('Case not implemented')
        
end