function Fe = TRI3_FE(ma,ex,ey,es)
%  function Fe = TRI3_FE(ma,ex,ey,es)
%--------------------------------------------------------------------------
% PURPOSE
%    Creation of an elementary internal force vector for TRI3
%--------------------------------------------------------------------------
% INPUT
%    ma       : material structutre
%    ex       : Elementary nodal X coordinates
%    ey       : Elementary nodal Y coordinates
%    es       : Stress at the integration point
%--------------------------------------------------------------------------
% OUTPUT
%    Fe       : Elementary internal force vector
%
%--------------------------------------------------------------------------
% REFERENCES
%     Benjamin RICHARD
%     10-02-2016
%--------------------------------------------------------------------------

%% Initialization of the Fe vector

% Thickness
t = ma.cara.epai;

%% On ne considère pas sigma_33 (soit nulle soit n'intervient pas dans Fe)
stress = es(:,[1 2 4]);

%% Exact integration
C = [ 1  ex(1,1) ey(1,1)   0          0          0
    0         0        0   1   ex(1,1)   ey(1,1)
    1  ex(1,2) ey(1,2)   0          0          0
    0         0        0   1   ex(1,2)   ey(1,2)
    1  ex(1,3) ey(1,3)   0          0          0
    0         0        0   1   ex(1,3)   ey(1,3)];

A=det([ones(3,1) ex(1,:)' ey(1,:)']);

B = [0 1 0 0 0 0;
    0 0 0 0 0 1;
    0 0 1 0 1 0]/C;

Fe = 0.5*A*t*B'*stress';
