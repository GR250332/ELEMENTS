function strain = BAR2_DE(ex,ey,ed)
% function strain = BAR2_DE(ex,ey,ed)
%--------------------------------------------------------------------------
% PURPOSE
%    Computation of the strain from a displacement vector ed for a BARRE
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

%% Computation of the length of the element
b = [ex(2) - ex(1) ; ey(2) - ey(1)];

L = sqrt(b'*b);

%% Computation of the expansion matrix
n = b'/L;

G = [n zeros(size(n)) ; ...
    zeros(size(n)) n];

%% Computation of the strains
strain = [-1 1] * G * ed' / L;
