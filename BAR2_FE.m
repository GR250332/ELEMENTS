function Fe = BAR2_FE(ma,ex,ey,es)
% function Fe = BAR2_FE(ma,ex,ey,es)
%--------------------------------------------------------------------------
% PURPOSE
%    Creation of an elementary internal force vector for BAR2
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
%    Benjamin RICHARD
%    26-02-2016
%--------------------------------------------------------------------------
%% Looking for thr cross section area
A = ma.cara.sect;

%% Computation of the length of the element
b = [ex(2) - ex(1) ; ey(2) - ey(1)];

L = sqrt(b'*b);

%% Computation of the expansion matrix
n = b'/L;

G = [n zeros(size(n)) ; ...
    zeros(size(n)) n];

%% Computation of the interface force vector
Fe = A * ([-1 1] * G)' * es;
