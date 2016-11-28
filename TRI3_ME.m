function Me = TRI3_ME(ex,ey,ma,varargin)
% function Me = TRI3_ME(ex,ey,ma,varargin)
%--------------------------------------------------------------------------
% PURPOSE
%    Calculate the elementary mass matrix for a TRI3 element
%--------------------------------------------------------------------------
% INPUT   
%    ex          : Elementary nodal X coordinates
%    ey          : Elementary nodal Y coordinates
%    ma          : Matarial structure
%    varargin    : element fields
%--------------------------------------------------------------------------
% OUTPUT
%    Me          : local mass matrix
%--------------------------------------------------------------------------
% REFERENCES
%     Benjamin RICHARD
%     11-02-2016
%--------------------------------------------------------------------------

%% Quelques informations récupérées
t = ma.cara.epai;
rho = ma.cara.rho;

%% Jacobien
A = 1/2*det([ones(3,1) ex' ey']);

%%
Me = A * rho * t / 12 * ...
    [2 0 1 0 1 0 ; ...
     0 2 0 1 0 1 ; ...
     1 0 2 0 1 0 ; ...
     0 1 0 2 0 1 ; ...
     1 0 1 0 2 0 ; ...
     0 1 0 1 0 2 ];
