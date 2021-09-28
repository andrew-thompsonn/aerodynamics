%% ASEN 3111 Aerodynamics Computational Assignment #4 - Main
% This script uses Prandtl lifting line theory to estimate the lift and
% drag on a finite wing, performs an error analysis on the estimate
% obtained using theory, and analyzes the effect taper ratio has on
% efficiency at varying aspect ratios. 
%
%   Author: Andrew Thompson
%   Created: 10/30/20 Edited: 11/13/20

%% Housekeeping
clc; clear all; close all; tic

%% Flow Parameters
vInf = 150*1.46667; % fps
rho = 23.77e-4;     % slug/ft^3

%% Wing Parameters
% Span
b = 100; 
% Chord length
cR = 15;  cT = 5;             
% Geometric AoA
geoR = 5; geoT = 0;           
% Planform Area
S = b*(cR + cT)/2;  

% Zero lift AoA & Lift slope NACA 2412
[xbR, ybR] = NACA_Airfoils(0.02, 0.4, 0.12, cR, 300);
[aeroR, a0R] = ComputeLiftProperties(xbR, ybR, vInf);
% Zero lift AoA & Lift slope NACA 0012
[xbT, ybT] = NACA_Airfoils(0, 0, 0.12, cT, 300);
[aeroT, a0T] = ComputeLiftProperties(xbT, ybT, vInf);

%% Parameter Struct
p = struct();
% So function inputs aren't disgusting
p.vInf  = vInf;     p.rho   = rho;
p.b     = 100;      p.S     = S;
p.cR    = cR;       p.cT    = cT;
p.geoR  = geoR;     p.geoT  = geoT;
p.aeroR = aeroR;    p.aeroT = aeroT; 
p.a0R   = a0R;      p.a0T   = a0T;  

%% Calculations
% Compute coefficients, lift force, drag force
[e, cL, cDi] = PLLT(b, a0T, a0R, cT, cR, aeroT, aeroR, geoT, geoR, 1000);
liftRef = cL *0.5*rho*vInf^2*S;
dragRef = cDi*0.5*rho*vInf^2*S;
% Error analysis 
[N5, N1, N01] = RelativeError(liftRef, dragRef, p);
% Taper analysis 
EfficiencyAnalysis(N01);

%% Results
fprintf("\n************* REFERENCE VALUES FROM 1000 ODD TERMS ************\n\n")
fprintf("\tCL: %.3f\n\tCDi: %.3f\n\te: %.3f\n\n",cL, cDi, e)
fprintf("\tL: %.3f [kip]\n\tDi: %.3f [kip]\n",liftRef/1000, dragRef/1000)
fprintf("\n******** NUMBER OF ODD TERMS FOR 5,1,0.1 PERCENT ERROR ********\n\n")
fprintf("\t5 percent: %i\n\t1 percent: %i\n\t1/10 percent: %i\n\n",N5,N1,N01)
fprintf("***************************************************************\n");toc