%% ASEN 3111 Aerodynamics Computational Assignment 2 - main
% This script calls the function Plot_Airfoil_Flow and Plot_RSME to
% produce the requested figures in the Computational Assignment #2 document
%
%       Author: Andrew Thompson
%       Created: 9/23/2020          Edited: 10/2/2020
%

clc; clear; close all;

%% Constants
rhoInf = 1.225;             % kg/m^3
pressureInf = 101.3*10^3;   % Pa
velocityInf = 68;           % m/s
alpha = 12;                 % deg
chord = 2;                  % m

% Parameter struct to be passed into Plot_RSME() and Vary_Parameters()
params = struct();
params.rhoInf = rhoInf;
params.pressureInf = pressureInf;
params.velocityInf = velocityInf;
params.alpha = alpha;
params.chord = chord;  

%% Problem 1 - Calculate and plot psi, phi, p
[pRef,psiRef,phiRef]=Plot_Airfoil_Flow(chord,alpha,velocityInf,pressureInf, ...
    rhoInf,350,true);
vRef = zeros(100,100);

%vRef = sqrt((2/rhoInf)*(pressureInf - pRef));
%% Problem 2 - Calculate and plot RSME of phi, psi, p for varying N
Plot_RSME(pRef, psiRef,phiRef,vRef, params);

%% Problem 3 - Characterize changes in psi and phi for varying parameters
Vary_Parameters(linspace(1,4,4),linspace(4,13,4),linspace(30,90,4),params);