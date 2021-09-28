%% ASEN 3111 - Computational Assignment 01 - Main
% The goal of this script is to demonstrate the use of numerical
% integration techniques to calculate aerodynamic forces and
% coefficients. The main file uses two functions: Trapezoidal_Integration.m
% and simpsons_integration.m
%
%       Author: Andrew Thompson
%       Collaborators: Jacob Pendergast
%       Created: 09/04/2020 Edited: 09/17/2020

% Clear cmd, workspace, figures
clc; clear all; close all; tic

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Problem 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------------
% Calculation
%-------------------------------------------------------------------------------
% An anonymous function representing the integrand of the cl equation
clIntegrand = @(theta) 0.5*(4*sin(theta)^2 + 4*sin(theta))*sin(theta);
% An anonymous function representing the integrand of the cd equation 
cdIntegrand = @(theta) 0.5*(4*sin(theta)^2*cos(theta) + sin(theta)*cos(theta));
% Array of different numbers of panels
numberOfPanelsArray = 1:5;
% Preallocate Arrays
clTraps=zeros(1,length(numberOfPanelsArray)); % c_L Trapezoidal
clSimps=zeros(1,length(numberOfPanelsArray)); % c_L Simpson's
cdTraps=zeros(1,length(numberOfPanelsArray)); % c_D Trapezoidal
cdSimps=zeros(1,length(numberOfPanelsArray)); % c_D Simpson's
% For every panel number to test 
for i = numberOfPanelsArray
   % Trapezoidally integrate for CL
   clTraps(i) = Trapezoidal_Integration(clIntegrand, 0, 2*pi, i);
   % Trapezoidally integrate for CD
   cdTraps(i) = Trapezoidal_Integration(cdIntegrand, 0, 2*pi, i);
   % Simpson's rule integration for Cl 
   clSimps(i) = simpsons_integration(clIntegrand, 0, 2*pi, i);
   % Simposn's rule integration for CD
   cdSimps(i) = simpsons_integration(cdIntegrand, 0, 2*pi, i);
end
fprintf("\n---------------------------PROBLEM 1---------------------------\n")
fprintf("\nAnalytically determined lift and drag coefficients:")
fprintf("\n\n\tc_l = 2*pi\t\tc_d = 0\n\n")
fprintf("Panels for Trapezoidal Rule 0.1%% Relative Error: \tN=3\n")
fprintf("Panels for Simpson's Rule 0.1%% Relative Error:   \tN=3\n\n")

%-------------------------------------------------------------------------------
% Plotting
%-------------------------------------------------------------------------------
% Figure for Cl vs N
figure; hold on;
% Plot trapezoidal and simpson's rule
plot(numberOfPanelsArray, clTraps, 'LineWidth', 2)
plot(numberOfPanelsArray, cdTraps, 'LineWidth', 2)
% Box, grid, labels, title, legend
box on; grid minor;
xlabel('N'); ylabel('c_L'); title("Coefficient of Lift and Drag using Trapezoidal Rule")
legend("c_l", "c_d")
% Figure for Cd vs N
figure; hold on;
% Plot trapezoidal and simpson's rule 
plot(numberOfPanelsArray, clSimps, 'LineWidth', 2) 
plot(numberOfPanelsArray, cdSimps, 'LineWidth', 2)
% Box, grid, labels, title, legend
box on; grid minor;
xlabel('N'); ylabel('c_d'); title("Coefficient of Lift and Drag using Simpson's Rule")
legend("c_l", "c_d")


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Problem 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------------
% Constants & Data
%-------------------------------------------------------------------------------
airSpeedInf = 30;          % m/s
densityInf = 1.225;        % kg/m^3
pressureInf = 101.3*10^3;  % Pa
t = 0.12;                  % m
c = 2;                     % m
alpha = 9;                 % degrees

% Load cp Data 
cpStruct = load('Cp.mat');
% Renaming cp variables for camelCase
cpUpperStruct = cpStruct.Cp_upper;
cpLowerStruct = cpStruct.Cp_lower;
% Dynamic Pressure of the flow
dynamicPressure = 0.5*densityInf*airSpeedInf^2;

%-------------------------------------------------------------------------------
% Integrands
%-------------------------------------------------------------------------------
% Anonymous fucntion to represent the shape of a NACA 0012
ytUp = @(x) (c*t/0.2)*(0.2969*sqrt(x/c)-0.126*(x/c)-0.3516*(x/c)^2+0.2843*(x/c)^3-0.1036*(x/c)^4);
% The lower side is the same but negative
ytLo = @(x) -1 * ytUp(x);
% Anonymous function to represent local pressure on the upper surface
pressureUp = @(x) dynamicPressure*fnval(cpUpperStruct, x/c) + pressureInf;
% Anonymous function to represent lcoal pressure on the lower surface
pressureLo = @(x) dynamicPressure*fnval(cpLowerStruct, x/c) + pressureInf;
% Anonymous function to represent the integrand of Normal force
normalIntegrand = @(x) (pressureLo(x) - pressureUp(x));

%-------------------------------------------------------------------------------
% Lift and Drag Computation
%-------------------------------------------------------------------------------
N=450;
% Integrating to find normal force per unit span 
normalUnitSpan = Trapezoidal_Integration(normalIntegrand, 0, 2, N);
% Integrating to find the upper axial force per unit span 
axialLo = Trapezoidal_Integration(pressureLo, 0, 2, N, ytLo);
% Integrating to find the lower axial force per unit span 
axialUp = Trapezoidal_Integration(pressureUp, 0, 2, N, ytUp);
% Take the difference to find axial force
axialUnitSpan = axialLo - axialUp;
% Lift equation
liftUnitSpan = normalUnitSpan * cosd(alpha) - axialUnitSpan * sind(alpha);
% Drag equation
dragUnitSpan = normalUnitSpan * sind(alpha) + axialUnitSpan * cosd(alpha);

%-------------------------------------------------------------------------------
% Relative Error Computation
%-------------------------------------------------------------------------------
% Initializing Status variables for while loop 
error = 1;
panels = 20;
fiveFound = false; oneFound = false; tenthFound = false;
% While relative error is greater than 1/10th of a percent
while error > 0.001
    % Calculate the normal force per unit span 
    normal = Trapezoidal_Integration(normalIntegrand, 0, 2, panels); 
    % Integrating to find the upper axial force per unit span 
    axialLo = Trapezoidal_Integration(pressureLo, 0, 2, N, ytLo);
    % Integrating to find the lower axial force per unit span 
    axialUp = Trapezoidal_Integration(pressureUp, 0, 2, N, ytUp);
    % Take the difference to find the axial force
    axial = axialLo - axialUp;
    % Calculate the lift force per unit span 
    lift = normal * cosd(alpha) - axial * sind(alpha);
    % Calculate the relative error 
    error = abs(liftUnitSpan - lift)/liftUnitSpan;
    % If the error is less than 5%
    if error < 0.05 && fiveFound == false
        % Document the number of panels 
        panelsFivePercentE = panels;
        fiveFound = true;
        panels = panels*2;
    % If the error is less than 1%
    elseif error < 0.01 && oneFound == false
        % Document the number of panels 
        panelsOnePercentE = panels;
        oneFound = true;
        panels = panels*2;
    % If the error is less than 1/10th percent 
    elseif error < 0.001 && tenthFound == false
        % Document the number of panels 
        panelsTenthPercentE = panels;
        tenthFound = true;
    % Otherwise,
    else
        % Increment the number of panels
        panels = panels + 1;
    end
end
fprintf("---------------------------PROBLEM 2---------------------------\n")
fprintf("Calculated lift and drag per unit span using Trapezoidal Rule:\n\n")
fprintf("\tN' = %i\t\tD' = %i\n\n",liftUnitSpan,dragUnitSpan)
fprintf("Panels for Trapezoidal Rule 5%% Relative Error:\t\tN=%i\n",panelsFivePercentE)
fprintf("Panels for Trapezoidal Rule 1%% Relative Error:\t\tN=%i\n",panelsOnePercentE)
fprintf("Panels for Trapezoidal Rule 0.1%% Relative Error:\tN=%i\n\n",panelsTenthPercentE)
toc