%% ASEN 3111 Aerodynamics Computational Assignment #4 - ComputeLiftProperties.m
% Computes the zero lift angle of attack and the lift slope curve for an
% airfoil using the vortex panel method.
%
%   Author: Andrew Thompson
%   Created: 11/6/2020 Edited: 11/19/2020
%
%   Parameters:     xb [Nx1] <double>
%                   yb [Nx1] <double>
%                   vInf <double>
%   Returns:        alpha0 <double>
%                   a <double>
%
function [alpha0, a] = ComputeLiftProperties(xb, yb, vInf)

    alpha = linspace(-5, 5, 300);
    for i = 1:300
        cl(i) = Vortex_Panel(xb, yb, vInf, alpha(i), false);
    end
    % Perform linear least squares fit
    coeffs = polyfit(alpha, cl, 1);
    % Get the slope 
    m = coeffs(1);
    % Get the y intercept 
    b = coeffs(2);
    % Solve for zero lift angle of attack
    alpha0 = -b/m;
    % Lift curve slope
    a = m;
end