%% ASEN 3111 Computational Assignment 3 - Zero_Lift_Alpha.m
% The purpose of this function is to compute the zero lfit angle of attack
% for an airfoil by using a least squares fit on a cl vs alpha curve.
%
%   Author: Andrew Thompson
%   Created: 10/27/2020 Edited: 10/27/2020
%
%   Parameters:     cl [1 x N] <double> - Coefficients of lift 
%                   alpha [1 x N] <double> - Angles of attack
%
%   Outputs:        alpha0 <Double> - Zero lift angle of attack

function alpha0 = Zero_Lift_Alpha(alpha, cl)
    % Perform linear least squares fit
    coeffs = polyfit(alpha, cl, 1);
    % Get the slope 
    m = coeffs(1);
    % Get the y intercept 
    b = coeffs(2);
    % Solve for zero lift angle of attack
    alpha0 = -b/m;
end