%% ASEN 3111 - Computational Assignment 01 - simpsons_integration.m
%
%       Author: Andrew Thompson
%       Created: 09/04/2020 Edited: 09/17/2020

%%
function output = simpsons_integration(integrand, a, b, numberOfPanels)
    %---------------------------------------------------------------------------
    % A function to approximate the area under a curve using simpsons's rule. 
    % Takes a start point, end point, and a number of panels to use when 
    % approximating the integral.
    %
    %   Inputs:     integrand <function_handle>
    %               a <double>
    %               b <double>
    %               numberOfPanels <int>
    %
    %   Outputs:    areaSum <double>
    %---------------------------------------------------------------------------
    
    % Calculate the stepsize
    stepSize = (b - a) / numberOfPanels;
    % Create a step array 
    stepArray = a:stepSize:b;
    % Initialize sum of areas
    areaSum = 0;
    % For all steps in the array excluding the last step
    for i = stepArray(1:numberOfPanels)
        % Get the function value between the two steps
        functionValue = integrand(i + (stepSize/2));
        % Calculate area of panel 
        area = stepSize * functionValue;
        % Add to area sum 
        areaSum = areaSum + area;
    end
    % Output total area 
    output = areaSum;
end