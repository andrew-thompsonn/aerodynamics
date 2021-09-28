%% ASEN 3111 - Computational Assignment 01 - Trapezoidal_Integration.m
%
%       Author: Andrew Thompson
%       Created: 09/04/2020 Edited: 09/17/2020

%%
function output = Trapezoidal_Integration(integrand, a, b, numberOfPanels, stepFunc)
    %---------------------------------------------------------------------------
    % A function to approximate the area under a curve using the
    % trapezoidal rule. Takes a start point, end point, and a number of
    % panels to use when approximating the integral.
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
    % Linearly spaced array of steps
    stepArray = a:stepSize:b;
    % Initialize a variable to store the solution
    areaSum = 0;
    % For every step in the array excluding the last step
    for i = stepArray(1:numberOfPanels)
        % Get left hand values of the integrand
        leftFunctionValue = integrand(i);
        % Get right hand values of the integrand
        rightFunctionValue = integrand(i+stepSize);
        % Average the two values
        averageFunctionValue = (leftFunctionValue + rightFunctionValue)/2;
        % If only four inputs 
        if nargin == 4
            % Multiply the average value of the integrand by the step size
            area = stepSize * averageFunctionValue;
        % Otherwise, 
        else 
            % Multiply the average value of the integrand by change in y
            area = (stepFunc(i+stepSize) - (stepFunc(i)))*averageFunctionValue;
        end 
        % Increment the sum with the calculated area
        areaSum = areaSum + area;
    end
    % Return the total area under the curve
    output = areaSum;
end