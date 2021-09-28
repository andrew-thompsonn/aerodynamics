%% ASEN 3111 Aerodynamics Computational Assignment 2 - Plot_Airfoil_Flow.m
% The purpose of this function is to use a discrete number of vortices
% to calculate and plot the stream function, velocity potential, and pressure 
% due to flow over a thin airfoil.
%
%       Author:     Andrew Thompson
%       Created:    9/18/2020          Edited:     10/2/2020
%
%       Inputs:     chord <double>, the chord length of the airfoil
%                   alpha <double>, the angle of attack of the airfoil
%                   velocityInf <double>, the free stream velocity of the flow
%                   pressureInf <double>, the free stream pressure of the flow
%                   N <double>, the number of panels (number of vorticess)
%                   plot <logical>, the option to enable plotting. If true,
%                                   plot all figures. If false, suppress
%                                   plots.
%
%       Outputs:    pOut [100x100] <double>
%                   psiOut [100x100] <double>
%                   phiOut [100x100] <double>
%

function [pOut,psiOut,phiOut]=Plot_Airfoil_Flow(chord,alpha,velocityInf,pressureInf,rhoInf,N,plotEnabled) 
    %% Grid parameters & Constants
    numPoints = 100;
    % Define domain & number of points(x)
    xMin = -2; xMax = 6; xN = numPoints;
    % Define domain & number of points(y)
    yMin = -2; yMax = 2;   yN = numPoints;
    % Create grid over domain
    [x, y] = meshgrid(linspace(xMin, xMax, xN), linspace(yMin, yMax, yN));
    % Get the step size
    stepSize = chord/N;
    % Convert angle of attack to rad
    alpha = (alpha/180)*pi;
    
    %% Initializing Stream function, Potential function and Velocity
    % Initialize a stream function of uniform flow
    airfoilStream = y*velocityInf*cos(alpha) - x*velocityInf*sin(alpha);
    % Initialize a potential function of uniform flow
    airfoilPotential = x*velocityInf*cos(alpha) + y*velocityInf*sin(alpha);
    % Initialize the free stream velocity in x direction
    velU = ones(100, 100)*velocityInf*cos(alpha);
    % Initialize the free stream velocity in y direction
    velV = ones(100, 100)*velocityInf*sin(alpha);
    
    %% Initialize helper functions for calculations
    % Initialize an anonymous function to represent the vortex strength
    getStrength = @(x) 2*alpha*velocityInf*sqrt((1-(x/chord))/(x/chord));
    % Function to represent the radius of the vortex
    radius = @(x, y, x1, y1) sqrt((x-x1).^2+(y-y1).^2);
    % Airfoil chord x and y values (Only used for plotting)
    airfoilChordX = linspace(0, chord, 20); airfoilChordY = zeros(1, 20);

    %% Computation of Stream, Potential, and Pressure
    for index = linspace(stepSize, chord - stepSize, N)
       % Calculate the strength of the vortex
       gamma = getStrength(index) * stepSize;
       % Calculate the stream function 
       psi = (gamma/(2*pi))*log(radius(x, y, index, 0));
       % Calculate the potential function 
       phi = (-gamma*atan2(-y, -(x - index)))/(2*pi);
       % Calculate velocity in the theta direction
       velTheta = gamma./(2*pi*radius(x, y, index, 0));
       
       % Add stream function to total stream function
       airfoilStream = airfoilStream + psi;
       % Add potential function to total potential function
       airfoilPotential = airfoilPotential + phi; 
       % Add velocity in the x direction to total x velocity
       velU = velU+velTheta.*sin(atan2(y,x-index));
       % Add velocity in the y direction to total y velocity
       velV = velV-velTheta.*cos(atan2(y,x-index));
    end
    % Take the magnitude of velocity 
    airfoilVelocity = sqrt(velU.^2 + velV.^2);
    % Calculate total pressure
    airfoilPressure = pressureInf - 0.5*rhoInf*airfoilVelocity.^2;
    % Return pressure, stream, and potential matrices
    pOut = airfoilPressure; psiOut = airfoilStream; phiOut = airfoilPotential;
    
    %% Plotting Pressure Contours, Stream, and Potential 
    if plotEnabled == true
        %-----------------------------------------------------------------------
        % Pressure Contours
        %-----------------------------------------------------------------------
        levelMin = min(airfoilPressure(:));
        levelMax = max(airfoilPressure(:));
        levels = linspace(levelMin, levelMax, 70);
        figure; hold on; box on; axis equal; colormap bone; colorbar
        contourf(x, y, airfoilPressure, levels); 
        plot(airfoilChordX, airfoilChordY, 'r', 'LineWidth', 3)
        xlabel('x [m]');ylabel('y [m]');title("Pressure [Pa] N="+num2str(N))
        legend('Pressure Contours', 'Chord')
        %-----------------------------------------------------------------------
        % Stream Lines
        %-----------------------------------------------------------------------
        figure;
        subplot(2, 1, 1); hold on; box on; axis equal
        contour(x, y, airfoilStream, 35, 'LineWidth', 1.2);
        plot(airfoilChordX, airfoilChordY, 'r', 'LineWidth', 2)
        xlabel('x [m]');ylabel('y [m]');
        title("Stream Function \psi(x,y) N="+num2str(N))
        legend('\psi(x, y)', 'Chord')
        %-----------------------------------------------------------------------
        % Equipotential Lines 
        %-----------------------------------------------------------------------
        subplot(2, 1, 2); hold on; box on; axis equal
        contour(x, y, airfoilPotential, 40, 'LineWIdth', 1.2);
        plot(airfoilChordX, airfoilChordY, 'r', 'LineWidth', 2)
        xlabel('x [m]');ylabel('y [m]');
        title("Potential Function \phi(x,y) N="+num2str(N))
        legend('\phi(x, y)', 'Chord')
    end
end
