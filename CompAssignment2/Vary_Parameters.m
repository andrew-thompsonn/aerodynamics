%% ASEN 3111 Aerodynamics Computational Assignment 2 - Vary_Parameters.m
% This function takes inputs for varying chord lengths, angles of attack,
% and freestream velocities to calculate and plot the potential and stream
% function for the resulting flows. Creates 3 figures each containing 8
% subplots.
%
%       Author: Andrew Thompson
%       Created: 10/4/2020      Edited: 10/4/2020
%
%       Inputs:     chordArray [1xN] <double>, an array of chord lengths
%                   alphaArray [1xN] <double>, an array of angles of attack
%                   velocityArray [1xN] <double>, an array of velocities
%                   params <struct>, Structure containing fields for parameters
%
%       Ouputs:     None

function Vary_Parameters(chordArray, alphaArray, velocityArray, params)  
    %% Generate a Constant Grid
    numPoints = 100;
    xMin = -2; xMax = 6; xN = numPoints;
    yMin = -2; yMax = 2;   yN = numPoints;
    [x, y] = meshgrid(linspace(xMin, xMax, xN), linspace(yMin, yMax, yN));
    
    % Baseline 
    [~,psi,phi]=Plot_Airfoil_Flow(params.chord,params.alpha,params.velocityInf, ...
        params.pressureInf, params.rhoInf,1000,false);   
    
    % Create consistent levels for all plots
    psiLevMin = min(psi(:))*2;
    psiLevMax = max(psi(:))*2;
    phiLevMin = min(phi(:))*2;
    phiLevMax = max(psi(:))*4;
    psiLevels = linspace(psiLevMin, psiLevMax, 40);
    phiLevels = linspace(phiLevMin, phiLevMax, 40);
    
    %% Vary Chord Length
    index = 1;
    figure;
    % For all chord values
    for chord = chordArray
        % Airfoil chord X and Y values (For plotting only)
        airfoilChordX = linspace(0, chord, 20); airfoilChordY = zeros(1, 20);
        % Calculate the velocity stream and potential
        [~,psi,phi]=Plot_Airfoil_Flow(chord,params.alpha,params.velocityInf, ...
            params.pressureInf, params.rhoInf,1000,false); 
        % Plot the Stream function
        subplot(length(chordArray), 2, index)
        contour(x, y, psi, psiLevels); hold on;
        plot(airfoilChordX, airfoilChordY, '-r', 'LineWidth', 2)
        xlabel("     x [m]"); ylabel("y [m]"); title("\psi c = "+num2str(chord))
        xlim([-2 6]); ylim([-2 2])
        index = index + 1;
        % Plot the potential function
        subplot(length(chordArray), 2, index)
        contour(x, y, phi, phiLevels); hold on
        plot(airfoilChordX, airfoilChordY, '-r', 'LineWidth', 2)
        xlabel("     x [m]"); ylabel("y [m]"); title("\phi c = "+num2str(chord))
        xlim([-2 6]); ylim([-2 2])
        index = index + 1;
    end
    
    %% Vary Angle of Attack
    index = 1;
    figure;
    % Airfoil chord X and Y values (For plotting only)
    airfoilChordX = linspace(0, params.chord, 20); airfoilChordY = zeros(1, 20);
    % For all angles of attack
    for alpha = alphaArray
        % Calculate the velocity stream and potential
        [~,psi,phi]=Plot_Airfoil_Flow(params.chord,alpha,params.velocityInf, ...
            params.pressureInf, params.rhoInf,1000,false); 
        % Plot the Stream function
        subplot(length(alphaArray), 2, index)
        contour(x, y, psi, psiLevels); hold on;
        plot(airfoilChordX, airfoilChordY, '-r', 'LineWidth', 2)
        xlabel("     x [m]"); ylabel("y [m]"); title("\psi \alpha = "+num2str(alpha)+char(176))
        xlim([-2 6]); ylim([-2 2])
        index = index + 1;
        % Plot the potential function
        subplot(length(alphaArray), 2, index)
        contour(x, y, phi, phiLevels); hold on
        plot(airfoilChordX, airfoilChordY, '-r', 'LineWidth', 2)
        xlabel("     x [m]"); ylabel("y [m]"); title("\phi \alpha = "+num2str(alpha)+char(176))
        xlim([-2 6]); ylim([-2 2])
        index = index + 1;
    end
    
    %% Vary Velocity
    index = 1;
    figure;
    % For all free stream velocities
    for velocityInf = velocityArray
        % Calculate the velocity stream and potential
        [~,psi,phi]=Plot_Airfoil_Flow(params.chord,params.alpha,velocityInf, ...
            params.pressureInf, params.rhoInf,1000,false); 
        % Plot the Stream function
        subplot(length(velocityArray), 2, index)
        contour(x, y, psi, psiLevels); hold on;
        plot(airfoilChordX, airfoilChordY, '-r', 'LineWidth', 2)
        xlabel("     x [m]"); ylabel("y [m]"); title("\psi V = "+num2str(velocityInf)+" m/s")
        xlim([-2 6]); ylim([-2 2])
        index = index + 1;
        % Plot the potential function
        subplot(length(velocityArray), 2, index)
        contour(x, y, phi, phiLevels); hold on
        plot(airfoilChordX, airfoilChordY, '-r', 'LineWidth', 2)
        xlabel("     x [m]"); ylabel("y [m]"); title("\phi V = "+num2str(velocityInf)+" m/s")
        xlim([-2 6]); ylim([-2 2])
        index = index + 1;
    end
end