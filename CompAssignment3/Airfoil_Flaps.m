%% ASEN 3111 Aerodynamics Computational Assignment 3 - Airfoil_Flaps.m
% This function creates boundary points for a NACA 0012 Airfoil with a
% flap. The geometry of the flap is based on the description in the
% computational assignment 3 document. Both the main airfoil and the flap are
% generated using NACA_Airfoils.m, but the flap is rotated using a rotation
% matrix. The vectors for the main airfoil and the flap are then
% restructured to be compatible with Vortex_Panel.m.
%
%   Author: Andrew Thompson
%   Created: 10/29/2020 Edited: 10/29/2020
%
%   Parameters:     deflection <double> - The deflection of the flap 
%                   chord <double> - The airfoil chord length
%
%   Outputs:        x [N x 1] <double> - X boundary values 
%                   y [N x 1] <double> - Y boundary values

function [x, y] = Airfoil_Flaps(deflection,chord)
    %% Init Geometries
    foilChord = 0.75*chord;
    flapChord = 0.25*chord;
    gap = 0.05*chord;
    
    %% Basic shapes
    [foilX,foilY] = NACA_Airfoils(0,0,0.12,foilChord,100);
    [flapX,flapY] = NACA_Airfoils(0,0,0.12,flapChord,100);

    %% Orienting flap
    % Rotation Matrix
    rotation = @(theta) [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
    % Concatenating flap coordinates
    flapCoords = [flapX flapY];
    % Perform rotation
    flap = flapCoords*rotation(deflection);
    % Shift flap 
    flap(:,1) = flap(:,1) + foilChord + gap;
    
    %% Output
    % Structure arrays in clockwise orientation around airfoil
    x = [flap(1:50,1); foilX; flap(51:end,1)];
    y = [flap(1:50,2); foilY; flap(51:end,2)];
    
end