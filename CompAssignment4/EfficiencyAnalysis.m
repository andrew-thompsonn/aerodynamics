%% ASEN 3111 Aerodynamics Computational Assignment #4 - EfficiencyAnalysis.m
% Analyizes the effect of varying taper ratio on span efficiency factor of
% a finite wing. Performs analysis on varying aspect ratios
%
%   Authors: Andrew Thompson
%   Created: 11/12/2020 Edited: 11/12/2020
%
%   Parameters:     p <struct> - Finite wing geometry and flow parameters
%                   N <int> - Number of odd terms in circulation expansion
%   Returns:        None

function EfficiencyAnalysis(N)
    % Thin Airfoil Theory
    a0R = deg2rad(2*pi);     
    a0T = deg2rad(2*pi);    
    % Taper and AR arrays
    aspect = 4:2:10;
    tapers = linspace(0, 1, 50);
    % Fix span 
    b = 150;
    % Initialize figure
    figure; hold on; grid minor;
    % For all aspect ratios
    for AR = aspect
        % Initialize efficiency array
        e = zeros(1, length(tapers));
        % Initialize index
        index = 1; 
        % For all taper ratios
        for taper = tapers
            % Solve for chord lengths
            aMat = [1, 1; 1, -taper];
            bVec = [2*b/AR; 0];
            chords = aMat\bVec;
            cT = chords(1);
            cR = chords(2);
            % Compute efficiency factor
            e(index) = PLLT(b,a0T,a0R,cT,cR,0,0,5,5,N);
            % Increment index
            index = index + 1;
        end
        % Plot the efficiency factor
        plot(tapers, e, 'LineWidth', 1.5)
    end
    xlabel("C_t/C_r");ylabel("e");
    title("Taper Ratio vs Efficiency Factor")
    legend("AR=4","AR=6","AR=8","AR=10")
end
