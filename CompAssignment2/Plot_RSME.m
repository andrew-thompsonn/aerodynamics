%% ASEN 3111 Aerodynamics Computational Assignment 2 - Plot_RSME.m
% This function calculates and plots the root-mean-square error of the pressure
% field, stream function, and potential function based on reference values
% for the three fields. 
%
%       Author: Andrew Thompson
%       Created: 10/2/2020      Edited: 10/2/2020
%
%       Inputs:     pRef [NxN] <double>
%                   psiRef [NxN] <double>
%                   phiRef [NxN] <double>
%
%       Outputs:    None

function Plot_RSME(pRef, psiRef, phiRef, vRef, params)
    % Create a figure
    figure; hold on; box on; grid minor
    xlabel('N'); ylabel('RSME');title('Error Analysis of Velocity and Pressure')
    panels = linspace(10, 300, 30);

    % For number of panels between 10 and 300
    for N = panels
        % Calculate p, psi, phi with plotting disabled for current N value
        [p,psi,phi]=Plot_Airfoil_Flow(params.chord,params.alpha, ...
            params.velocityInf,params.pressureInf,params.rhoInf,N,false);
        
        % Calculate magnitude of velocity using bernoulli's
        v = sqrt((2/params.rhoInf)*(params.pressureInf-p));
        
        % Calculate mean RMSE
        rmsP = mean2(sqrt((pRef - p).^2)/N);
        rmsV = mean2(sqrt((vRef - v).^2)/N);
        rmsPsi = mean2(sqrt((psiRef - psi).^2)/N);
        rmsPhi = mean2(sqrt((phiRef - phi).^2)/N);
        
        % Plot Points
        plot(N, rmsP, 'ks', 'MarkerFaceColor', 'r')
        plot(N, rmsV, 'kv', 'MarkerFaceColor', 'b')
        
    end
    legend('Pressure', 'Velocity')
    
end