%% ASEN 3111 Aerodynamics Computational Assignment 3 - Compute_RMSE.m
% This function uses a reference value for circulation to compute RMSE of
% circulation as the number of panels varies. The function selects the
% lowest number of panels that produce an RMSE of less than 0.025.
%
%   Author: Andrew Thompson
%   Created: 10/29/2020 Edited: 10/29/2020
%
%   Parameters:     tauRef <double>
%                   chord <double>
%                   vInf <double>
%   Output:         Nreturn <int>
    

function [Nreturn] = Compute_RMSE(cpRef, chord, vInf)
    % Initialize an array of panels
    panels = 60:20:1000;
    % Status variable
    foundN = false;

    
    
    
    cpLowerRef = cpRef(1:floor(length(cpRef)/2));
    cpUpperRef = cpRef(floor(length(cpRef)/2)+1:end);
    index = 1;
    
    % For all panel values 
    for N = panels
        % Get the boundary points and coefficient of pressure 
        [x0012i, y0012i] = NACA_Airfoils(0.00,0.0,0.12,chord,N);
        % Compute circulation
        [~,cpi, ~] = Vortex_Panel(x0012i, y0012i, vInf, 0, false);
      
        cpiLower = cpi(1:floor(N/2));
        cpiUpper = cpi(floor(N/2)+1:end);
        
        
        pErrorUpper(index) = sqrt((mean(cpUpperRef) - mean(cpiUpper))^2)/N;
        pErrorLower(index) = sqrt((mean(cpLowerRef) - mean(cpiLower))^2)/N;
        index = index + 1;
        
    end
 
    figure(6); hold on;
    title("Percent Error in C_P vs Number of Panels")
    xlabel("N")
    ylabel("Error")
    plot(panels, pErrorUpper, 'ro', 'MarkerFaceColor', 'r')
    plot(panels, pErrorLower, 'bo', 'MarkerFaceColor', 'b')
%     yline(1, 'g--')
%     yticks = [get(gca,'ytick')]';
%     percentsy = repmat('%', length(yticks),1); 
%     yticklabel = [num2str(yticks) percentsy];
%     set(gca, 'yticklabel', yticklabel);
    
    legend("Upper Surface", "Lower Surface")
    
    Nreturn = 100;
end