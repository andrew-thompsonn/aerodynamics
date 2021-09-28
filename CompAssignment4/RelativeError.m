%% ASEN 3111 Aerodynamics Computational Assignment #4 - RelativeError.m
% Determines the number of odd terms in the circulation series expansion to
% obtain values for lift and drag within 5%, 1%, and 0.1% relative error to
% a given reference value. 
%
%   Author: Andrew Thompson
%   Created: 11/12/2020 Edited: 11/12/2020
%
%   Parameters:     liftRef <double> - Reference value for lift 
%                   dragRef <double> - Reference value for drag
%                   p <struct> - Finite wing geometry and flow parameters
%   Returns:        N5  <int> - Number of odd terms for 5% error
%                   N1  <int> - Number of odd terms for 1% error
%                   N10 <int> - Number of odd terms for 0.1% error

function [N5, N1, N10] = RelativeError(liftRef, dragRef, p)
    % Initialize error for loop
    err1Found = false;
    err2Found = false;
    err3Found = false;
    N = 5;
    % While 1/10th percent relative error has not been found
    while err3Found == false
        % Compute lift and drag with N terms
        [~, cL, cDi] = PLLT(p.b,p.a0T,p.a0R,p.cT,p.cR,p.aeroT,p.aeroR,p.geoT,p.geoR,N);
        lift = cL *0.5*p.rho*p.vInf^2*p.S;
        drag = cDi*0.5*p.rho*p.vInf^2*p.S;
        % Calculate error for lift and drag 
        liftError = 100*abs(lift-liftRef)/liftRef;
        dragError = 100*abs(drag-dragRef)/dragRef;
        % If error less than 5 percent 
        if liftError < 5 && dragError < 5 && err1Found == false
            % Document number of terms
            N5 = N;
            err1Found = true;
        % If error less than 1 percent
        elseif liftError < 1 && dragError < 1 && err2Found == false
            % Document number of terms
            N1 = N;
            err2Found = true;
        % If error less than 1/10 percent
        elseif liftError < 0.1 && dragError < 0.1 && err3Found == false
            % Document number of terms
            N10 = N;
            err3Found = true;
        end
        % Increment number of terms
        N = N + 1;
    end
end