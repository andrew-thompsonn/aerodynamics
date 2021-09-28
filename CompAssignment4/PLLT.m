%% ASEN 3111 Aerodynamics CA4 - PLLT.m
% This function solves the fundamental equation of Prandtl Lifting Line
% Theory for finite wings with thick airfoils. Linearly varying arrays of
% lift curve slope, chord lenght, zero lift angle of attack, and
% geometric angle of attack are generated and used to develop a system of
% equations. This system is solved for A_0 - A_n to calculate lift, induced
% drag, and wing efficiency. balls
%
%   Author: Andrew Thompson
%   Collaborator: Jacob Pendergast
%   Created: 11/06/20 Edited: 11/06/20
%
%   Parameters:     b     <double> - wing span
%                   a0T   <double> - lift curve slope (tip) [/deg]
%                   a0R   <double> - lift curve slope (root)[/deg]
%                   cT    <double> - chord length (tip)
%                   cR    <double> - chord length (root) [deg]
%                   aeroT <double> - zero lift AoA (tip)
%                   aeroR <double> - zero lift AoA (root)
%                   geoT  <double> - geometric AoA (tip)
%                   geoR  <double> - geometric AoA (root)
%                   N     <int>    - number of odd terms 
%   Returns:        e     <double> - span efficiency
%                   cL    <double> - coefficient of lift 
%                   cDi   <double> - induced coefficient of drag
%
function [e,cL,cDi] = PLLT(b,a0T,a0R,cT,cR,aeroT,aeroR,geoT,geoR,N)
    %% Unit conversions
    a0R   = rad2deg(a0R);       a0T   = rad2deg(a0T);
    aeroR = deg2rad(aeroR);     aeroT = deg2rad(aeroT);
    geoR  = deg2rad(geoR);      geoT  = deg2rad(geoT);
    
    %% Spanwise functions of theta
    a0    = @(t) a0R   + (a0T   - a0R)  * cos(t);
    chord = @(t) cR    + (cT    -  cR)  * cos(t);
    aero  = @(t) aeroR + (aeroT - aeroR)* cos(t);
    geo   = @(t) geoR  + (geoT  - geoR) * cos(t);
    
    %% Constants
    S = b*(cR + cT)/2;
    AR = b^2/S;
    
    %% Build matrix of coefficients
    % Initialize matrix
    aMatrix = zeros(N);
    % Initialize RHS 
    RHS = zeros(N, 1);
    % For every row 
    for i = 1:N
         % Define current theta
         theta = i*pi/(2*N);
         % Fill RHS of the equation
         RHS(i) = geo(theta) - aero(theta);
         % Get the chord length from current theta
         c = chord(theta);
         % Get the lift slope from current theta
         a = a0(theta);
         % For every odd term
         for j = 1:N
             % Define current step
             n = 2*j-1;
             % Fill coefficients 
             aMatrix(i,j) =  4*b*sin(n*theta)/(a*c)+n*(sin(n*theta))/sin(theta);
         end
    end
    
    %% Solving for coefficients & Computing lift, drag, efficiency
    % Solve the matrix equation
    An = aMatrix\RHS; 
    % Compute Coefficient of Lift
    cL = pi*AR*An(1);
    % Initialize delta
    delta = 0;
    % For all odd terms 
    for i = 2:N
        n = 2*i-1;
        % Increment delta
        delta = delta + n*(An(i)/An(1))^2;
    end
    % Compute Span efficiency 
    e = 1/(1+delta);
    % Compute Coefficient of Induced Drag
    cDi = cL^2/(pi*AR*e);
end
