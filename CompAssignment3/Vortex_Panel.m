%% ASEN 3111 Aerodynamics Computational Assignment 3 - Vortex_Panel.m
% The purpose of this function is to use the vortex panel method to
% determine coefficient of lift, and plot coefficient of pressure around an
% airfoil at a specific angle of attack and velocity. 
%
% This function was developed by following the process in Keuthe & Chow's 
% Fortran Vortex Panel implementation.
%
%   Author: Andrew Thompson
%   Created: 10/9/2020  Edited: 10/26/2020
%   Note: With plotting enabled, a figure must be created prior to the
%         function call
%
%   Parameters:     xb [N x 1] <double> - x locations of boundary points
%                   yb [N x 1] <double> - y locations of boundary points
%                   vInf <double> - Free stream velocity
%                   alpha <double> - Angle of attack
%                   plotStatus <bool> - Option to supress output plots
%
%   Outputs:        cl <double> - Coefficient of lift
%                   cp [N x 1] <double> - Array of coefficients of pressure


function [cl,cp,tau] = Vortex_Panel(xb, yb, vInf, alpha, plotStatus)
    %% Initializing variables and constants
    m = length(xb)-1;
    mp1 = m+1;
    alpha = deg2rad(alpha);
    
    x =     zeros(m,1);
    y =     zeros(m,1);
    s =     zeros(m,1);
    sine =  zeros(m,1);
    cosine= zeros(m,1);
    theta = zeros(m,1);
    RHS =   zeros(m,1);
    cn1 =   zeros(m,m);
    cn2 =   zeros(m,m);
    ct1 =   zeros(m,m);
    ct2 =   zeros(m,m);
    an =    zeros(mp1,mp1);
    at =    zeros(m,mp1);

    %% Compute Panel Information
    for index = 1:m
        % Get 2cnd point 
        indexP1 = index + 1;
        % Get x reference value 
        x(index) = 0.5*(xb(index)+xb(indexP1));
        % Get y reference value 
        y(index) = 0.5*(yb(index)+yb(indexP1));
        % Compute panel length
        s(index) = sqrt((xb(indexP1)-xb(index))^2+(yb(indexP1)-yb(index))^2);
        % Compute angle theta at boundaries
        theta(index) = atan2(yb(indexP1)-yb(index),xb(indexP1)-xb(index));
        % Create array of sin values 
        sine(index) = sin(theta(index));
        % Create array of cosine values
        cosine(index) = cos(theta(index));
        % Create RHS of equation
        RHS(index) = sin(theta(index)-alpha);
    end
    
    %% Constructing 'cn2', 'cn1', 'ct2', 'ct1' Matricies
    for i = 1:m
        for j = 1:m
            % If on the diagonal
            if i == j
                cn1(i,j) = -1;
                cn2(i,j) = 1;
                ct1(i,j) = 0.5*pi;
                ct2(i,j) = 0.5*pi;
            % If not on diagonal
            else
                % Compute constants
                A = -(x(i)-xb(j))*cosine(j)-(y(i)-yb(j))*sine(j);
                B = (x(i)-xb(j))^2 + (y(i)-yb(j))^2;
                C = sin(theta(i)-theta(j));
                D = cos(theta(i)-theta(j));
                E = (x(i)-xb(j))*sine(j)-(y(i)-yb(j))*cosine(j);
                F = log(1+s(j)*(s(j)+2.*A)/B);
                G = atan2(E*s(j), B+A*s(j));
                P = (x(i)-xb(j))*sin(theta(i)-2.*theta(j))+(y(i)-yb(j))* ...
                    cos(theta(i)-2.*theta(j));
                Q = (x(i)-xb(j))*cos(theta(i)-2.*theta(j))-(y(i)-yb(j))* ...
                    sin(theta(i)-2.*theta(j));
                % Use constants to compute values for cn1,cn2,ct1,ct2
                cn2(i,j) = D+0.5*Q*F/s(j)-(A*C+D*E)*G/s(j);
                cn1(i,j) = 0.5*D*F+C*G -cn2(i,j);
                ct2(i,j) = C+0.5*P*F/s(j)+(A*D-C*E)*G/s(j);
                ct1(i,j) = 0.5*C*F-D*G-ct2(i,j);
            end  
        end
    end
    %% Constructing 'An' Matrix in Equation (5.47)
    for i = 1:m
        % Build matrix An & At with cn1,cn2,ct1,ct2
        an(i,1)   = cn1(i,1);
        an(i,mp1) = cn2(i,m);
        at(i,1)   = ct1(i,1);
        at(i,mp1) = ct2(i,m);
        for j = 2:m
            an(i,j) = cn1(i,j)+cn2(i,j-1);
            at(i,j) = ct1(i,j)+ct2(i,j-1);
        end
    end
    
    an(mp1,1) = 1;
    an(mp1,mp1) = 1;

    for j = 2:m
        an(mp1,j) = 0;
    end

    RHS(mp1) = 0;
    %% Solving Equation (5.47) & Computing Pressure/Lift
    % Compute non dimensional sheet strength by solving the matrix equation
    gamma = an\RHS;
    % Compute non dimensional tangential velocity
    v = cos(theta - alpha) + at*gamma;
    % Compute coefficient of pressure 
    cp = 1-v.^2;
    % Get the chord length
    chord = max(xb)-min(xb);
    % Compute the total, dimensionalized circulation
    circulation = v.*s*vInf;
    sum(circulation)
    % Use circulation to compute coefficient of lift 
    cl = sum(2*circulation)/(vInf*chord);
    tau = sum(circulation);
    theta

    %% Plotting Coefficient of Pressure
    N = length(cp);
    if plotStatus == true
        hold on;
        grid minor;
        % Reverse y axis
        set(gca, 'YDir','reverse')
        % Plot x normalized by chord length
        p1 = plot(x(1:ceil(N/2))/chord, cp(1:ceil(N/2)), 'b', 'LineWidth', 1.5);
        p2 = plot(x(ceil(N/2):end)/chord, cp(ceil(N/2):end), 'r', 'LineWidth', 1.5);
        yline(0, '--');
        legend([p1,p2],["Lower Surface", "Upper Surface"])
        ylim([-4.5 1])
        xlabel("x/c")
        ylabel("C_p")
        title("Coefficient of Pressure \alpha = "+num2str(rad2deg(alpha))+char(176))
    end

end

