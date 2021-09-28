%% ASEN 3111 Computational Assignment 3 - main.m
% The purpose of this script is to perform a vortex panel method on several
% different airfoils to calculate the coefficient of lift. The results for
% airfoils of varying camber and thickness are then compared. To select an
% appropriate number of panels, the RMSE of total circulation
% was compared to a reference value computed  using 1000 panels.

%% Clean
clc; close all; clear all;

%% Constants
alpha = linspace(-10, 10, 21);
vInf = 30;
chord = 1;
N = 1000;
plotDisabled = false;
plotEnabled = true;

%% Select N Using RMSE Analysis, Reference Value = 1000 panels
[x0012, y0012] = NACA_Airfoils(0.00, 0.0, 0.12, chord, N);
[~,cpRef, ~] = Vortex_Panel(x0012,y0012,vInf,0,plotDisabled);
N = Compute_RMSE(cpRef, chord, vInf);
fprintf("Selected N for circulation RMSE < 0.025: %i\n",N);

%% Airfoil Boundary Points 
[x0012, y0012] = NACA_Airfoils(0.00, 0.0, 0.12, chord, N);
[x2412, y2412] = NACA_Airfoils(0.02, 0.4, 0.12, chord, N);
[x4412, y4412] = NACA_Airfoils(0.04, 0.4, 0.12, chord, N);
[x2424, y2424] = NACA_Airfoils(0.02, 0.4, 0.24, chord, N);

%% Plot NACA 0012 Coefficient of Pressure at -5, 0, 5, 10 AoA
figure(1);
subplot(2,2,1); Vortex_Panel(x0012, y0012, vInf,-5, plotEnabled);
subplot(2,2,2); Vortex_Panel(x0012, y0012, vInf, 0, plotEnabled);
subplot(2,2,3); Vortex_Panel(x0012, y0012, vInf, 5, plotEnabled);
subplot(2,2,4); Vortex_Panel(x0012, y0012, vInf,10, plotEnabled);

%% Initialize Arrays
cL0012 = zeros(1, length(alpha));
cL2412 = zeros(1, length(alpha));
cL4412 = zeros(1, length(alpha));
cL2424 = zeros(1, length(alpha));

%% Calculate Lift at Varying AoA
for index = 1:length(alpha)
    cL0012(index) = Vortex_Panel(x0012,y0012,vInf,alpha(index),plotDisabled);
    cL2412(index) = Vortex_Panel(x2412,y2412,vInf,alpha(index),plotDisabled);
    cL4412(index) = Vortex_Panel(x4412,y4412,vInf,alpha(index),plotDisabled);
    cL2424(index) = Vortex_Panel(x2424,y2424,vInf,alpha(index),plotDisabled);
end
%% Calculate Zero Lift AoA
alpha0012 = Zero_Lift_Alpha(alpha, cL0012);
alpha2412 = Zero_Lift_Alpha(alpha, cL2412);
alpha4412 = Zero_Lift_Alpha(alpha, cL4412);
alpha2424 = Zero_Lift_Alpha(alpha, cL2424);

%% Plot CL vs AoA & Zero lift AoA for 0012, 2412, 4412, 2424
figure(2); hold on; grid minor; 
% Coefficient of lift vs alpha 
p1 = plot(alpha, cL0012, 'b', 'LineWidth', 1.2);
p2 = plot(alpha, cL2412, 'r', 'LineWidth', 1.2);
p3 = plot(alpha, cL4412, 'k', 'LineWidth', 1.2);
p4 = plot(alpha, cL2424, 'm', 'LineWidth', 1.2);
% Zero lift angles of attack
a1 = plot(alpha0012, 0, 'bs', 'MarkerFaceColor', 'b');
a2 = plot(alpha2412, 0, 'rs', 'MarkerFaceColor', 'r');
a3 = plot(alpha4412, 0, 'ks', 'MarkerFaceColor', 'k');
a4 = plot(alpha2424, 0, 'ms', 'MarkerFaceColor', 'm');
% Title, labels, legend, axes
l1 = xline(0, '-.'); 
l2 = yline(0, '-.');
xlabel("\alpha"+char(176)); ylabel("C_L");
title("Coefficient of Lift vs Angle of Attack")
legend([p1,p2,p3,p4],["NACA 0012","NACA 2412","NACA 4412","NACA 2424"])

%% BONUS
% Get boundary points for 0012 with 30 deg deflection flap
[xFlap,yFlap] = Airfoil_Flaps(30,1);
% Initialize array of cl
clFlap = zeros(1, length(alpha));
% For all alpha 
for index = 1:length(alpha)
    % Get cl with flap 
    clFlap(index) = Vortex_Panel(xFlap,yFlap,30,alpha(index), false);
end
figure(4); hold on; grid minor;
% Plotting cls
p5 = plot(alpha, clFlap, 'r', 'LineWidth', 1.2);
p6 = plot(alpha, cL0012, 'b', 'LineWidth', 1.2);
% Title, labels, legend, axes
l3 = xline(0, '-.'); 
l4 = yline(0, '-.');
xlabel("\alpha"+char(176)); ylabel("C_L")
title("Coefficient of Lift vs Angle of Attack")
legend([p5, p6],["NACA 0012 With Flap", "NACA 0012 Without Flap"])
% Plotting airfoil with flap (Just for fun)
figure(5); hold on; grid minor; axis equal;
plot(xFlap, yFlap, 'b', 'LineWidth', 1.2)
xlabel("x [m]"); ylabel("y [m]")
title("NACA 0012 With Flap")

%%



