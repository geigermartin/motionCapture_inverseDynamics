%% INFORMATION

%   This script can be used to calaculate the ankle, knee, and hip joint forces and torques
%   with measured:
%   - body weight [kg]
%   - joint (toe (D), ankle (A), knee (B), hip (C)) coordinates [m]
%   - COP coordinates [m]
%   - ground reaction force [N]

% ******************************************************************************
%%
clear; close all; clc;

%% Fix parameters 

% gravitational acceleration [m/s^2]
g = 10;

% attachment point of spring mass damper 1 (SMD1) proximal from segment COM
att_prox_COM = 0.3; 

%% Parameters changing with subject

% body mass [kg] 
m = 60; 

% relative COM-position from proximal joint (based on de Leva 1996)
% female
COM_rel_f = 0.4014 ; 
COM_rel_s = 0.4416;
COM_rel_t = 0.3612;
% male
% COM_rel_f = 0.4415;
% COM_rel_s = 0.4459;
% COM_rel_t = 0.4095;

% relative segment mass of one foot/shank/thigh (based on Zatsiorsky 1990a, 1990b and 1993)
% female
m_rel_f = 0.0129;
m_rel_s = 0.0481;
m_rel_t = 0.1478;
% male
% m_rel_f = 0.0137;
% m_rel_s = 0.0433;
% m_rel_t = 0.1416;

%% Parameter changing with measurement settings

% sampling rates [Hz]
sr_Kinematics = 250;
sr_ForcePlate = 1000;

%% Parameters depending on time

test = 0; % 1... fixed values of time-dependent parameters for testing

if test == 1
    
    % joint and COP coordinates [m]
    x_A = 0.15;
    y_A = 0.15;
    x_B = 0.21;
    y_B = 0.55;
    x_C = 0.15;
    y_C = 0.93;
    x_D = 0.35;
    y_D = 0;
    x_COP_long = 0.33;
    y_COP_long = 0;
    
    % COM acceleration [m/s^2]
    a_x_COM_f = 3;
    a_y_COM_f = 2;
    a_x_COM_s = -1;
    a_y_COM_s = 4;
    a_x_COM_t = 0;
    a_y_COM_t = 4.5;
    
    % angular acceleration [rad/r^2]
    alpha_f = -50;
    alpha_s = 25;
    alpha_t = -15;
    
    % ground reaction force [N]
    F_x_long = 200;
    F_y_long = 800;   
    
else
    
    % import model joint data and divide by 1000 -> [m]
    % left leg (right leg --> column K:R)
    x_C = xlsread('Gelenksdaten','C:C')/1000;
    y_C = xlsread('Gelenksdaten','D:D')/1000;
    x_B = xlsread('Gelenksdaten','E:E')/1000;
    y_B = xlsread('Gelenksdaten','F:F')/1000;
    x_A = xlsread('Gelenksdaten','G:G')/1000;
    y_A = xlsread('Gelenksdaten','H:H')/1000;
    x_D = xlsread('Gelenksdaten','I:I')/1000;
    y_D = xlsread('Gelenksdaten','J:J')/1000;
    
    % import GRF [N] and divide by 2 (only one leg), import COP data and divide by 1000 -> [m]
    F_x_long = xlsread('GRFundCoP','C:C')/2;
    F_y_long = xlsread('GRFundCoP','D:D')/2;
    x_COP_long = xlsread('GRFundCoP','E:E')/1000;
    y_COP_long = xlsread('GRFundCoP','F:F')/1000;
    
    % filter force and COP data 
    F_x_long = movmean(F_x_long,11);
    F_y_long = movmean(F_y_long, 11);
    x_COP_long = movmean(x_COP_long, 11);
    y_COP_long = movmean(y_COP_long, 11);
    
    % preallocate empty variables
    F_x = zeros(length(x_A), 1);
    F_y = zeros(length(x_A), 1);
    x_COP = zeros(length(x_A), 1);
    y_COP = zeros(length(x_A), 1);
    
    % reduce number of force and COP data points 
    sr_ratio = sr_ForcePlate/sr_Kinematics;
    for ii = 1 : length(x_A) 
        F_x(ii, 1) = F_x_long(sr_ratio * ii - (sr_ratio-1));
        F_y(ii, 1) = F_y_long(sr_ratio * ii - (sr_ratio-1));
        x_COP(ii, 1) = x_COP_long(sr_ratio * ii - (sr_ratio-1));
        y_COP(ii, 1) = y_COP_long(sr_ratio * ii - (sr_ratio-1));
    end    
    
    % pick first ground contact phase (threshold: nearly 5 % of the body weight)
    e_cp = find(diff(F_y > 0.5*m));
    e_touchdown = e_cp(1); 
    e_takeoff = e_cp(2); 
    
    % shorten model joint data (first touchdown to first takeoff)
    x_C = x_C(e_touchdown:e_takeoff);
    y_C = y_C(e_touchdown:e_takeoff);
    x_B = x_B(e_touchdown:e_takeoff);
    y_B = y_B(e_touchdown:e_takeoff);
    x_A = x_A(e_touchdown:e_takeoff);
    y_A = y_A(e_touchdown:e_takeoff);
    x_D = x_D(e_touchdown:e_takeoff);
    y_D = y_D(e_touchdown:e_takeoff);
    
    % shorten GRF and COP data (first touchdown to first takeoff)
    F_x = F_x(e_touchdown:e_takeoff);
    F_y = F_y(e_touchdown:e_takeoff);
    x_COP = x_COP(e_touchdown:e_takeoff);
    y_COP = y_COP(e_touchdown:e_takeoff);
    
end

% ******************************************************************************
%%      CALCULATIONS

%% Calculate segment mass, COM coordinates and moments of inertia
    
% calculate segment mass of one foot (m_f), shank (m_s), and thigh (m_t) [kg]
m_f = m_rel_f * m;
m_s = m_rel_s * m;
m_t = m_rel_t * m;

% calculate COM coordinates from measured joint coordinates [m]
x_COM_f = x_A + COM_rel_f * (x_D - x_A);
y_COM_f = y_A + COM_rel_f * (y_D - y_A);
x_COM_s = x_B + COM_rel_s * (x_A - x_B);
y_COM_s = y_B + COM_rel_s * (y_A - y_B);
x_COM_t = x_C + COM_rel_t * (x_B - x_C);
y_COM_t = y_C + COM_rel_t * (y_B - y_C);

% calculate segment lengths from measured joint coordinates [m]
l_f = sqrt( (x_A - x_D).^2 + (y_A - y_D).^2 );
l_s = sqrt( (x_B - x_A).^2 + (y_B - y_A).^2 );
l_t = sqrt( (x_C - x_B).^2 + (y_C - y_B).^2 );

% calculate moments of inertia [kg * m^2]
J_f = (l_f * (COM_rel_f/2)).^2 * 0.5 * m_f + (l_f * ( (1- COM_rel_f)/2)).^2 * 0.5 * m_f;
J_s = (l_s * (COM_rel_s/2)).^2 * 0.5 * m_s + (l_s * ( (1- COM_rel_s)/2)).^2 * 0.5 * m_s;
J_t= (l_t * (COM_rel_t/2)).^2 * 0.5 * m_t + (l_t * ( (1- COM_rel_t)/2)).^2 * 0.5 * m_t;

%% Calculate COM and angular acceleration  

if test  ~= 1

    % calculate vertical and horizontal COM acceleration (numerical approximation)
    a_x_COM_f = diff( diff(x_COM_f)/(1/sr_Kinematics) )/(1/sr_Kinematics);
    a_y_COM_f = diff( diff(y_COM_f)/(1/sr_Kinematics) )/(1/sr_Kinematics);
    a_x_COM_s = diff( diff(x_COM_s)/(1/sr_Kinematics) )/(1/sr_Kinematics);
    a_y_COM_s = diff( diff(y_COM_s)/(1/sr_Kinematics) )/(1/sr_Kinematics);
    a_x_COM_t = diff( diff(x_COM_t)/(1/sr_Kinematics) )/(1/sr_Kinematics);
    a_y_COM_t = diff( diff(y_COM_t)/(1/sr_Kinematics) )/(1/sr_Kinematics);

    % determine angles between segments
    phi_f = atan2( (y_D - y_A),  (x_D - x_A) );
    phi_s = atan2( (x_A - x_B),  (y_B - y_A) );
    phi_t = atan2( (x_B - x_C),  (y_C - y_B) );

    % calculate angular acceleration
    alpha_f = diff( diff(phi_f)/(1/sr_Kinematics) )/(1/sr_Kinematics);
    alpha_s = diff( diff(phi_s)/(1/sr_Kinematics) )/(1/sr_Kinematics);
    alpha_t = diff( diff(phi_t)/(1/sr_Kinematics) )/(1/sr_Kinematics);

end 

% --------------------------------------------------------------------
%%       Wobbling mass model

%% Calculate attachment points 

% calculate attachment points of SMD1 for shank and thigh
x_att_s = x_COM_s + att_prox_COM * (x_B - x_A);
y_att_s = x_COM_s + att_prox_COM * (y_B - y_A);
x_att_t = x_COM_t + att_prox_COM * (x_C - x_B);
y_att_t = x_COM_t + att_prox_COM * (y_C - y_B);

%% Calculate spring and damper forces of SMD1 and SMD2 for shank and thigh

% calculate simulation time and define time vector for from workspace block
sim_time = length(x_A)/sr_Kinematics; 
time = linspace(0, sim_time, length(x_A))';

% call WM model and calculate spring and damper forces for shank
p.xy_mass_0 = [x_COM_s(1) - 0.01, y_COM_s(1) - 0.01];
p.xy_mass_dot_0 = [0 0]; 
p.k = 2500;
p.c = 200;
p.m_wob = 1/2 * m_s; 
p.xy_att1 = [time, x_att_s, y_att_s];
p.xy_COM = [time, x_COM_s, y_COM_s];
out = sim('WobblingMassModel', [0 sim_time]);
time_s= out.tout;
FK1_s = interp1(time_s, out.FK1, time);
FD1_s = interp1(time_s, out.FD1, time);
FK2_s = interp1(time_s, out.FK2, time);
FD2_s = interp1(time_s, out.FD2, time);
xy_m_wob_s = interp1(time_s, out.xy_m_wob, time);

% call WM model and calculate spring and damper forces for thigh
p.xy_mass_0 = [x_COM_t(1) - 0.01, y_COM_t(1) - 0.01];
p.xy_mass_dot_0 = [0 0]; 
p.m_wob = 2/3 * m_t;
p.xy_att1 = [time, x_att_t, y_att_t];
p.xy_COM = [time, x_COM_t, y_COM_t];
out = sim ('WobblingMassModel', [0 sim_time]);
time_t= out.tout;
FK1_t = interp1(time_t, out.FK1, time);
FD1_t = interp1(time_t, out.FD1, time);
FK2_t = interp1(time_t, out.FK2, time);
FD2_t = interp1(time_t, out.FD2, time);
xy_m_wob_t = interp1(time_t, out.xy_m_wob, time);
    
%% Calculate ankle, knee and hip joint forces [N] and torques [Nm]

% separate x- and y-component of spring, damper force & wobbling-mass-postion vectors
% shank
FK1_s_x = FK1_s(:,1);
FK1_s_y = FK1_s(:,2);
FD1_s_x = FD1_s(:,1);
FD1_s_y = FD1_s(:,2);
FK2_s_x = FK2_s(:,1);
FK2_s_y = FK2_s(:,2);
FD2_s_x = FD2_s(:,1);
FD2_s_y = FD2_s(:,2);
x_m_wob_s = xy_m_wob_s(:,1);
y_m_wob_s = xy_m_wob_s(:,2);
% thigh
FK1_t_x = FK1_t(:,1);
FK1_t_y = FK1_t(:,2);
FD1_t_x = FD1_t(:,1);
FD1_t_y = FD1_t(:,2);
FK2_t_x = FK2_t(:,1);
FK2_t_y = FK2_t(:,2);
FD2_t_x = FD2_t(:,1);
FD2_t_y = FD2_t(:,2);
x_m_wob_t = xy_m_wob_t(:,1);
y_m_wob_t = xy_m_wob_t(:,2);

% preallocate empty variables
FA_x = zeros(length(a_x_COM_f),1);
FA_y = zeros(length(a_x_COM_f), 1);
FB_x = zeros(length(a_x_COM_f), 1);
FB_y = zeros(length(a_x_COM_f), 1);
FC_x = zeros(length(a_x_COM_f), 1);
FC_y = zeros(length(a_x_COM_f), 1);
MA = zeros(length(a_x_COM_f), 1);
MB = zeros(length(a_x_COM_f), 1);
MC = zeros(length(a_x_COM_f), 1);

% calculate joint forces and torques
for t = 1 : length(a_x_COM_f) 
    
    % ankle joint
    FA_x(t) = m_f * a_x_COM_f(t) - F_x(t);
    FA_y(t) = m_f * a_y_COM_f(t) - F_y(t) + m_f * g;
    MA(t) = J_f(t) * alpha_f(t) - ( (x_COP(t) - x_COM_f(t)) * F_y(t) - (y_COP(t) - y_COM_f(t)) * F_x(t) ) ...
        - ( (x_A(t) - x_COM_f(t)) * FA_y(t) - (y_A(t) - y_COM_f(t)) * FA_x(t) );

    % knee joint
    FB_x(t) = m_s * a_x_COM_s(t) + FA_x(t) + FK1_s_x(t) + FD1_s_x(t) + FK2_s_x(t) + FD2_s_x(t);
    FB_y(t) = m_s * a_y_COM_s(t) + FA_y(t) + m_s * g + FK1_s_y(t) + FD1_s_y(t) + FK2_s_y(t) + FD2_s_y(t);
    MB(t) = J_s(t) * alpha_s(t) + MA(t) - ( (x_A(t) - x_COM_s(t)) * (-FA_y(t)) - (y_A(t) - y_COM_s(t)) * (-FA_x(t)) ) ...
        - ( (x_B(t) - x_COM_s(t)) * FB_y(t) - (y_B(t) - y_COM_s(t)) * FB_x(t) ) ...
        - ( (x_att_s(t) - x_COM_s(t)) * -FK1_s_y(t) - (y_att_s(t) - y_COM_s(t)) * -FK1_s_x(t) ) ... % spring
        - ( (x_att_s(t) - x_COM_s(t)) * -FD1_s_y(t) - (y_att_s(t) - y_COM_s(t)) * -FD1_s_x(t) );    % damper

    % hip joint
    FC_x(t) = m_t * a_x_COM_t(t) + FB_x(t) + FK1_t_x(t) + FD1_t_x(t) + FK2_t_x(t) + FD2_t_x(t);
    FC_y(t) = m_t * a_y_COM_t(t) + FB_y(t) + m_t * g + FK1_t_y(t) + FD1_t_y(t) + FK2_t_y(t) + FD2_t_y(t);
    MC(t) = J_t(t) * alpha_t(t) + MB(t) - ( (x_B(t) - x_COM_t(t)) * (-FB_y(t)) - (y_B(t) - y_COM_t(t)) * (-FB_x(t)) ) ...
        - ( (x_C(t) - x_COM_t(t)) * FC_y(t) - (y_C(t) - y_COM_t(t)) * FC_x(t) ) ...
        - ( (x_att_t(t) - x_COM_t(t)) * -FK1_t_y(t) - (y_att_t(t) - y_COM_t(t)) * -FK1_t_x(t) ) ... % spring
        - ( (x_att_t(t) - x_COM_t(t)) * -FD1_t_y(t) - (y_att_t(t) - y_COM_t(t)) * -FD1_t_x(t) );    % damper

end

% ******************************************************************************
%%      OUTPUT

%% Save joint forces and torques

% create matrices with forces and torques for each segment within structure indyn 
invdyn.ankle_Fx_Fy_M = [FA_x, FA_y, MA];
invdyn.knee_Fx_Fy_M = [FB_x, FB_y, MB];
invdyn.hip_Fx_Fy_M = [FC_x, FC_y, MC];

% save data in invdyn.mat 
save invdyn.mat invdyn

%% Plot forces and torques of wobbling masses

% close open figures
close

% use latex-Interpreter
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% use desired colors 
b1 = [0.301 0.745 0.933]; % light blue
b2 = [0.08 0.25 0.652]; % dark blue

% --------------------------------------------------------------------
% Figure 1: Stick figure with wobbling masses

% define stick figure
st_x = [x_D, x_A, x_B, x_C];
st_y = [y_D, y_A, y_B, y_C];

% define the filename for the GIF
gif_filename = 'wobblingMass.gif';

% create a figure
fig = figure('Position', [100, 100, 800, 600]);

% loop through each time step
for ii = 1:length(x_A)
    % Wipe the slate clear
    clf
    hold on 
    
    % extract and plot data at the current time
    plot(st_x(ii, :), st_y(ii, :), 'k', 'LineWidth', 1.5)
    scatter(x_m_wob_t(ii), y_m_wob_t(ii), 60, 'filled', 'MarkerFaceColor', b2)
    scatter(x_m_wob_s(ii), y_m_wob_s(ii), 60, 'filled', 'MarkerFaceColor', b1)
    
    % set axis limits and describe plot
    axis([-0.5 0.5 0 1.0]);
    xlabel({'x [m]'}, 'FontSize', 20)
    ylabel({'y [m]'}, 'FontSize', 20)
    legend('Left leg', 'Wobbling mass thigh', 'Wobbling mass shank', 'FontSize', 15)
    
    % force MATLAB to draw the current image
    drawnow;

    % capture the current frame
    frame = getframe(fig);
    
    % convert the frame to an indexed image for GIF
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    
    % save the frame to the GIF file
    if ii == 1
        imwrite(imind, cm, gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.1);
    else
        imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
    end
end
