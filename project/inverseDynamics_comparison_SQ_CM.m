%% INFORMATION

%   This script can be used to calaculate the ankle, knee, and hip joint 
%   forces and torques with measured:
%   - body weight [kg]
%   - joint (toe (D), ankle (A), knee (B), hip (C)) coordinates [m]
%   - COP coordinates [m]
%   - ground reaction force [N]

%   It can also be used to calculate jump height, takeoff power,
%   acceleration and velocity with recorded GRF and model data.

%   In combination with the Simulink wobbling mass model, ankle, knee and 
%   hip joint forces can be calculated considering wobbling masses.

%% 
clear; close all; clc;

% *************************************************************************
%% ***************************SQUAT JUMP***********************************

%% -----------------------------INPUT--------------------------------------

%% Fix parameters 

% gravitational acceleration [m/s^2]
    g = 10;

% attachment point of spring mass damper 1 (SMD1) proximal from segment COM
    att_prox_COM = 0.3; 

%% Parameters changing with subject

% body mass [kg] 
    m = 70; 

% relative COM-position from proximal joint (based on de Leva 1996)
% male
    COM_rel_f = 0.4415;
    COM_rel_s = 0.4459;
    COM_rel_t = 0.4095;

% relative segment mass of one foot/shank/thigh 
% (based on Zatsiorsky 1990a, 1990b and 1993)
% male
    m_rel_f = 0.0137;
    m_rel_s = 0.0433;
    m_rel_t = 0.1416;

%% Parameter changing with measurement settings

% sampling rates [Hz]
    sr_Kinematics = 100;
    sr_ForcePlate = 1000;

%% Parameters depending on time                                  
    
% import model joint data (left leg) and divide by 1000 -> [m]
    x_C_SQ = xlsread('Gelenksdaten06','C:C')/1000;               %hip (HJC)
    y_C_SQ = xlsread('Gelenksdaten06','D:D')/1000;
    x_B_SQ = xlsread('Gelenksdaten06','E:E')/1000;               %knee (KJC)
    y_B_SQ = xlsread('Gelenksdaten06','F:F')/1000;
    x_A_SQ = xlsread('Gelenksdaten06','G:G')/1000;               %ankle (AJC)
    y_A_SQ = xlsread('Gelenksdaten06','H:H')/1000;
    x_D_SQ = xlsread('Gelenksdaten06','I:I')/1000;               %toe (TO)
    y_D_SQ = xlsread('Gelenksdaten06','J:J')/1000;
    
% import GRF [N] and divide by 2 (one leg) 
% import COP data and divide by 1000 -> [m]
    F_x_long_SQ = xlsread('GRFundCoP06','C:C')/2 ;
    F_y_long_SQ = xlsread('GRFundCoP06','D:D')/2 ;
    x_COP_long_SQ = xlsread('GRFundCoP06','E:E')/1000;
    y_COP_long_SQ = xlsread('GRFundCoP06','F:F')/1000;
    
% Since Vicon data sets m*g=0 shift F_y data by 0.5*m*g(one leg)
    F_y_long_SQ = F_y_long_SQ + 0.5*m*g;
    
% filter force and COP data 
    F_x_long_SQ = movmean(F_x_long_SQ,11);
    F_y_long_SQ = movmean(F_y_long_SQ, 11);
    x_COP_long_SQ = movmean(x_COP_long_SQ, 11);
    y_COP_long_SQ = movmean(y_COP_long_SQ, 11);
    
% preallocate empty variables
    F_x_SQ = zeros(length(x_A_SQ), 1);
    F_y_SQ = zeros(length(x_A_SQ), 1);
    x_COP_SQ = zeros(length(x_A_SQ), 1);
    y_COP_SQ = zeros(length(x_A_SQ), 1);
    
% reduce number of force and COP data points 
    sr_ratio_SQ = sr_ForcePlate/sr_Kinematics;
    for ii = 1 : length(x_A_SQ) 
        F_x_SQ(ii, 1) = F_x_long_SQ(sr_ratio_SQ * ii - (sr_ratio_SQ-1));
        F_y_SQ(ii, 1) = F_y_long_SQ(sr_ratio_SQ * ii - (sr_ratio_SQ-1));
        x_COP_SQ(ii, 1) = x_COP_long_SQ(sr_ratio_SQ * ii - (sr_ratio_SQ-1));
        y_COP_SQ(ii, 1) = y_COP_long_SQ(sr_ratio_SQ * ii - (sr_ratio_SQ-1));
    end
    
%--------------------------------------------------------------------------    
%% shorten data at the end, make vectors same length
    
% shorten standing phase at the end of data. y_C (hip) indicates that
% legs are extended (0.88 m).
        
    x_C_SQ = x_C_SQ(1:200);
    y_C_SQ = y_C_SQ(1:200);
    x_B_SQ = x_B_SQ(1:200);
    y_B_SQ = y_B_SQ(1:200);
    x_A_SQ = x_A_SQ(1:200);
    y_A_SQ = y_A_SQ(1:200);
    x_D_SQ = x_D_SQ(1:200);
    y_D_SQ = y_D_SQ(1:200);
  
% shorten GRF and COP data
    F_x_SQ = F_x_SQ(1:200);
    F_y_SQ = F_y_SQ(1:200);
    x_COP_SQ = x_COP_SQ(1:200);
    y_COP_SQ = y_COP_SQ(1:200);
        
%% Plot: GRF data full length and split into takeoff, flying and landing
    
%split GRF vector into takeoff, flying phase and landing
    fp_SQ = find(F_y_SQ <= 0.5*m*g);
    e_takeoff_SQ = fp_SQ(4);                                % takeoff: first data point after heels lift
    e_landing_SQ = fp_SQ(58);                               % landing: last data point before heels touch ground

    F_y_1_SQ = F_y_SQ(1:e_takeoff_SQ);                      % GRF until takeoff
    F_y_2_SQ = F_y_SQ(e_takeoff_SQ:e_landing_SQ);           % GRF in flight phase
    F_y_3_SQ = F_y_SQ(e_landing_SQ:end);                    % GRF from landing until end

    time_SQ = (1:1:length(F_y_SQ))';                        % time vector for F_y_SQ
    time_SQ = (time_SQ)* 0.01;                              % [s]
    
    time_1_SQ = time_SQ(1:e_takeoff_SQ);                    % time vector until takeoff
    time_2_SQ = time_SQ(e_takeoff_SQ:e_landing_SQ);         % time vector in flight phase
    time_3_SQ = time_SQ(e_landing_SQ:end);                  % time vector from landing until end
    
%Plot GRF data at full length and split into takeoff, flight phase and
%landing
    figure (1);
    sgtitle('GRF Squat Jump')
    
    subplot(2,3,1:3)
        plot(time_SQ,F_y_SQ, 'LineWidth', 5);
        axis([0.01 max(time_SQ) -100 1500]);grid on
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[ N ]'}, 'FontSize', 15)
        
    subplot(2,3,4);
        plot(time_1_SQ,F_y_1_SQ, 'LineWidth', 5); title('Takeoff')
        axis([0.01 time_SQ(e_takeoff_SQ) -100 1500]);grid on
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[ N ]'}, 'FontSize', 15)
        
    subplot(2,3,5);
        plot(time_2_SQ,F_y_2_SQ, 'LineWidth', 5);title('Flight Phase')
        axis([time_SQ(e_takeoff_SQ) time_SQ(e_landing_SQ) -100 1500]);grid on
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[ N ]'}, 'FontSize', 15)
        
    subplot(2,3,6);
        plot(time_3_SQ,F_y_3_SQ, 'LineWidth', 5);title('Landing')
        axis([time_SQ(e_landing_SQ) max(time_3_SQ) -100 1500]);grid on
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[ N ]'}, 'FontSize', 15)
    
%--------------------------------------------------------------------------    
%% ----------------------CALCULATIONS:SQUAT JUMP---------------------------

%% Calculate segment mass, COM coordinates and moments of inertia
    
% calculate segment mass of one foot (m_f), shank (m_s), and thigh (m_t) [kg]
    m_f = m_rel_f * m;
    m_s = m_rel_s * m;
    m_t = m_rel_t * m;

% calculate COM coordinates from measured joint coordinates [m]
    x_COM_f_SQ = x_A_SQ + COM_rel_f * (x_D_SQ - x_A_SQ);
    y_COM_f_SQ = y_A_SQ + COM_rel_f * (y_D_SQ - y_A_SQ);
    x_COM_s_SQ = x_B_SQ + COM_rel_s * (x_A_SQ - x_B_SQ);
    y_COM_s_SQ = y_B_SQ + COM_rel_s * (y_A_SQ - y_B_SQ);
    x_COM_t_SQ = x_C_SQ + COM_rel_t * (x_B_SQ - x_C_SQ);
    y_COM_t_SQ = y_C_SQ + COM_rel_t * (y_B_SQ - y_C_SQ);  
        
% calculate segment lengths from measured joint coordinates [m]
    l_f_SQ = sqrt( (x_A_SQ - x_D_SQ).^2 + (y_A_SQ - y_D_SQ).^2 );
    l_s_SQ = sqrt( (x_B_SQ - x_A_SQ).^2 + (y_B_SQ - y_A_SQ).^2 );
    l_t_SQ = sqrt( (x_C_SQ - x_B_SQ).^2 + (y_C_SQ - y_B_SQ).^2 );

% calculate moments of inertia [kg * m^2]
    J_f_SQ = (l_f_SQ * (COM_rel_f/2)).^2 * 0.5 * m_f + (l_f_SQ * ( (1- COM_rel_f)/2)).^2 * 0.5 * m_f;
    J_s_SQ = (l_s_SQ * (COM_rel_s/2)).^2 * 0.5 * m_s + (l_s_SQ * ( (1- COM_rel_s)/2)).^2 * 0.5 * m_s;
    J_t_SQ = (l_t_SQ * (COM_rel_t/2)).^2 * 0.5 * m_t + (l_t_SQ * ( (1- COM_rel_t)/2)).^2 * 0.5 * m_t;
      
%% Calculate COM and angular acceleration  

% calculate vertical and horizontal COM acceleration (numerical approximation)
    a_x_COM_f_SQ = diff( diff(x_COM_f_SQ)/(1/sr_Kinematics) )/(1/sr_Kinematics);
    a_y_COM_f_SQ = diff( diff(y_COM_f_SQ)/(1/sr_Kinematics) )/(1/sr_Kinematics);
    a_x_COM_s_SQ = diff( diff(x_COM_s_SQ)/(1/sr_Kinematics) )/(1/sr_Kinematics);
    a_y_COM_s_SQ = diff( diff(y_COM_s_SQ)/(1/sr_Kinematics) )/(1/sr_Kinematics);
    a_x_COM_t_SQ = diff( diff(x_COM_t_SQ)/(1/sr_Kinematics) )/(1/sr_Kinematics);
    a_y_COM_t_SQ = diff( diff(y_COM_t_SQ)/(1/sr_Kinematics) )/(1/sr_Kinematics);
    
% determine angles between segments
    DADA_SQ = sqrt((x_D_SQ-x_A_SQ).^2 + (y_D_SQ-y_A_SQ).^2);    %||D-A||
    ABAB_SQ = sqrt((x_A_SQ-x_B_SQ).^2 + (y_A_SQ-y_B_SQ).^2);    %||A-B||
    BCBC_SQ = sqrt((x_B_SQ-x_C_SQ).^2 + (y_B_SQ-y_C_SQ).^2);    %||B-C||

    phi_f_SQ = asin( (y_D_SQ - y_A_SQ)./  (DADA_SQ) );
    phi_s_SQ = asin( (x_A_SQ - x_B_SQ)./  (ABAB_SQ) );
    phi_t_SQ = asin( (x_B_SQ - x_C_SQ)./  (BCBC_SQ) );
    
% calculate angular acceleration
    alpha_f_SQ = diff( diff(phi_f_SQ)/(1/sr_Kinematics) )/(1/sr_Kinematics);
    alpha_s_SQ = diff( diff(phi_s_SQ)/(1/sr_Kinematics) )/(1/sr_Kinematics);
    alpha_t_SQ = diff( diff(phi_t_SQ)/(1/sr_Kinematics) )/(1/sr_Kinematics);
    
%% Calculate jump height

% by flight phase: data points between takeoff and landing
    flightphase_SQ = find(diff(F_y_long_SQ <= 0.5*m*g));     
    takeoff_SQ = flightphase_SQ(2); 
    landing_SQ = flightphase_SQ(3);

    flight_time_SQ = (landing_SQ - takeoff_SQ) * 0.001;     % [s]
    
    v_takeoff_SQ   = g*flight_time_SQ/2;                    % [m/s]
    
    jump_height_SQ = v_takeoff_SQ^2/(2*g);                  % [m]
    
% by impulse
    F_y_SQ_ks = F_y_SQ *2 -(m*g);                           % y-component of GRF for 2 legs shifted back to 0=m*g
    
    time_SQ = (1:1:length(F_y_SQ_ks))';                     % time vector for F_y_SQ_ks
    time_SQ = (time_SQ)* 0.01;                              % [s]
    
    Q_SQ = cumtrapz(time_SQ,(F_y_SQ_ks));                   % Integral of force over time for F_y_SQ_ks
    
    v_SQ = Q_SQ/m;                                          % calculation of velocity
    v_takeoff_ks_SQ = (max(v_SQ) + abs(min(v_SQ))) /2;      % max(v_SQ) = v(takeoff); min(v_SQ) = v(landing)
                                                            % since absolute values of velocity at takeoff and landing are equal, 
                                                            % calculate mean of both velocities                                                         
    jump_height_ks_SQ = ((v_takeoff_ks_SQ)^2)/(2*g);
    
    a_SQ = diff(v_SQ)./diff(time_SQ);                       % calculation of acceleration
    time_SQ_a = time_SQ(1:end-1);                           % time vector for acceleration

% by model data
    jump_height_md_SQ = max(y_C_SQ)-0.88;                   % hip

%% Calculate Take-off Power
    
    Takeoff_Power_SQ = (m*g*jump_height_md_SQ)/ (0.88 - min(y_C_SQ));   % TP = mgh/s (s: displacement between hip data for extended and flexed legs)

%% Calculate ankle, knee and hip joint forces [N] and torques [Nm] w/o wobbling mass

% preallocate empty variables
    FA_x_SQ = zeros(length(a_x_COM_f_SQ),1);        
    FA_y_SQ = zeros(length(a_x_COM_f_SQ), 1);
    FB_x_SQ = zeros(length(a_x_COM_f_SQ), 1);
    FB_y_SQ = zeros(length(a_x_COM_f_SQ), 1);
    FC_x_SQ = zeros(length(a_x_COM_f_SQ), 1);
    FC_y_SQ = zeros(length(a_x_COM_f_SQ), 1);
    MA_SQ = zeros(length(a_x_COM_f_SQ), 1);
    MB_SQ = zeros(length(a_x_COM_f_SQ), 1);
    MC_SQ = zeros(length(a_x_COM_f_SQ), 1);

% calculate joint forces and torques
for t = 1 : length(a_x_COM_f_SQ) 
    
    % ankle joint
    FA_x_SQ(t) = m_f * a_x_COM_f_SQ(t) - F_x_SQ(t);
    FA_y_SQ(t) = m_f * a_y_COM_f_SQ(t) + F_y_SQ(t) + m_f * g;
    MA_SQ(t) = J_f_SQ(t) * alpha_f_SQ(t) - ( (x_COP_SQ(t) - x_COM_f_SQ(t)) * F_y_SQ(t) - (y_COP_SQ(t) - y_COM_f_SQ(t)) * F_x_SQ(t) ) ...
        - ( (x_A_SQ(t) - x_COM_f_SQ(t)) * FA_y_SQ(t) - (y_A_SQ(t) - y_COM_f_SQ(t)) * FA_x_SQ(t) );

    % knee joint
    FB_x_SQ(t) = m_s * a_x_COM_s_SQ(t) + FA_x_SQ(t);
    FB_y_SQ(t) = m_s * a_y_COM_s_SQ(t) + FA_y_SQ(t) + m_s * g;
    MB_SQ(t)   = J_s_SQ(t) * alpha_s_SQ(t) + MA_SQ(t) - ( (x_A_SQ(t) - x_COM_s_SQ(t)) * (-FA_y_SQ(t)) - (y_A_SQ(t) - y_COM_s_SQ(t)) * (-FA_x_SQ(t)) ) ...
                 - ( (x_B_SQ(t) - x_COM_s_SQ(t)) * FB_y_SQ(t) - (y_B_SQ(t) - y_COM_s_SQ(t)) * FB_x_SQ(t) );
        
    % hip joint
    FC_x_SQ(t) = m_t * a_x_COM_t_SQ(t) + FB_x_SQ(t);
    FC_y_SQ(t) = m_t * a_y_COM_t_SQ(t) + FB_y_SQ(t) + m_t * g;
    MC_SQ(t)   = J_t_SQ(t) * alpha_t_SQ(t) + MB_SQ(t) - ( (x_B_SQ(t) - x_COM_t_SQ(t)) * (-FB_y_SQ(t)) - (y_B_SQ(t) - y_COM_t_SQ(t)) * (-FB_x_SQ(t)) ) ...
                 - ( (x_C_SQ(t) - x_COM_t_SQ(t)) * FC_y_SQ(t) - (y_C_SQ(t) - y_COM_t_SQ(t)) * FC_x_SQ(t) );

end
    
%% Split calculated joint-forces and torques into takeoff and landing
    
    fp_SQ = find(F_y_SQ <= 0.5*m*g); 
    e_takeoff_SQ = fp_SQ(4);
    e_landing_SQ = fp_SQ(58); 
        
%----------------------takeoff w/o wobbling mass---------------------------
    FA_x_tk_SQ = FA_x_SQ(1:e_takeoff_SQ);           % ankle
    FA_y_tk_SQ = FA_y_SQ(1:e_takeoff_SQ);
    MA_tk_SQ = MA_SQ(1:e_takeoff_SQ);
    
    FB_x_tk_SQ = FB_x_SQ(1:e_takeoff_SQ);           % knee
    FB_y_tk_SQ = FB_y_SQ(1:e_takeoff_SQ);
    MB_tk_SQ = MB_SQ(1:e_takeoff_SQ);
    
    FC_x_tk_SQ = FC_x_SQ(1:e_takeoff_SQ);           % hip
    FC_y_tk_SQ = FC_y_SQ(1:e_takeoff_SQ);
    MC_tk_SQ = MC_SQ(1:e_takeoff_SQ);
    
%----------------------landing w/o wobbling mass---------------------------
    FA_x_lg_SQ = FA_x_SQ(e_landing_SQ:end);         % ankle
    FA_y_lg_SQ = FA_y_SQ(e_landing_SQ:end);
    MA_lg_SQ = MA_SQ(e_landing_SQ:end);
    
    FB_x_lg_SQ = FB_x_SQ(e_landing_SQ:end);         % knee
    FB_y_lg_SQ = FB_y_SQ(e_landing_SQ:end);
    MB_lg_SQ = MB_SQ(e_landing_SQ:end);
    
    FC_x_lg_SQ = FC_x_SQ(e_landing_SQ:end);         % hip
    FC_y_lg_SQ = FC_y_SQ(e_landing_SQ:end);
    MC_lg_SQ = MC_SQ(e_landing_SQ:end);
    
%--------------------------------------------------------------------------    
%% -----------------Calculations for SJ with wobbling mass-----------------

% --------------------------Wobbling mass model---------------------------- 

%% Calculate attachment points

% calculate attachment points of SMD1 for shank and thigh
    x_att_s_SQ = x_COM_s_SQ + att_prox_COM * (x_B_SQ - x_A_SQ);
    y_att_s_SQ = x_COM_s_SQ + att_prox_COM * (y_B_SQ - y_A_SQ);
    x_att_t_SQ = x_COM_t_SQ + att_prox_COM * (x_C_SQ - x_B_SQ);
    y_att_t_SQ = x_COM_t_SQ + att_prox_COM * (y_C_SQ - y_B_SQ);
  

%% Calculate spring and damper forces of SMD1 and SMD2 for shank and thigh

% calculate simulation time and define time vector for from workspace block
    sim_time_SQ = length(x_A_SQ)/sr_Kinematics; 
    time_SQ = linspace(0, sim_time_SQ, length(x_A_SQ))';

% call WM model and calculate spring and damper forces for shank
    p.xy_mass_0 = [x_COM_s_SQ(1)-0.02 , y_COM_s_SQ(1)-0.02];
    p.xy_mass_dot_0 = [0 0]; 
    p.k = 2500;
    p.c = 200;
    p.m_wob = 1/2 * m_s; 
    p.xy_att1 = [time_SQ, x_att_s_SQ, y_att_s_SQ];
    p.xy_COM = [time_SQ, x_COM_s_SQ, y_COM_s_SQ];
    out = sim('../lecture/WobblingMassModel', [0 sim_time_SQ]);
    time_s_SQ= out.tout;

    FK1_s_SQ = interp1(time_s_SQ, out.FK1, time_SQ);
    FD1_s_SQ = interp1(time_s_SQ, out.FD1, time_SQ);
    FK2_s_SQ = interp1(time_s_SQ, out.FK2, time_SQ);
    FD2_s_SQ = interp1(time_s_SQ, out.FD2, time_SQ);
    xy_m_wob_s_SQ = interp1(time_s_SQ, out.xy_m_wob, time_SQ);

% call WM model and calculate spring and damper forces for thigh
    p.xy_mass_0 = [x_COM_t_SQ(1)-0.04, y_COM_t_SQ(1)-0.04];
    p.xy_mass_dot_0 = [0 0]; 
    p.m_wob = 2/3 * m_t;
    p.xy_att1 = [time_SQ, x_att_t_SQ, y_att_t_SQ];
    p.xy_COM = [time_SQ, x_COM_t_SQ, y_COM_t_SQ];
    out = sim ('../lecture/WobblingMassModel', [0 sim_time_SQ]);
    time_t_SQ= out.tout;

    FK1_t_SQ = interp1(time_t_SQ, out.FK1, time_SQ);
    FD1_t_SQ = interp1(time_t_SQ, out.FD1, time_SQ);
    FK2_t_SQ = interp1(time_t_SQ, out.FK2, time_SQ);
    FD2_t_SQ = interp1(time_t_SQ, out.FD2, time_SQ);
    xy_m_wob_t_SQ = interp1(time_t_SQ, out.xy_m_wob, time_SQ);
    
%% Calculate ankle, knee and hip joint forces [N] and torques [Nm] with wobbling mass

% separate x- and y-component of spring, damper force & wobbling-mass-postion vectors
% shank
    FK1_s_x_SQ = FK1_s_SQ(:,1);
    FK1_s_y_SQ = FK1_s_SQ(:,2);
    FD1_s_x_SQ = FD1_s_SQ(:,1);
    FD1_s_y_SQ = FD1_s_SQ(:,2);
    FK2_s_x_SQ = FK2_s_SQ(:,1);
    FK2_s_y_SQ = FK2_s_SQ(:,2);
    FD2_s_x_SQ = FD2_s_SQ(:,1);
    FD2_s_y_SQ = FD2_s_SQ(:,2);
    x_m_wob_s_SQ = xy_m_wob_s_SQ(:,1);
    y_m_wob_s_SQ = xy_m_wob_s_SQ(:,2);

% thigh
    FK1_t_x_SQ = FK1_t_SQ(:,1);
    FK1_t_y_SQ = FK1_t_SQ(:,2);
    FD1_t_x_SQ = FD1_t_SQ(:,1);
    FD1_t_y_SQ = FD1_t_SQ(:,2);
    FK2_t_x_SQ = FK2_t_SQ(:,1);
    FK2_t_y_SQ = FK2_t_SQ(:,2);
    FD2_t_x_SQ = FD2_t_SQ(:,1);
    FD2_t_y_SQ = FD2_t_SQ(:,2);
    x_m_wob_t_SQ = xy_m_wob_t_SQ(:,1);
    y_m_wob_t_SQ = xy_m_wob_t_SQ(:,2);

% preallocate empty variables
    FA_x_wob_SQ = zeros(length(a_x_COM_f_SQ), 1);      
    FA_y_wob_SQ = zeros(length(a_x_COM_f_SQ), 1);
    FB_x_wob_SQ = zeros(length(a_x_COM_f_SQ), 1);
    FB_y_wob_SQ = zeros(length(a_x_COM_f_SQ), 1);
    FC_x_wob_SQ = zeros(length(a_x_COM_f_SQ), 1);
    FC_y_wob_SQ = zeros(length(a_x_COM_f_SQ), 1);
    MA_wob_SQ = zeros(length(a_x_COM_f_SQ), 1);
    MB_wob_SQ = zeros(length(a_x_COM_f_SQ), 1);
    MC_wob_SQ = zeros(length(a_x_COM_f_SQ), 1);

% calculate joint forces and torques
for t = 1 : length(a_x_COM_f_SQ) 
    
    % ankle joint
    FA_x_wob_SQ(t) = m_f * a_x_COM_f_SQ(t) - F_x_SQ(t);
    FA_y_wob_SQ(t) = m_f * a_y_COM_f_SQ(t) + F_y_SQ(t) + m_f * g;
    MA_wob_SQ(t) = J_f_SQ(t) * alpha_f_SQ(t) - ( (x_COP_SQ(t) - x_COM_f_SQ(t)) * F_y_SQ(t) - (y_COP_SQ(t) - y_COM_f_SQ(t)) * F_x_SQ(t) ) ...
        - ( (x_A_SQ(t) - x_COM_f_SQ(t)) * FA_y_wob_SQ(t) - (y_A_SQ(t) - y_COM_f_SQ(t)) * FA_x_wob_SQ(t) );

    % knee joint
    FB_x_wob_SQ(t) = m_s * a_x_COM_s_SQ(t) + FA_x_wob_SQ(t) + FK1_s_x_SQ(t) + FD1_s_x_SQ(t) + FK2_s_x_SQ(t) + FD2_s_x_SQ(t);
    FB_y_wob_SQ(t) = m_s * a_y_COM_s_SQ(t) + FA_y_wob_SQ(t) + m_s * g + FK1_s_y_SQ(t) + FD1_s_y_SQ(t) + FK2_s_y_SQ(t) + FD2_s_y_SQ(t);
    MB_wob_SQ(t) = J_s_SQ(t) * alpha_s_SQ(t) + MA_wob_SQ(t) - ( (x_A_SQ(t) - x_COM_s_SQ(t)) * (-FA_y_wob_SQ(t)) - (y_A_SQ(t) - y_COM_s_SQ(t)) * (-FA_x_wob_SQ(t)) ) ...
        - ( (x_B_SQ(t) - x_COM_s_SQ(t)) * FB_y_wob_SQ(t) - (y_B_SQ(t) - y_COM_s_SQ(t)) * FB_x_wob_SQ(t) ) ...
        - ( (x_att_s_SQ(t) - x_COM_s_SQ(t)) * -FK1_s_y_SQ(t) - (y_att_s_SQ(t) - y_COM_s_SQ(t)) * -FK1_s_x_SQ(t) ) ... % spring
        - ( (x_att_s_SQ(t) - x_COM_s_SQ(t)) * -FD1_s_y_SQ(t) - (y_att_s_SQ(t) - y_COM_s_SQ(t)) * -FD1_s_x_SQ(t) );    % damper

    % hip joint
    FC_x_wob_SQ(t) = m_t * a_x_COM_t_SQ(t) + FB_x_wob_SQ(t) + FK1_t_x_SQ(t) + FD1_t_x_SQ(t) + FK2_t_x_SQ(t) + FD2_t_x_SQ(t);
    FC_y_wob_SQ(t) = m_t * a_y_COM_t_SQ(t) + FB_y_wob_SQ(t) + m_t * g + FK1_t_y_SQ(t) + FD1_t_y_SQ(t) + FK2_t_y_SQ(t) + FD2_t_y_SQ(t);
    MC_wob_SQ(t) = J_t_SQ(t) * alpha_t_SQ(t) + MB_wob_SQ(t) - ( (x_B_SQ(t) - x_COM_t_SQ(t)) * (-FB_y_wob_SQ(t)) - (y_B_SQ(t) - y_COM_t_SQ(t)) * (-FB_x_wob_SQ(t)) ) ...
        - ( (x_C_SQ(t) - x_COM_t_SQ(t)) * FC_y_wob_SQ(t) - (y_C_SQ(t) - y_COM_t_SQ(t)) * FC_x_wob_SQ(t) ) ...
        - ( (x_att_t_SQ(t) - x_COM_t_SQ(t)) * -FK1_t_y_SQ(t) - (y_att_t_SQ(t) - y_COM_t_SQ(t)) * -FK1_t_x_SQ(t) ) ... % spring
        - ( (x_att_t_SQ(t) - x_COM_t_SQ(t)) * -FD1_t_y_SQ(t) - (y_att_t_SQ(t) - y_COM_t_SQ(t)) * -FD1_t_x_SQ(t) );    % damper

end

%% Split calculated joint-forces and torques into takeoff and landing 
    
    fp_SQ = find(F_y_SQ <= 0.5*m*g); 
    e_takeoff_SQ = fp_SQ(4);
    e_landing_SQ = fp_SQ(58); 
    
%------------------------takeoff with wobbling mass------------------------
    FA_x_tk_wob_SQ = FA_x_wob_SQ(1:e_takeoff_SQ);       % ankle
    FA_y_tk_wob_SQ = FA_y_wob_SQ(1:e_takeoff_SQ);
    MA_tk_wob_SQ = MA_wob_SQ(1:e_takeoff_SQ);
    
    FB_x_tk_wob_SQ = FB_x_wob_SQ(1:e_takeoff_SQ);       % knee
    FB_y_tk_wob_SQ = FB_y_wob_SQ(1:e_takeoff_SQ);
    MB_tk_wob_SQ = MB_wob_SQ(1:e_takeoff_SQ);
    
    FC_x_tk_wob_SQ = FC_x_wob_SQ(1:e_takeoff_SQ);       % hip
    FC_y_tk_wob_SQ = FC_y_wob_SQ(1:e_takeoff_SQ);
    MC_tk_wob_SQ = MC_wob_SQ(1:e_takeoff_SQ);
    
%------------------------landing with wobbling mass------------------------
    FA_x_lg_wob_SQ = FA_x_wob_SQ(e_landing_SQ:end);     % ankle
    FA_y_lg_wob_SQ = FA_y_wob_SQ(e_landing_SQ:end);
    MA_lg_wob_SQ = MA_wob_SQ(e_landing_SQ:end);
    
    FB_x_lg_wob_SQ = FB_x_wob_SQ(e_landing_SQ:end);     % knee
    FB_y_lg_wob_SQ = FB_y_wob_SQ(e_landing_SQ:end);
    MB_lg_wob_SQ = MB_wob_SQ(e_landing_SQ:end);
    
    FC_x_lg_wob_SQ = FC_x_wob_SQ(e_landing_SQ:end);     % hip
    FC_y_lg_wob_SQ = FC_y_wob_SQ(e_landing_SQ:end);
    MC_lg_wob_SQ = MC_wob_SQ(e_landing_SQ:end);
 
% *************************************************************************
%%  OUTPUT - save joint forces and torques of SJ into structure

    % matrices with forces and torques for each segment with wobbling mass
    invdyn.SQ_wm_ankle_Fx_Fy_M = [FA_x_wob_SQ, FA_y_wob_SQ, MA_wob_SQ];
    invdyn.SQ_wm_knee_Fx_Fy_M = [FB_x_wob_SQ, FB_y_wob_SQ, MB_wob_SQ];
    invdyn.SQ_wm_hip_Fx_Fy_M = [FC_x_wob_SQ, FC_y_wob_SQ, MC_wob_SQ];  
    
    % matrices with forces and torques for each segment w/o wobbling mass
    invdyn.SQ_ankle_Fx_Fy_M = [FA_x_SQ, FA_y_SQ, MA_SQ];
    invdyn.SQ_knee_Fx_Fy_M = [FB_x_SQ, FB_y_SQ, MB_SQ];
    invdyn.SQ_hip_Fx_Fy_M = [FC_x_SQ, FC_y_SQ, MC_SQ];
    
    % save data in invdyn.mat 
    save invdynSquatJump.mat invdyn
% *************************************************************************
% ----------------------End of Calculations for SJ-------------------------

%% ------------------------------------------------------------------------
%% *************************COUNTERMOVEMENT JUMP***************************

%% -------------------------------INPUT------------------------------------

% Fix parameters, Parameters changing with subject and Parameter changing
% with measurement setting remain the same as in the Squat Jump.

% Parameters depending on time
    
% import model joint data (left legt) and divide by 1000 -> [m]
    x_C_CM = xlsread('Gelenksdaten08','C:C')/1000;               % hip (HJC)
    y_C_CM = xlsread('Gelenksdaten08','D:D')/1000;
    x_B_CM = xlsread('Gelenksdaten08','E:E')/1000;               % knee (KJC)
    y_B_CM = xlsread('Gelenksdaten08','F:F')/1000;
    x_A_CM = xlsread('Gelenksdaten08','G:G')/1000;               % ankle (AJC)
    y_A_CM = xlsread('Gelenksdaten08','H:H')/1000;
    x_D_CM = xlsread('Gelenksdaten08','I:I')/1000;               % toe (TO)
    y_D_CM = xlsread('Gelenksdaten08','J:J')/1000;
    
% import GRF [N] and divide by 2 (one leg) 
% import COP data and divide by 1000 -> [m]
    F_x_long_CM = xlsread('GRFundCoP08','C:C')/2 ;
    F_y_long_CM = xlsread('GRFundCoP08','D:D')/2 ;
    x_COP_long_CM = xlsread('GRFundCoP08','E:E')/1000;
    y_COP_long_CM = xlsread('GRFundCoP08','F:F')/1000;
    
% Since Vicon data sets m*g=0 shift F_y data by 0.5*m*g(one leg)
    F_y_long_CM = F_y_long_CM + 0.5*m*g;
    
% filter force and COP data 
    F_x_long_CM = movmean(F_x_long_CM,11);
    F_y_long_CM = movmean(F_y_long_CM, 11);
    x_COP_long_CM = movmean(x_COP_long_CM, 11);
    y_COP_long_CM = movmean(y_COP_long_CM, 11);
    
% preallocate empty variables
    F_x_CM = zeros(length(x_A_CM), 1);
    F_y_CM = zeros(length(x_A_CM), 1);
    x_COP_CM = zeros(length(x_A_CM), 1);
    y_COP_CM = zeros(length(x_A_CM), 1);
    
% reduce number of force and COP data points 
    sr_ratio_CM = sr_ForcePlate/sr_Kinematics;
    for ii = 1 : length(x_A_CM) 
        F_x_CM(ii, 1) = F_x_long_CM(sr_ratio_CM * ii - (sr_ratio_CM-1));
        F_y_CM(ii, 1) = F_y_long_CM(sr_ratio_CM * ii - (sr_ratio_CM-1));
        x_COP_CM(ii, 1) = x_COP_long_CM(sr_ratio_CM * ii - (sr_ratio_CM-1));
        y_COP_CM(ii, 1) = y_COP_long_CM(sr_ratio_CM * ii - (sr_ratio_CM-1));
    end
    
%------------------------------------------------------------------------    
%% shorten data at the end, make vectors same length
    
% shorten standing phase at the end of data. y_C (hip) indicates that
% legs are extended (0.88 m).
  
    x_C_CM = x_C_CM(1:325);
    y_C_CM = y_C_CM(1:325);
    x_B_CM = x_B_CM(1:325);
    y_B_CM = y_B_CM(1:325);
    x_A_CM = x_A_CM(1:325);
    y_A_CM = y_A_CM(1:325);
    x_D_CM = x_D_CM(1:325);
    y_D_CM = y_D_CM(1:325);
    
% shorten GRF and COP data
    F_x_CM = F_x_CM(1:325);
    F_y_CM = F_y_CM(1:325);
    x_COP_CM = x_COP_CM(1:325);
    y_COP_CM = y_COP_CM(1:325);
    
%% Plot: GRF data full length and split into takeoff, flying and landing
    
%split GRF vector into takeoff, flying phase and landing
    fp_CM = find(F_y_CM <= 0.5*m*g);                        
    e_takeoff_CM = fp_CM(34);                               % takeoff: first data point after heels lift
    e_landing_CM = fp_CM(90);                               % landing: last data point before heels touch ground

    F_y_1_CM = F_y_CM(1:e_takeoff_CM);                      % GRF until takeoff
    F_y_2_CM = F_y_CM(e_takeoff_CM:e_landing_CM);           % GRF in flight phase
    F_y_3_CM = F_y_CM(e_landing_CM:end);                    % GRF from landing until end
      
    time_CM = (1:1:length(F_y_CM))';                        % time vector for F_y_CM
    time_CM = (time_CM)*0.01;                               % [s]
    
    time_1_CM = time_CM(1:e_takeoff_CM);                    % time vector until takeoff
    time_2_CM = time_CM(e_takeoff_CM:e_landing_CM);         % time vector in flight phase
    time_3_CM = time_CM(e_landing_CM:end);                  % time vector from landing until end
    
%Plot GRF data at full length and split into takeoff, flight phase and
%landing
    figure (2);
    sgtitle('GRF Countermovement Jump')
    
    subplot(2,3,1:3)
        plot(time_CM,F_y_CM, 'LineWidth', 5);
        axis([0.01 max(time_CM) -100 1550]);grid on
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[ N ]'}, 'FontSize', 15)
        
    subplot(2,3,4);
        plot(time_1_CM,F_y_1_CM, 'LineWidth', 5);title('Takeoff')
        axis([0.01 time_CM(e_takeoff_CM) -100 1550]);grid on
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[ N ]'}, 'FontSize', 15)
        
    subplot(2,3,5);
        plot(time_2_CM,F_y_2_CM, 'LineWidth', 5);title('Flight Phase')
        axis([time_CM(e_takeoff_CM) time_CM(e_landing_CM) -100 1550]);grid on
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[ N ]'}, 'FontSize', 15)
        
    subplot(2,3,6);
        plot(time_3_CM,F_y_3_CM, 'LineWidth', 5);title('Landing')
        axis([time_CM(e_landing_CM) max(time_3_CM) -100 1550]);grid on
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[ N ]'}, 'FontSize', 15)
    
%--------------------------------------------------------------------------    
%% ------------------CALCULATIONS:COUNTERMOVEMENT JUMP---------------------

%% Calculate segment mass, COM coordinates and moments of inertia

% calculate segment mass of one foot (m_f), shank (m_s), and thigh (m_t) [kg]
    m_f = m_rel_f * m;
    m_s = m_rel_s * m;
    m_t = m_rel_t * m;

% calculate COM coordinates from measured joint coordinates [m]
    x_COM_f_CM = x_A_CM + COM_rel_f * (x_D_CM - x_A_CM);
    y_COM_f_CM = y_A_CM + COM_rel_f * (y_D_CM - y_A_CM);
    x_COM_s_CM = x_B_CM + COM_rel_s * (x_A_CM - x_B_CM);
    y_COM_s_CM = y_B_CM + COM_rel_s * (y_A_CM - y_B_CM);
    x_COM_t_CM = x_C_CM + COM_rel_t * (x_B_CM - x_C_CM);
    y_COM_t_CM = y_C_CM + COM_rel_t * (y_B_CM - y_C_CM);
    
% calculate segment lengths from measured joint coordinates [m]
    l_f_CM = sqrt( (x_A_CM - x_D_CM).^2 + (y_A_CM - y_D_CM).^2 );
    l_s_CM = sqrt( (x_B_CM - x_A_CM).^2 + (y_B_CM - y_A_CM).^2 );
    l_t_CM = sqrt( (x_C_CM - x_B_CM).^2 + (y_C_CM - y_B_CM).^2 );

% calculate moments of inertia [kg * m^2]
    J_f_CM_V = (l_f_CM * (COM_rel_f/2)).^2 * 0.5 * m_f + (l_f_CM * ( (1- COM_rel_f)/2)).^2 * 0.5 * m_f;
    J_s_CM_V = (l_s_CM * (COM_rel_s/2)).^2 * 0.5 * m_s + (l_s_CM * ( (1- COM_rel_s)/2)).^2 * 0.5 * m_s;
    J_t_CM_V = (l_t_CM * (COM_rel_t/2)).^2 * 0.5 * m_t + (l_t_CM * ( (1- COM_rel_t)/2)).^2 * 0.5 * m_t;
    
%% Calculate COM and angular acceleration  

% calculate vertical and horizontal COM acceleration (numerical approximation)
    a_x_COM_f_CM = diff( diff(x_COM_f_CM)/(1/sr_Kinematics) )/(1/sr_Kinematics);
    a_y_COM_f_CM = diff( diff(y_COM_f_CM)/(1/sr_Kinematics) )/(1/sr_Kinematics);
    a_x_COM_s_CM = diff( diff(x_COM_s_CM)/(1/sr_Kinematics) )/(1/sr_Kinematics);
    a_y_COM_s_CM = diff( diff(y_COM_s_CM)/(1/sr_Kinematics) )/(1/sr_Kinematics);
    a_x_COM_t_CM = diff( diff(x_COM_t_CM)/(1/sr_Kinematics) )/(1/sr_Kinematics);
    a_y_COM_t_CM = diff( diff(y_COM_t_CM)/(1/sr_Kinematics) )/(1/sr_Kinematics);

% determine angles between segments
    DADA_CM = sqrt((x_D_CM-x_A_CM).^2 + (y_D_CM-y_A_CM).^2);
    ABAB_CM = sqrt((x_A_CM-x_B_CM).^2 + (y_A_CM-y_B_CM).^2);
    BCBC_CM = sqrt((x_B_CM-x_C_CM).^2 + (y_B_CM-y_C_CM).^2);

    phi_f_CM = asin( (y_D_CM - y_A_CM)./  (DADA_CM) );
    phi_s_CM = asin( (x_A_CM - x_B_CM)./  (ABAB_CM) );
    phi_t_CM = asin( (x_B_CM - x_C_CM)./  (BCBC_CM) );

% calculate angular acceleration
    alpha_f_CM = diff( diff(phi_f_CM)/(1/sr_Kinematics) )/(1/sr_Kinematics);
    alpha_s_CM = diff( diff(phi_s_CM)/(1/sr_Kinematics) )/(1/sr_Kinematics);
    alpha_t_CM = diff( diff(phi_t_CM)/(1/sr_Kinematics) )/(1/sr_Kinematics);

%% Calculate jump height

% by flight phase: data points between takeoff and landing
    flightphase_CM = find(diff(F_y_long_CM <= 0.5*m*g));
    takeoff_CM = flightphase_CM(3); 
    landing_CM = flightphase_CM(4);
    
    flight_time_CM = (landing_CM - takeoff_CM) * 0.001;     % [s]   
    
    v_takeoff_CM = g*flight_time_CM/2;                      % [m/s] 
    
    jump_height_CM = v_takeoff_CM^2/(2*g);                  % [m]

% by impulse
    F_y_CM_ks = F_y_CM *2 -(m*g);                           % y-component of GRF for 2 legs shifted back to 0=m*g
        
    time_CM = (1:1:length(F_y_CM))';                        % time vector for F_y_CM_ks
    time_CM = (time_CM)*0.01;                               % [s]
    
    Q_CM = cumtrapz(time_CM,F_y_CM_ks);                     % Integral of force over time for F_y_CM_ks
    
    v_CM = Q_CM/m;                                          % calculation of velocity
    v_takeoff_ks_CM = (max(v_CM) + abs(min(v_CM))) /2;      % max(v_SQ) = v(takeoff); min(v_SQ) = v(landing)
                                                            % since absolute values of velocity at takeoff and landing are equal, 
                                                            % calculate mean of both velocities 
    jump_height_ks_CM = ((v_takeoff_ks_CM)^2)/(2*g);  
    
    a_CM = diff(v_CM)./diff(time_CM);                       % calculation of acceleration 
    time_CM_a = time_CM(1:end-1);                           % time vector for acceleration
    
% by model data
    jump_height_md_CM = max(y_C_CM)-y_C_CM(1);              % hip
    
    
%% Calculate takeoff power

    Takeoff_Power_CM = (m*g*jump_height_md_CM)/ (0.88 - min(y_C_CM));   % TP = mgh/s (s: displacement between hip data for extended and flexed legs)
   
%--------------------------------------------------------------------------    
%% Calculate ankle, knee and hip joint forces [N] and torques [Nm] w/o wobbling masses

% preallocate empty variables
    FA_x_CM = zeros(length(a_x_COM_f_CM),1);        
    FA_y_CM = zeros(length(a_x_COM_f_CM), 1);
    FB_x_CM = zeros(length(a_x_COM_f_CM), 1);
    FB_y_CM = zeros(length(a_x_COM_f_CM), 1);
    FC_x_CM = zeros(length(a_x_COM_f_CM), 1);
    FC_y_CM = zeros(length(a_x_COM_f_CM), 1);
    MA_CM = zeros(length(a_x_COM_f_CM), 1);
    MB_CM = zeros(length(a_x_COM_f_CM), 1);
    MC_CM = zeros(length(a_x_COM_f_CM), 1);
    
for t = 1 : length(a_x_COM_f_CM) 
    
    % ankle joint
    FA_x_CM(t) = m_f * a_x_COM_f_CM(t) - F_x_CM(t);
    FA_y_CM(t) = m_f * a_y_COM_f_CM(t) + F_y_CM(t) + m_f * g;
    MA_CM(t) = J_f_CM_V(t) * alpha_f_CM(t) - ( (x_COP_CM(t) - x_COM_f_CM(t)) * F_y_CM(t) - (y_COP_CM(t) - y_COM_f_CM(t)) * F_x_CM(t) ) ...
        - ( (x_A_CM(t) - x_COM_f_CM(t)) * FA_y_CM(t) - (y_A_CM(t) - y_COM_f_CM(t)) * FA_x_CM(t) );

    % knee joint
    FB_x_CM(t) = m_s * a_x_COM_s_CM(t) + FA_x_CM(t);
    FB_y_CM(t) = m_s * a_y_COM_s_CM(t) + FA_y_CM(t) + m_s * g;
    MB_CM(t)   = J_s_CM_V(t) * alpha_s_CM(t) + MA_CM(t) - ( (x_A_CM(t) - x_COM_s_CM(t)) * (-FA_y_CM(t)) - (y_A_CM(t) - y_COM_s_CM(t)) * (-FA_x_CM(t)) ) ...
                 - ( (x_B_CM(t) - x_COM_s_CM(t)) * FB_y_CM(t) - (y_B_CM(t) - y_COM_s_CM(t)) * FB_x_CM(t) );
        

    % hip joint
    FC_x_CM(t) = m_t * a_x_COM_t_CM(t) + FB_x_CM(t);
    FC_y_CM(t) = m_t * a_y_COM_t_CM(t) + FB_y_CM(t) + m_t * g;
    MC_CM(t)   = J_t_CM_V(t) * alpha_t_CM(t) + MB_CM(t) - ( (x_B_CM(t) - x_COM_t_CM(t)) * (-FB_y_CM(t)) - (y_B_CM(t) - y_COM_t_CM(t)) * (-FB_x_CM(t)) ) ...
                 - ( (x_C_CM(t) - x_COM_t_CM(t)) * FC_y_CM(t) - (y_C_CM(t) - y_COM_t_CM(t)) * FC_x_CM(t) );

end

%% Split calculated joint-forces and torques into takeoff and landing
        
    fp_CM = find(F_y_CM <= 0.5*m*g); 
    e_takeoff_CM = fp_CM(34);
    e_landing_CM = fp_CM(90);
    
%----------------------takeoff w/o wobbling mass---------------------------
    FA_x_tk_CM = FA_x_CM(1:e_takeoff_CM);           %ankle
    FA_y_tk_CM = FA_y_CM(1:e_takeoff_CM);
    MA_tk_CM = MA_CM(1:e_takeoff_CM);
    
    FB_x_tk_CM = FB_x_CM(1:e_takeoff_CM);           %knee
    FB_y_tk_CM = FB_y_CM(1:e_takeoff_CM);
    MB_tk_CM = MB_CM(1:e_takeoff_CM);
    
    FC_x_tk_CM = FC_x_CM(1:e_takeoff_CM);           %hip
    FC_y_tk_CM = FC_y_CM(1:e_takeoff_CM);
    MC_tk_CM = MC_CM(1:e_takeoff_CM);
    
%----------------------landing w/o wobbling mass---------------------------
    FA_x_lg_CM = FA_x_CM(e_landing_CM:end);         %ankle
    FA_y_lg_CM = FA_y_CM(e_landing_CM:end);
    MA_lg_CM = MA_CM(e_landing_CM:end);
    
    FB_x_lg_CM = FB_x_CM(e_landing_CM:end);         %knee
    FB_y_lg_CM = FB_y_CM(e_landing_CM:end);
    MB_lg_CM = MB_CM(e_landing_CM:end);
    
    FC_x_lg_CM = FC_x_CM(e_landing_CM:end);         %hip
    FC_y_lg_CM = FC_y_CM(e_landing_CM:end);
    MC_lg_CM = MC_CM(e_landing_CM:end);
    
%--------------------------------------------------------------------------
%% -----------------Calculations for CMJ with wobbling mass----------------

% --------------------------Wobbling mass model---------------------------- 
 
 %% Calculate attachment points

% calculate attachment points of SMD1 for shank and thigh
    x_att_s_CM = x_COM_s_CM + att_prox_COM * (x_B_CM - x_A_CM);
    y_att_s_CM = x_COM_s_CM + att_prox_COM * (y_B_CM - y_A_CM);
    x_att_t_CM = x_COM_t_CM + att_prox_COM * (x_C_CM - x_B_CM);
    y_att_t_CM = x_COM_t_CM + att_prox_COM * (y_C_CM - y_B_CM);

%% Calculate spring and damper forces of SMD1 and SMD2 for shank and thigh

% calculate simulation time and define time vector for from workspace block
    sim_time_CM = length(x_A_CM)/sr_Kinematics; 
    time_CM = linspace(0, sim_time_CM, length(x_A_CM))';

% call WM model and calculate spring and damper forces for shank
    p.xy_mass_0 = [x_COM_s_CM(1) - 0.02, y_COM_s_CM(1) - 0.02];
    p.xy_mass_dot_0 = [0 0]; 
    p.k = 2500;
    p.c = 200;
    p.m_wob = 1/2 * m_s; 
    p.xy_att1 = [time_CM, x_att_s_CM, y_att_s_CM];
    p.xy_COM = [time_CM, x_COM_s_CM, y_COM_s_CM];
    out = sim('../lecture/WobblingMassModel', [0 sim_time_CM]);
    time_s_CM= out.tout;

    FK1_s_CM = interp1(time_s_CM, out.FK1, time_CM);
    FD1_s_CM = interp1(time_s_CM, out.FD1, time_CM);
    FK2_s_CM = interp1(time_s_CM, out.FK2, time_CM);
    FD2_s_CM = interp1(time_s_CM, out.FD2, time_CM);
    xy_m_wob_s_CM = interp1(time_s_CM, out.xy_m_wob, time_CM);

% call WM model and calculate spring and damper forces for thigh
    p.xy_mass_0 = [x_COM_t_CM(1) - 0.04, y_COM_t_CM(1) - 0.04];
    p.xy_mass_dot_0 = [0 0]; 
    p.m_wob = 2/3 * m_t;
    p.xy_att1 = [time_CM, x_att_t_CM, y_att_t_CM];
    p.xy_COM = [time_CM, x_COM_t_CM, y_COM_t_CM];
    out = sim ('../lecture/WobblingMassModel', [0 sim_time_CM]);
    time_t_CM= out.tout;

    FK1_t_CM = interp1(time_t_CM, out.FK1, time_CM);
    FD1_t_CM = interp1(time_t_CM, out.FD1, time_CM);
    FK2_t_CM = interp1(time_t_CM, out.FK2, time_CM);
    FD2_t_CM = interp1(time_t_CM, out.FD2, time_CM);
    xy_m_wob_t_CM = interp1(time_t_CM, out.xy_m_wob, time_CM);
    
%% Calculate ankle, knee and hip joint forces [N] and torques [Nm] with wobbling masses

% separate x- and y-component of spring, damper force & wobbling-mass-postion vectors
% shank
    FK1_s_x_CM = FK1_s_CM(:,1);
    FK1_s_y_CM = FK1_s_CM(:,2);
    FD1_s_x_CM = FD1_s_CM(:,1);
    FD1_s_y_CM = FD1_s_CM(:,2);
    FK2_s_x_CM = FK2_s_CM(:,1);
    FK2_s_y_CM = FK2_s_CM(:,2);
    FD2_s_x_CM = FD2_s_CM(:,1);
    FD2_s_y_CM = FD2_s_CM(:,2);
    x_m_wob_s_CM = xy_m_wob_s_CM(:,1);
    y_m_wob_s_CM = xy_m_wob_s_CM(:,2);
    
% thigh
    FK1_t_x_CM = FK1_t_CM(:,1);
    FK1_t_y_CM = FK1_t_CM(:,2);
    FD1_t_x_CM = FD1_t_CM(:,1);
    FD1_t_y_CM = FD1_t_CM(:,2);
    FK2_t_x_CM = FK2_t_CM(:,1);
    FK2_t_y_CM = FK2_t_CM(:,2);
    FD2_t_x_CM = FD2_t_CM(:,1);
    FD2_t_y_CM = FD2_t_CM(:,2);
    x_m_wob_t_CM = xy_m_wob_t_CM(:,1);
    y_m_wob_t_CM = xy_m_wob_t_CM(:,2);

% preallocate empty variables
    FA_x_wob_CM = zeros(length(a_x_COM_f_CM),1);        
    FA_y_wob_CM = zeros(length(a_x_COM_f_CM), 1);
    FB_x_wob_CM = zeros(length(a_x_COM_f_CM), 1);
    FB_y_wob_CM = zeros(length(a_x_COM_f_CM), 1);
    FC_x_wob_CM = zeros(length(a_x_COM_f_CM), 1);
    FC_y_wob_CM = zeros(length(a_x_COM_f_CM), 1);
    MA_wob_CM = zeros(length(a_x_COM_f_CM), 1);
    MB_wob_CM = zeros(length(a_x_COM_f_CM), 1);
    MC_wob_CM = zeros(length(a_x_COM_f_CM), 1);

% calculate joint forces and torques
for t = 1 : length(a_x_COM_f_CM) 
    
    % ankle joint
    FA_x_wob_CM(t) = m_f * a_x_COM_f_CM(t) - F_x_CM(t);
    FA_y_wob_CM(t) = m_f * a_y_COM_f_CM(t) + F_y_CM(t) + m_f * g;
    MA_wob_CM(t) = J_f_CM_V(t) * alpha_f_CM(t) - ( (x_COP_CM(t) - x_COM_f_CM(t)) * F_y_CM(t) - (y_COP_CM(t) - y_COM_f_CM(t)) * F_x_CM(t) ) ...
        - ( (x_A_CM(t) - x_COM_f_CM(t)) * FA_y_wob_CM(t) - (y_A_CM(t) - y_COM_f_CM(t)) * FA_x_wob_CM(t) );

    % knee joint
    FB_x_wob_CM(t) = m_s * a_x_COM_s_CM(t) + FA_x_wob_CM(t) + FK1_s_x_CM(t) + FD1_s_x_CM(t) + FK2_s_x_CM(t) + FD2_s_x_CM(t);
    FB_y_wob_CM(t) = m_s * a_y_COM_s_CM(t) + FA_y_wob_CM(t) + m_s * g + FK1_s_y_CM(t) + FD1_s_y_CM(t) + FK2_s_y_CM(t) + FD2_s_y_CM(t);
    MB_wob_CM(t) = J_s_CM_V(t) * alpha_s_CM(t) + MA_wob_CM(t) - ( (x_A_CM(t) - x_COM_s_CM(t)) * (-FA_y_wob_CM(t)) - (y_A_CM(t) - y_COM_s_CM(t)) * (-FA_x_wob_CM(t)) ) ...
        - ( (x_B_CM(t) - x_COM_s_CM(t)) * FB_y_wob_CM(t) - (y_B_CM(t) - y_COM_s_CM(t)) * FB_x_wob_CM(t) ) ...
        - ( (x_att_s_CM(t) - x_COM_s_CM(t)) * -FK1_s_y_CM(t) - (y_att_s_CM(t) - y_COM_s_CM(t)) * -FK1_s_x_CM(t) ) ... % spring
        - ( (x_att_s_CM(t) - x_COM_s_CM(t)) * -FD1_s_y_CM(t) - (y_att_s_CM(t) - y_COM_s_CM(t)) * -FD1_s_x_CM(t) );    % damper

    % hip joint
    FC_x_wob_CM(t) = m_t * a_x_COM_t_CM(t) + FB_x_wob_CM(t) + FK1_t_x_CM(t) + FD1_t_x_CM(t) + FK2_t_x_CM(t) + FD2_t_x_CM(t);
    FC_y_wob_CM(t) = m_t * a_y_COM_t_CM(t) + FB_y_wob_CM(t) + m_t * g + FK1_t_y_CM(t) + FD1_t_y_CM(t) + FK2_t_y_CM(t) + FD2_t_y_CM(t);
    MC_wob_CM(t) = J_t_CM_V(t) * alpha_t_CM(t) + MB_wob_CM(t) - ( (x_B_CM(t) - x_COM_t_CM(t)) * (-FB_y_wob_CM(t)) - (y_B_CM(t) - y_COM_t_CM(t)) * (-FB_x_wob_CM(t)) ) ...
        - ( (x_C_CM(t) - x_COM_t_CM(t)) * FC_y_wob_CM(t) - (y_C_CM(t) - y_COM_t_CM(t)) * FC_x_wob_CM(t) ) ...
        - ( (x_att_t_CM(t) - x_COM_t_CM(t)) * -FK1_t_y_CM(t) - (y_att_t_CM(t) - y_COM_t_CM(t)) * -FK1_t_x_CM(t) ) ... % spring
        - ( (x_att_t_CM(t) - x_COM_t_CM(t)) * -FD1_t_y_CM(t) - (y_att_t_CM(t) - y_COM_t_CM(t)) * -FD1_t_x_CM(t) );    % damper

end

%% Split calculated joint-forces and torques into takeoff and landing
        
    fp_CM = find(F_y_CM <= 0.5*m*g); 
    e_takeoff_CM = fp_CM(34);
    e_landing_CM = fp_CM(90);
    
%----------------------takeoff with wobbling mass--------------------------
    FA_x_tk_wob_CM = FA_x_wob_CM(1:e_takeoff_CM);       % ankle
    FA_y_tk_wob_CM = FA_y_wob_CM(1:e_takeoff_CM);
    MA_tk_wob_CM = MA_wob_CM(1:e_takeoff_CM);
    
    FB_x_tk_wob_CM = FB_x_wob_CM(1:e_takeoff_CM);       % knee
    FB_y_tk_wob_CM = FB_y_wob_CM(1:e_takeoff_CM);
    MB_tk_wob_CM = MB_wob_CM(1:e_takeoff_CM);
    
    FC_x_tk_wob_CM = FC_x_wob_CM(1:e_takeoff_CM);       % hip
    FC_y_tk_wob_CM = FC_y_wob_CM(1:e_takeoff_CM);
    MC_tk_wob_CM = MC_wob_CM(1:e_takeoff_CM);
    
%----------------------landing with wobbling mass--------------------------
    FA_x_lg_wob_CM = FA_x_wob_CM(e_landing_CM:end);     % ankle
    FA_y_lg_wob_CM = FA_y_wob_CM(e_landing_CM:end);
    MA_lg_wob_CM = MA_wob_CM(e_landing_CM:end);
    
    FB_x_lg_wob_CM = FB_x_wob_CM(e_landing_CM:end);     % knee
    FB_y_lg_wob_CM = FB_y_wob_CM(e_landing_CM:end);
    MB_lg_wob_CM = MB_wob_CM(e_landing_CM:end);
    
    FC_x_lg_wob_CM = FC_x_wob_CM(e_landing_CM:end);     % hip
    FC_y_lg_wob_CM = FC_y_wob_CM(e_landing_CM:end);
    MC_lg_wob_CM = MC_wob_CM(e_landing_CM:end);
    
% *************************************************************************
%%  OUTPUT - Save joint forces and torques of CMJ in a structure

    % matrices with forces and torques for each segment with wobbling mass
    invdyn.CM_wm_ankle_Fx_Fy_M = [FA_x_wob_CM, FA_y_wob_CM, MA_wob_CM];
    invdyn.CM_wm_knee_Fx_Fy_M = [FB_x_wob_CM, FB_y_wob_CM, MB_wob_CM];
    invdyn.CM_wm_hip_Fx_Fy_M = [FC_x_wob_CM, FC_y_wob_CM, MC_wob_CM];  
    
    % matrices with forces and torques for each segment w/o wobbling mass
    invdyn.CM_ankle_Fx_Fy_M = [FA_x_CM, FA_y_CM, MA_CM];
    invdyn.CM_knee_Fx_Fy_M = [FB_x_CM, FB_y_CM, MB_CM];
    invdyn.CM_hip_Fx_Fy_M = [FC_x_CM, FC_y_CM, MC_CM];
    
    % save data in invdyn.mat 
    save invdynCountermovementJump.mat invdyn
% *************************************************************************
% ----------------------End of Calculations for SJ-------------------------

%% Create time vectors for plots

    time_SQ = time_SQ(1:length(FA_y_SQ));
    time_SQ_tk = time_SQ(1:e_takeoff_SQ);
    time_SQ_lg = time_SQ(e_landing_SQ:length(a_x_COM_f_SQ));
    
    time_CM = time_CM(1:length(FA_y_CM));
    time_CM_tk = time_CM(1:e_takeoff_CM);
    time_CM_lg = time_CM(e_landing_CM:length(a_x_COM_f_CM));
    
%% --------------------Plots to compare CMJ and SJ-------------------------

% Plot 1: forces in joints (w/o wobbling mass) for SJ CMJ
    figure(3)
    sgtitle('Joint Forces w/o wobbling mass - full jump length', 'FontSize', 20)
    
    subplot(3,2,1); plot(time_SQ,FC_x_SQ,'g', 'LineWidth', 3); title('Hip - SJ', 'FontSize', 15)
        grid on
        hold on;    
        plot(time_SQ,FC_y_SQ,'b', 'LineWidth', 3)
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        axis([0 max(time_SQ) -1500 2500]);
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[N]'}, 'FontSize', 15)
        legend('x-component', 'y-component', 'FontSize', 11)
        
    subplot(3,2,2); plot(time_CM,FC_x_CM,'g', 'LineWidth', 3); title('Hip - CMJ', 'FontSize', 15)
        grid on
        hold on;    
        plot(time_CM,FC_y_CM,'b', 'LineWidth', 3)
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        axis([0 max(time_CM) -1500 2500]);
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[N]'}, 'FontSize', 15)
        legend('x-component', 'y-component', 'FontSize', 11)
    
    subplot(3,2,3); plot(time_SQ,FB_x_SQ,'g', 'LineWidth', 3); title('Knee - SJ', 'FontSize', 15)
        grid on
        hold on;    
        plot(time_SQ,FB_y_SQ,'b', 'LineWidth', 3)
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        axis([0 max(time_SQ) -1500 2500]);
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[N]'}, 'FontSize', 15)
        legend('x-component', 'y-component', 'FontSize', 11)
    
    subplot(3,2,4); plot(time_CM,FB_x_CM,'g', 'LineWidth', 3); title('Knee - CMJ', 'FontSize', 15)
        grid on
        hold on;    
        plot(time_CM,FB_y_CM,'b', 'LineWidth', 3)
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        axis([0 max(time_CM) -1500 2500]);
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[N]'}, 'FontSize', 15)
        legend('x-component', 'y-component', 'FontSize', 11)
    
    subplot(3,2,5); plot(time_SQ,FA_x_SQ,'g', 'LineWidth', 3); title('Ankle - SJ', 'FontSize', 15)
        grid on
        hold on;    
        plot(time_SQ,FA_y_SQ,'b', 'LineWidth', 3)
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        axis([0 max(time_SQ) -1500 2500]);
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[N]'}, 'FontSize', 15)
        legend('x-component', 'y-component', 'FontSize', 11)
    
    subplot(3,2,6); plot(time_CM,FA_x_CM,'g', 'LineWidth', 3); title('Ankle - CMJ', 'FontSize', 15)
        grid on
        hold on;    
        plot(time_CM,FA_y_CM,'b', 'LineWidth', 3)
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        axis([0 max(time_CM) -1500 2500]);
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[N]'}, 'FontSize', 15)
        legend('x-component', 'y-component', 'FontSize', 11)
    
%% Plot2: forces and moments in takeoff (w/o wobbling mass)
% Compare y_components of joint forces 

    figure(4)
    sgtitle('y-components of Forces in Joints for Takeoff', 'FontSize', 20)
    
    subplot(3,2,1); plot(time_SQ_tk,FC_y_tk_SQ, 'LineWidth', 3); title('Hip - SJ', 'FontSize', 15)
        grid on
        hold on
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        axis([0.01 max(time_SQ_tk)+0.01 100 1200]);
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[N]'}, 'FontSize', 15)
    
    subplot(3,2,2); plot(time_CM_tk,FC_y_tk_CM, 'LineWidth', 3); title('Hip - CMJ', 'FontSize', 15)
        grid on
        hold on
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        axis([0.01 max(time_CM_tk)+0.01 100 1200]);
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[N]'}, 'FontSize', 15)
        
    subplot(3,2,3); plot(time_SQ_tk,FB_y_tk_SQ, 'LineWidth', 3); title('Knee - SJ', 'FontSize', 15)
        grid on
        hold on
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        axis([0.01 max(time_SQ_tk)+0.01 100 1200]);
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[N]'}, 'FontSize', 15)
        
    subplot(3,2,4); plot(time_CM_tk,FB_y_tk_CM, 'LineWidth', 3); title('Knee - CMJ', 'FontSize', 15)
        grid on
        hold on
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        axis([0.01 max(time_CM_tk)+0.01 100 1200]);
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[N]'}, 'FontSize', 15)
        
    subplot(3,2,5); plot(time_SQ_tk,FA_y_tk_SQ, 'LineWidth', 3); title('Ankle - SJ', 'FontSize', 15)
        grid on
        hold on
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        axis([0.01 max(time_SQ_tk)+0.01 100 1200]);
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[N]'}, 'FontSize', 15)
        
    subplot(3,2,6); plot(time_CM_tk,FA_y_tk_CM, 'LineWidth', 3); title('Ankle - CMJ', 'FontSize', 15)
        grid on
        hold on
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        axis([0.01 max(time_CM_tk)+0.01 100 1200]);
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[N]'}, 'FontSize', 15)

% Compare moments in joints of SJ and CMJ

    figure(5)
    sgtitle('Moments in Joints for Takeoff', 'FontSize', 20)
    
    subplot(3,2,1); plot(time_SQ_tk,MC_tk_SQ); title('Hip - SJ', 'FontSize', 15)
        grid on
        hold on
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        axis([0.01 max(time_SQ_tk) -1000 1000]);
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[Nm]'}, 'FontSize', 15)
        
    subplot(3,2,2); plot(time_CM_tk,MC_tk_CM); title('Hip - CMJ', 'FontSize', 15)
        grid on
        hold on
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        axis([0 max(time_CM_tk) -1000 1000]);
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[Nm]'}, 'FontSize', 15)
        
    subplot(3,2,3); plot(time_SQ_tk,MB_tk_SQ); title('Knee - SJ', 'FontSize', 15)
        grid on
        hold on
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        axis([0.01 max(time_SQ_tk) -1000 1000]);
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[Nm]'}, 'FontSize', 15)
        
    subplot(3,2,4); plot(time_CM_tk,MB_tk_CM); title('Knee - CMJ', 'FontSize', 15)
        grid on
        hold on
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        axis([0 max(time_CM_tk) -1000 1000]);
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[Nm]'}, 'FontSize', 15)
        
    subplot(3,2,5); plot(time_SQ_tk,MA_tk_SQ); title('Ankle - SJ', 'FontSize', 15)
        grid on
        hold on
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        axis([0.01 max(time_SQ_tk) -1000 1000]);
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[Nm]'}, 'FontSize', 15)
    
    subplot(3,2,6); plot(time_CM_tk,MA_tk_CM); title('Ankle - CMJ', 'FontSize', 15)
        grid on
        hold on
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        axis([0 max(time_CM_tk) -1000 1000]);
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[Nm]'}, 'FontSize', 15)

%% Plot: forces and moments in landing (without wobbling mass)

% Compare forces in joints while landing for SJ and CMJ
    
    figure(6)
    sgtitle('y-components of Forces in Joints for Landing', 'FontSize', 20)
    
    subplot(3,2,1); plot(time_SQ_lg,FC_y_lg_SQ, 'LineWidth', 3); title('Hip - SJ', 'FontSize', 15)
        grid on
        hold on
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        axis([min(time_SQ_lg) max(time_SQ_lg) 0 2500]);
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[N]'}, 'FontSize', 15)
        
    subplot(3,2,2); plot(time_CM_lg,FC_y_lg_CM, 'LineWidth', 3); title('Hip - CMJ', 'FontSize', 15)
        grid on
        hold on
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        axis([min(time_CM_lg) max(time_CM_lg) 0 2500]);
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[N]'}, 'FontSize', 15)
        
    subplot(3,2,3); plot(time_SQ_lg,FB_y_lg_SQ, 'LineWidth', 3); title('Knee - SJ', 'FontSize', 15)
        grid on
        hold on
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        axis([min(time_SQ_lg) max(time_SQ_lg) 0 2500]);
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[N]'}, 'FontSize', 15)
        
    subplot(3,2,4); plot(time_CM_lg,FB_y_lg_CM, 'LineWidth', 3); title('Knee - CMJ', 'FontSize', 15)
        grid on
        hold on
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        axis([min(time_CM_lg) max(time_CM_lg) 0 2500]);
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[N]'}, 'FontSize', 15)
        
    subplot(3,2,5); plot(time_SQ_lg,FA_y_lg_SQ, 'LineWidth', 3); title('Ankle - SJ', 'FontSize', 15)
        grid on
        hold on
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        axis([min(time_SQ_lg) max(time_SQ_lg) 0 2500]);
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[N]'}, 'FontSize', 15)
        
    subplot(3,2,6); plot(time_CM_lg,FA_y_lg_CM, 'LineWidth', 3); title('Ankle - CMJ', 'FontSize', 15)
        grid on
        hold on
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        axis([min(time_CM_lg) max(time_CM_lg) 0 2500]);
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[N]'}, 'FontSize', 15)

% Compare moments in joints for landing (without wobbling mass)
    
    figure(7)
    sgtitle('Moments in Joints for Landing', 'FontSize', 20)
    
    subplot(3,2,1); plot(time_SQ_lg,MC_lg_SQ); title('Hip - SJ', 'FontSize', 15)
        grid on
        hold on
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        axis([min(time_SQ_lg) max(time_SQ_lg) -900 1000]);
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[Nm]'}, 'FontSize', 15)
        
    subplot(3,2,2); plot(time_CM_lg,MC_lg_CM); title('Hip - CMJ', 'FontSize', 15)
        grid on
        hold on
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        axis([min(time_CM_lg) max(time_CM_lg) -900 1000]);
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[Nm]'}, 'FontSize', 15)
        
    subplot(3,2,3); plot(time_SQ_lg,MB_lg_SQ); title('Knee - SJ', 'FontSize', 15)
        grid on
        hold on
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        axis([min(time_SQ_lg) max(time_SQ_lg) -900 1000]);
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[Nm]'}, 'FontSize', 15)
        
    subplot(3,2,4); plot(time_CM_lg,MB_lg_CM); title('Knee - CMJ', 'FontSize', 15)
        grid on
        hold on
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        axis([min(time_CM_lg) max(time_CM_lg) -900 1000]);
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[Nm]'}, 'FontSize', 15)
        
    subplot(3,2,5); plot(time_SQ_lg,MA_lg_SQ); title('Ankle - SJ', 'FontSize', 15)
        grid on
        hold on
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        axis([min(time_SQ_lg) max(time_SQ_lg) -900 1000]);
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[Nm]'}, 'FontSize', 15)
        
    subplot(3,2,6); plot(time_CM_lg,MA_lg_CM); title('Ankle - CMJ', 'FontSize', 15)
        grid on
        hold on
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        axis([min(time_CM_lg) max(time_CM_lg) -900 1000]);
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[Nm]'}, 'FontSize', 15)

%% Plot GRF, acceleration and velocity 
% *************************************************************************
% Calculation of key points
% key points for CMJ
    a_x_CM = 60;                                        % a: start of jump
    [b_y_CM, b_x_CM] = min(a_CM(1:e_takeoff_CM-10));    % b: maximal acceleration downwards
    [c_y_CM, c_x_CM] = min(v_CM(1:e_takeoff_CM-10));    % c: maximal velocity downwards
    d_y_CM = max(v_CM)-v_takeoff_ks_CM;                 % d: lowest point of jump 
    d_x_CM = 131;
    [e_y_CM, e_x_CM] = max(v_CM);                       % e: takeoff
    f_CM = find(F_y_CM <= 0);                           % f: start of flight (toe-off)
    f_x_CM = f_CM(1);
    [g_y_CM, g_x_CM] = min(v_CM);                       % g: landing
   
% key points for SJ
    a_x_SQ = 9;
    [e_y_SQ, e_x_SQ] = max(v_SQ);
    f_SQ = find(F_y_SQ <= 0);
    f_x_SQ = f_SQ(1);
    [g_y_SQ, g_x_SQ] = min(v_SQ);

% *************************************************************************
% ---------------------------End of Calculations---------------------------

% y-component of GRF forces shorted at the end by 2 due to 2 derivation
    F_y_SQ_plot = F_y_SQ(1:198)*2;              
    F_y_CM_plot = F_y_CM(1:323)*2;
    
    figure(8)
    sgtitle('GRF & Kinematic Curves', 'FontSize', 20)
    
    subplot(3,2,1); plot(time_SQ,F_y_SQ_plot, 'LineWidth', 2); title('GRF - SJ', 'FontSize', 15)
        hold on; 
        plot(a_x_SQ*0.01,F_y_SQ_plot(a_x_SQ),'r.', 'MarkerSize', 35)
        plot((e_x_SQ)*0.01,F_y_SQ_plot(e_x_SQ),'b.', 'MarkerSize', 35)
        plot(f_x_SQ*0.01,F_y_SQ_plot(f_x_SQ),'c.', 'MarkerSize', 35)
        plot((g_x_SQ-3)*0.01,F_y_SQ_plot(g_x_SQ-3),'y.', 'MarkerSize', 35)
        legend('GRF','a','e','f','g', 'FontSize',10,'Location','northwest')
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        axis([0 max(time_SQ) -100 3000]);grid on
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[N]'}, 'FontSize', 15)
        
    subplot(3,2,2); plot(time_CM,F_y_CM_plot, 'LineWidth', 2); title('GRF - CMJ', 'FontSize', 15)
        hold on; 
        plot(a_x_CM*0.01,F_y_CM_plot(a_x_CM),'r.', 'MarkerSize', 35)
        plot(b_x_CM*0.01,F_y_CM_plot(b_x_CM),'k.', 'MarkerSize', 35)
        plot(c_x_CM*0.01,F_y_CM_plot(c_x_CM),'m.', 'MarkerSize', 35)
        plot(d_x_CM*0.01,F_y_CM_plot(d_x_CM),'g.', 'MarkerSize', 35)
        plot(e_x_CM*0.01,F_y_CM_plot(e_x_CM),'b.', 'MarkerSize', 35)
        plot(f_x_CM*0.01,F_y_CM_plot(f_x_CM),'c.', 'MarkerSize', 35)
        plot((g_x_CM-3)*0.01,F_y_CM_plot(g_x_CM-3),'y.', 'MarkerSize', 35)
        legend('GRF','a', 'b', 'c','d','e','f','g', 'FontSize',10,'Location','northwest')
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        axis([0 max(time_CM) -100 3000]);grid on
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[N]'}, 'FontSize', 15)
        
    subplot(3,2,3); plot(time_SQ_a,a_SQ, 'LineWidth', 2); title('Acceleration - SJ','FontSize', 15)
        hold on; 
        plot(a_x_SQ*0.01,a_SQ(a_x_SQ),'r.', 'MarkerSize', 35)
        plot((e_x_SQ)*0.01,a_SQ(e_x_SQ),'b.', 'MarkerSize', 35)
        plot(f_x_SQ*0.01,a_SQ(f_x_SQ),'c.', 'MarkerSize', 35)
        plot((g_x_SQ-3)*0.01,a_SQ(g_x_SQ-3),'y.', 'MarkerSize', 35)
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        axis([0 max(time_SQ) -20 40]);grid on
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[ m/s² ]'}, 'FontSize', 15)
        
    subplot(3,2,4); plot(time_CM_a,a_CM, 'LineWidth', 2); title('Acceleration - CMJ', 'FontSize', 15)
        hold on; 
        plot(a_x_CM*0.01,a_CM(a_x_CM),'r.', 'MarkerSize', 35)
        plot(b_x_CM*0.01,a_CM(b_x_CM),'k.', 'MarkerSize', 35)
        plot(c_x_CM*0.01,a_CM(c_x_CM),'m.', 'MarkerSize', 35)
        plot(d_x_CM*0.01,a_CM(d_x_CM),'g.', 'MarkerSize', 35)
        plot(e_x_CM*0.01,a_CM(e_x_CM),'b.', 'MarkerSize', 35)
        plot(f_x_CM*0.01,a_CM(f_x_CM),'c.', 'MarkerSize', 35)
        plot((g_x_CM-3)*0.01,a_CM(g_x_CM-3),'y.', 'MarkerSize', 35)
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        axis([0 max(time_CM) -20 40]);grid on
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[ m/s² ]'}, 'FontSize', 15)
        
    subplot(3,2,5); plot(time_SQ,v_SQ(1:length(time_SQ)), 'LineWidth', 2); title('Velocity - SJ', 'FontSize', 15)
        hold on; 
        plot(a_x_SQ*0.01,v_SQ(a_x_SQ),'r.', 'MarkerSize', 35)
        plot((e_x_SQ-1)*0.01,v_SQ(e_x_SQ),'b.', 'MarkerSize', 35)
        plot(f_x_SQ*0.01,v_SQ(f_x_SQ),'c.', 'MarkerSize', 35)
        plot(g_x_SQ*0.01,v_SQ(g_x_SQ),'y.', 'MarkerSize', 35)
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        axis([0 max(time_SQ) -3 4]);grid on
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[ m/s ]'}, 'FontSize', 15)
        
    subplot(3,2,6); plot(time_CM,v_CM(1:length(time_CM)), 'LineWidth', 2); title('Velocity - CMJ', 'FontSize', 15)
        hold on; 
        plot(a_x_CM*0.01,v_CM(a_x_CM),'r.', 'MarkerSize', 35)  
        plot(b_x_CM*0.01,v_CM(b_x_CM),'k.', 'MarkerSize', 35)
        plot(c_x_CM*0.01,v_CM(c_x_CM),'m.', 'MarkerSize', 35)
        plot(d_x_CM*0.01,v_CM(d_x_CM),'g.', 'MarkerSize', 35)
        plot(e_x_CM*0.01,v_CM(e_x_CM),'b.', 'MarkerSize', 35)
        plot(f_x_CM*0.01,v_CM(f_x_CM),'c.', 'MarkerSize', 35)
        plot(g_x_CM*0.01,v_CM(g_x_CM),'y.', 'MarkerSize', 35)
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        axis([0 max(time_CM) -3 4]);grid on
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[ m/s ]'}, 'FontSize', 15)
  
%% Plot jump height

% highest point in hip model data
    [y_SQ, x_SQ] = max(y_C_SQ);
    [y_CM, x_CM] = max(y_C_CM);
    
    figure(9)
    sgtitle('Jump Height', 'FontSize', 20)
    
    subplot(1,2,1); plot(time_SQ,y_C_SQ(1:length(time_SQ))-0.88,'b', 'LineWidth', 3); title('Jump Height - SJ', 'FontSize', 15)
        grid on
        hold on
        axis([0 max(time_SQ) -0.5 0.5]);
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[m]'}, 'FontSize', 15)
        plot(x_SQ/100,jump_height_SQ,'rx','MarkerSize', 25, 'LineWidth', 3)
        plot(x_SQ/100,jump_height_ks_SQ,'kx','MarkerSize', 25, 'LineWidth', 3)
        legend('by Hip Marker: 0.4370 m', 'by Flying Time: 0.3754 m', 'by Impulse: 0.3533 m', 'FontSize', 23) 
        
    subplot(1,2,2); plot(time_CM,y_C_CM(1:length(time_CM))-0.88,'b', 'LineWidth', 3); title('Jump Height - CMJ', 'FontSize', 15) %-0.88(marker) -0.04(Platte)
        grid on
        hold on
        axis([0 max(time_CM) -0.5 0.5]);
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[m]'}, 'FontSize', 15)
        plot(x_CM/100,jump_height_CM,'rx','MarkerSize', 25, 'LineWidth', 3)
        plot(x_CM/100,jump_height_ks_CM,'kx','MarkerSize', 25, 'LineWidth', 3)
        legend('by Hip Marker: 0.4495 m', 'by Flying Time: 0.4104 m' , 'by Impulse: 0.3866 m', 'FontSize', 23,'Location','northwest')
        
%% Plot: Compare y-components of Joint forces with and without wobbling mass     
 
% y-components of Joint forces in hip & knee joint with and without wobbling mass - SJ
    figure(10)
    sgtitle('Joint Forces SJ with & without wobbling mass', 'FontSize', 20)
 
    subplot(1,2,1); plot(time_SQ,FC_y_SQ +(m*g),'r','LineWidth',1); title('Hip - SJ', 'FontSize', 15)
        hold on;  
        plot(time_SQ,FC_y_wob_SQ +(m*g),'k','LineWidth',1)
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        axis([0 max(time_SQ) 0 3500]);grid on
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[N]'}, 'FontSize', 15)
        legend('without wobbling mass', 'with wobbling mass', 'FontSize', 13,'Location','northwest')
        
    subplot(1,2,2); plot(time_SQ,FB_y_SQ +(m*g),'r','LineWidth',1); title('Knee - SJ', 'FontSize', 15)
        hold on;  
        plot(time_SQ,FB_y_wob_SQ +(m*g),'k','LineWidth',1)
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        axis([0 max(time_SQ) 0 3500]);grid on
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[N]'}, 'FontSize', 15)
        legend('without wobbling mass', 'with wobbling mass', 'FontSize', 13,'Location','northwest')
        
% y-components of Joint forces in hip & knee joint with and without wobbling mass - CMJ
    figure(11)
    sgtitle('Joint Forces CMJ with & without wobbling mass', 'FontSize', 20)
 
    subplot(1,2,1); plot(time_CM,FC_y_CM +(m*g),'r'); title('Hip - CMJ', 'FontSize', 15)
        hold on;  
        plot(time_CM,FC_y_wob_CM +(m*g),'k')
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        axis([0 max(time_CM) 250 3500]);grid on
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[N]'}, 'FontSize', 15)
        legend('without wobbling mass', 'with wobbling mass', 'FontSize', 13,'Location','northwest')
        
    subplot(1,2,2); plot(time_CM,FB_y_CM +(m*g),'r'); title('Knee - CMJ', 'FontSize', 15)
        hold on;  
        plot(time_CM,FB_y_wob_CM +(m*g),'k')
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        axis([0 max(time_CM) 250 3500]);grid on
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[N]'}, 'FontSize', 15)
        legend('without wobbling mass', 'with wobbling mass', 'FontSize', 13,'Location','northwest')
        
%% Plot: x-components & y-components of Joint Forces in SJ with and without wobbling mass

% x-components & y-components of Joint forces in hip & knee joint with and without wobbling mass - SJ
    figure(12)
    sgtitle('Joint Forces during SJ with and without wobbling mass', 'FontSize', 20)
 
    subplot(2,2,1); plot(time_SQ,FC_x_SQ,'r'); title('Hip x - SJ', 'FontSize', 15)
        hold on;  
        plot(time_SQ,FC_x_wob_SQ,'k')
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        axis([0 max(time_SQ) -1250 1250]);grid on
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[N]'}, 'FontSize', 15)
        legend('without wobbling mass', 'with wobbling mass', 'FontSize', 13)
        
    subplot(2,2,2); plot(time_SQ,FC_y_SQ,'r'); title('Hip y - SJ', 'FontSize', 15)
        hold on;  
        plot(time_SQ,FC_y_wob_SQ,'k')
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        axis([0 max(time_SQ) -250 2750]);grid on
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[N]'}, 'FontSize', 15)
        legend('without wobbling mass', 'with wobbling mass', 'FontSize', 13)
        
    subplot(2,2,3); plot(time_SQ,FB_x_SQ,'r'); title('Knee x - SJ', 'FontSize', 15)
        hold on;  
        plot(time_SQ,FB_x_wob_SQ,'k')
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        axis([0 max(time_SQ) -1250 1250]);grid on
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[N]'}, 'FontSize', 15)
        legend('without wobbling mass', 'with wobbling mass', 'FontSize', 13)
        
    subplot(2,2,4); plot(time_SQ,FB_y_SQ,'r'); title('Knee y - SJ', 'FontSize', 15)
        hold on;  
        plot(time_SQ,FB_y_wob_SQ,'k')
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        axis([0 max(time_SQ) -250 2750]);grid on
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[N]'}, 'FontSize', 15)
        legend('without wobbling mass', 'with wobbling mass', 'FontSize', 13)

%% Plot: x-components & y-components of Joint Forces in CMJ with and without wobbling mass

% x-components & y-components of Joint forces in hip & knee joint with and without wobbling mass - CMJ
    figure(13)
    sgtitle('Joint Forces during CMJ with and without wobbling mass', 'FontSize', 20)
 
    subplot(2,2,1); plot(time_CM,FC_x_CM,'r'); title('Hip x - CMJ', 'FontSize', 15)
        hold on;  
        plot(time_CM,FC_x_wob_CM,'k')
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        axis([0 max(time_CM) -1200 1500]);
        grid on
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[N]'}, 'FontSize', 15)
        legend('without wobbling mass', 'with wobbling mass', 'FontSize', 13,'Location','northwest')
        
    subplot(2,2,2); plot(time_CM,FC_y_CM,'r'); title('Hip y - CMJ', 'FontSize', 15)
        hold on;  
        plot(time_CM,FC_y_wob_CM,'k')
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        axis([0 max(time_CM) -250 2700]);
        grid on
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[N]'}, 'FontSize', 15)
        legend('without wobbling mass', 'with wobbling mass', 'FontSize', 13,'Location','northwest')
        
    subplot(2,2,3); plot(time_CM,FB_x_CM,'r'); title('Knee x - CMJ', 'FontSize', 15)
        hold on;  
        plot(time_CM,FB_x_wob_CM,'k')
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        axis([0 max(time_CM) -1200 1500]);
        grid on
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[N]'}, 'FontSize', 15)
        legend('without wobbling mass', 'with wobbling mass', 'FontSize', 13,'Location','northwest')
        
    subplot(2,2,4); plot(time_CM,FB_y_CM,'r'); title('Knee y - CMJ', 'FontSize', 15)
        hold on;  
        plot(time_CM,FB_y_wob_CM,'k')
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        axis([0 max(time_CM) -250 2700]);
        grid on
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[N]'}, 'FontSize', 15)
        legend('without wobbling mass', 'with wobbling mass', 'FontSize', 13,'Location','northwest')

%% Plot: moments in joints of SJ and CMJ

    figure(14)
    sgtitle('Moments in Joints', 'FontSize', 20)
    
    subplot(3,2,1); plot(time_SQ,MC_SQ, 'LineWidth', 3); title('Hip - SJ', 'FontSize', 15)
        grid on
        hold on
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        axis([0.01 max(time_SQ) -1000 1000]);
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[Nm]'}, 'FontSize', 15)
        
    subplot(3,2,2); plot(time_CM,MC_CM, 'LineWidth', 3); title('Hip - CMJ', 'FontSize', 15)
        grid on
        hold on
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        axis([0 max(time_CM) -1000 1000]);
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[Nm]'}, 'FontSize', 15)
        
    subplot(3,2,3); plot(time_SQ,MB_SQ, 'LineWidth', 3); title('Knee - SJ', 'FontSize', 15)
        grid on
        hold on
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        axis([0.01 max(time_SQ) -1000 1000]);
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[Nm]'}, 'FontSize', 15)
        
    subplot(3,2,4); plot(time_CM,MB_CM, 'LineWidth', 3); title('Knee - CMJ', 'FontSize', 15)
        grid on
        hold on
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        axis([0 max(time_CM) -1000 1000]);
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[Nm]'}, 'FontSize', 15)
        
    subplot(3,2,5); plot(time_SQ,MA_SQ, 'LineWidth', 3); title('Ankle - SJ', 'FontSize', 15)
        grid on
        hold on
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        axis([0.01 max(time_SQ) -1000 1000]);
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[Nm]'}, 'FontSize', 15)
    
    subplot(3,2,6); plot(time_CM,MA_CM, 'LineWidth', 3); title('Ankle - CMJ', 'FontSize', 15)
        grid on
        hold on
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        axis([0 max(time_CM) -1000 1000]);
        xlabel({'[s]'}, 'FontSize', 15)
        ylabel({'[Nm]'}, 'FontSize', 15)
        