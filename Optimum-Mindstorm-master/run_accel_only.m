%% Run ACTUAL Lap_Sim_CSAPS Sections 0-2 + accel envelope + accel sim
%% Skips the cornering envelope (which can hang) — tests powertrain+tire accel only
clear; clc;

global r_max accel grip deccel lateral cornering gear shift_points...
    top_speed r_min path_boundaries tire_radius shift_time...
    powertrainpackage track_width path_boundaries_ax

%% Section 1: Tire Model (ACTUAL CODE from Lap_Sim_CSAPS.m)
disp('Loading Tire Model')
filename = 'Hoosier_R25B_18.0x7.5-10_FY_12psi.mat';
load(filename)
filename = 'Hoosier_R25B_18.0x7.5-10_FX_12psi.mat';
load(filename)
tire_radius = 9.05/12;
tyreRadius = tire_radius/3.28;
sf_x = .6;
sf_y = .47;

%% Section 2: Powertrain (ACTUAL CODE)
disp('Loading Engine Model')
engineSpeed = [6200:100:14100];
engineTq = [41.57 42.98 44.43 45.65 46.44 47.09 47.52 48.58 49.57 50.41 51.43 51.48 51 49.311 48.94 48.66 49.62 49.60 47.89 47.91 48.09 48.57 49.07 49.31 49.58 49.56 49.84 50.10 50.00 50.00 50.75 51.25 52.01 52.44 52.59 52.73 53.34 53.72 52.11 52.25 51.66 50.5 50.34 50.50 50.50 50.55 50.63 50.17 50.80 49.73 49.35 49.11 48.65 48.28 48.28 47.99 47.68 47.43 47.07 46.67 45.49 45.37 44.67 43.8 43.0 42.3 42.00 41.96 41.70 40.43 39.83 38.60 38.46 37.56 36.34 35.35 33.75 33.54 32.63 31.63];
primaryReduction = 76/36;
gear = [33/12, 32/16, 30/18, 26/18, 30/23, 29/24];
finalDrive = 40/12;
shiftpoint = 14000;
drivetrainLosses = .85;
shift_time = .25;
T_lock = 90;
gearTot = gear(end)*finalDrive*primaryReduction;
VMAX = floor(3.28*shiftpoint/(gearTot/tyreRadius*60/(2*pi)));
T_lock = T_lock/100;
powertrainpackage = {engineSpeed engineTq primaryReduction gear finalDrive shiftpoint drivetrainLosses};

%% Section 3: Vehicle Architecture (ACTUAL CODE)
disp('Loading Vehicle Characteristics')
LLTD = 51.5;
W = 660;
WDF = 50;
cg = 13.2/12;
l = 60.5/12;
twf = 46/12;
twr = 44/12;
LLTD = LLTD/100;
WDF = WDF/100;
m = W/32.2;
WF = W*WDF;
WR = W*(1-WDF);
a = l*(1-WDF);
b = l*WDF;
tw = twf;
track_width = twf;

%% Section 4: Suspension (ACTUAL — zeroed out as in original)
rg_f = 0; rg_r = 0; pg = 0;
WRF = 180; WRR = 180;
IA_staticf = 0; IA_staticr = 0;
IA_compensationf = 10; IA_compensationr = 20;
casterf = 0; KPIf = 0; casterr = 4.1568; KPIr = 0;
IA_staticf = deg2rad(IA_staticf);
IA_staticr = deg2rad(IA_staticr);
IA_compensationf = IA_compensationf/100;
IA_compensationr = IA_compensationr/100;
casterf = deg2rad(casterf);
KPIf = deg2rad(KPIf);
casterr = deg2rad(casterr);
KPIr = deg2rad(KPIr);
IA_roll_inducedf = asin(2/twf/12);
IA_roll_inducedr = asin(2/twr/12);
IA_gainf = IA_roll_inducedf*IA_compensationf;
IA_gainr = IA_roll_inducedr*IA_compensationr;

%% Section 5: Aero (ACTUAL CODE)
Cl = .0418;
Cd = .0184;
CoP = 48;
CoP = CoP/100;

%% Section 6: Build Acceleration Envelope ONLY (ACTUAL CODE from Lap_Sim_CSAPS lines 140-238)
disp('Building Acceleration Envelope (ACTUAL tire spline code)')
deltar = 0; deltaf = 0;
velocity = 15:5:130;
g_gear = 1;
spcount = 1;
shift_points(1) = 0;

for i = 1:1:length(velocity)
    gp = g_gear;
    V = velocity(i);
    DF = Cl*V^2;
    dxf = DF*CoP/2/WRF;
    dxr = DF*(1-CoP)/2/WRR;
    IA_0f = IA_staticf - dxf*IA_gainf;
    IA_0r = IA_staticr - dxr*IA_gainr;
    wf = (WF+DF*CoP)/2;
    wr = (WR+DF*(1-CoP))/2;
    Ax = 0;
    WS = W/2;
    pitch = -Ax*pg*pi/180;
    wf = wf-Ax*cg*WS/l;
    wr = wr+Ax*cg*WS/l;
    IA_f = -l*12*sin(pitch)/2*IA_gainf + IA_0f;
    IA_r = l*12*sin(pitch)/2*IA_gainr + IA_0r;
    sl = [0:.01:.11];
    for k = 1:length(sl)
        fxf(k) = fnval([sl(k);-wf;rad2deg(-IA_f)],full_send_x)*sf_x;
        fxr(k) = fnval([sl(k);-wr;rad2deg(-IA_r)],full_send_x)*sf_x;
    end
    fxf(find(abs(fxf) > 1000)) = [];
    fxr(find(abs(fxr) > 1000)) = [];
    FXF = max(fxf);
    FXR = max(fxr);
    FX = abs(2*FXR);
    AX = FX/W;
    AX_diff = AX-Ax;
    while AX_diff>0
        Ax = Ax+.01;
        WS = W/2;
        pitch = -Ax*pg*pi/180;
        wf = (WF+DF*CoP)/2;
        wr = (WR+DF*(1-CoP))/2;
        wf = wf-Ax*cg*WS/l/24;
        wr = wr+Ax*cg*WS/l/24;
        IA_f = -l*12*sin(pitch)/2*IA_gainf + IA_0f;
        IA_r = l*12*sin(pitch)/2*IA_gainr + IA_0r;
        sl = [0:.01:.11];
        for k = 1:length(sl)
            fxf(k) = fnval([sl(k);-wf;rad2deg(-IA_f)],full_send_x)*sf_x;
            fxr(k) = fnval([sl(k);-wr;rad2deg(-IA_r)],full_send_x)*sf_x;
        end
        fxf(find(abs(fxf) > 1000)) = [];
        fxr(find(abs(fxr) > 1000)) = [];
        FXF = max(fxf);
        FXR = max(fxr);
        FX = abs(2*FXR);
        AX = FX/W;
        AX_diff = AX-Ax;
    end
    A_xr(i) = AX;
    output = Powertrainlapsim(max(7.5,V/3.28));
    FX_pwr = output(1)*.2248;
    FX_pwr = FX_pwr-Cd*V^2;
    fx_pwr(i) = FX_pwr/W;
    AX_i = min(FX_pwr/W, A_xr(i));
    output = Powertrainlapsim(V/3.28);
    g_gear = output(2);
    gear(i) = g_gear;
    if g_gear>gp
        spcount = spcount+1;
        shift_points(spcount) = V;
    end
    A_Xr(i) = AX_i;
end
A_Xr(A_Xr < 0) = 0;
accel = csaps(velocity,A_Xr);
grip = csaps(velocity,A_xr);
top_speed = VMAX;

%% Dump acceleration envelope
fprintf('\n=== MATLAB Acceleration Envelope (ACTUAL tire model) ===\n');
v_eval = 15:5:130;
accel_g = fnval(accel, v_eval);
grip_g = fnval(grip, v_eval);
fprintf('%10s %12s %12s %12s\n', 'v_ftps', 'v_ms', 'accel_g', 'grip_g');
for i = 1:length(v_eval)
    fprintf('%10.1f %12.2f %12.4f %12.4f\n', v_eval(i), v_eval(i)*0.3048, accel_g(i), grip_g(i));
end

% Save envelope
env_out = [v_eval'*0.3048, accel_g', grip_g'];
csvwrite('matlab_accel_envelope.csv', env_out);
fprintf('\nSaved to matlab_accel_envelope.csv (v_ms, accel_g, grip_g)\n');

%% Run 247ft accel sim (ACTUAL code from Section 16)
fprintf('\n=== 247ft Acceleration Sim (ACTUAL code) ===\n');
count = 0;
v_s = 0; vel_s = v_s;
vmax_s = top_speed;
interval_s = 1;
segment_s = 1:1:247;
time_shifting_s = 0;
gear_s = 1;
shifting_s = 0;

for i_s = 1:length(segment_s)
    d_s = 1;
    gears_s = find((shift_points - vel_s) > 0);
    if isempty(gears_s)
        newgear_s = length(shift_points);
    else
        newgear_s = gears_s(1) - 1;
    end
    if newgear_s > gear_s
        shifting_s = 1;
    end
    AX_s = fnval(accel, max(15, vel_s));
    dd_s = d_s / interval_s;
    for j_s = 1:interval_s
        count = count + 1;
        if shifting_s == 1 && vel_s < vmax_s
            dt_f_s(count) = dd_s / max(vel_s, 0.01);
            time_shifting_s = time_shifting_s + dt_f_s(count);
            ax_f_s(count) = 0;
            v_f_s(count) = vel_s;
        elseif vel_s < vmax_s && vel_s > 0
            ax_f_s(count) = AX_s;
            tt_s = roots([0.5*32.2*ax_f_s(count), vel_s, -dd_s]);
            tt_s = tt_s(tt_s > 0);
            if isempty(tt_s); dt_f_s(count) = dd_s/max(vel_s,0.01);
            else; dt_f_s(count) = max(tt_s); end
            dv_s = 32.2 * ax_f_s(count) * dt_f_s(count);
            dv_s = min(dv_s, vmax_s - vel_s);
            v_f_s(count) = vel_s + dv_s;
            vel_s = v_f_s(count);
        elseif vel_s == 0
            ax_f_s(count) = AX_s;
            dt_f_s(count) = sqrt(2*dd_s/(32.2*max(ax_f_s(count),0.01)));
            vel_s = 32.2 * ax_f_s(count) * dt_f_s(count);
            v_f_s(count) = vel_s;
        else
            dt_f_s(count) = dd_s / vel_s;
            ax_f_s(count) = 0;
            v_f_s(count) = vmax_s;
        end
        if time_shifting_s > shift_time
            shifting_s = 0;
            time_shifting_s = 0;
            gear_s = newgear_s;
        end
    end
end

accel_time = sum(dt_f_s) + 0.1;
time_trace = cumsum(dt_f_s);
dist_trace = cumsum(ones(size(dt_f_s)));

fprintf('247ft (75.3m) time: %.3fs (%.3fs without reaction)\n', accel_time, accel_time-0.1);
fprintf('Final speed: %.1f ft/s = %.1f m/s = %.1f km/h\n', v_f_s(end), v_f_s(end)*0.3048, v_f_s(end)*0.3048*3.6);
fprintf('Max accel: %.4f g\n', max(ax_f_s));
fprintf('Shift points (ft/s): '); fprintf('%.1f ', shift_points); fprintf('\n');

for tgt = [49.2, 98.4, 164.0, 247.0]
    idx = find(dist_trace >= tgt, 1);
    if ~isempty(idx)
        fprintf('  %5.1fft (%4.1fm): t=%.3fs, v=%.1fft/s (%.1fm/s)\n', ...
            tgt, tgt*0.3048, time_trace(idx), v_f_s(idx), v_f_s(idx)*0.3048);
    end
end

% Save
out = [dist_trace'*0.3048, time_trace', v_f_s'*0.3048, ax_f_s'*9.81];
csvwrite('matlab_actual_accel.csv', out);
fprintf('\nSaved to matlab_actual_accel.csv\n');
fprintf('Done.\n');
