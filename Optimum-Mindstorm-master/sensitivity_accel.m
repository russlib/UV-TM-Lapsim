%% MATLAB Longitudinal Sensitivity: 75m accel time vs CG height, wheelbase
%% Uses the ACTUAL Lap_Sim_CSAPS GGV build + accel sim code
clear; clc;

global r_max accel grip deccel lateral cornering gear shift_points...
    top_speed r_min path_boundaries tire_radius shift_time...
    powertrainpackage track_width path_boundaries_ax

%% Fixed params
load('Hoosier_R25B_18.0x7.5-10_FY_12psi.mat');
load('Hoosier_R25B_18.0x7.5-10_FX_12psi.mat');
tire_radius = 9.05/12;
sf_x = .6; sf_y = .47;

engineSpeed = [6200:100:14100];
engineTq = [41.57 42.98 44.43 45.65 46.44 47.09 47.52 48.58 49.57 50.41 51.43 51.48 51 49.311 48.94 48.66 49.62 49.60 47.89 47.91 48.09 48.57 49.07 49.31 49.58 49.56 49.84 50.10 50.00 50.00 50.75 51.25 52.01 52.44 52.59 52.73 53.34 53.72 52.11 52.25 51.66 50.5 50.34 50.50 50.50 50.55 50.63 50.17 50.80 49.73 49.35 49.11 48.65 48.28 48.28 47.99 47.68 47.43 47.07 46.67 45.49 45.37 44.67 43.8 43.0 42.3 42.00 41.96 41.70 40.43 39.83 38.60 38.46 37.56 36.34 35.35 33.75 33.54 32.63 31.63];
primaryReduction = 76/36;
gear_ratios = [33/12, 32/16, 30/18, 26/18, 30/23, 29/24];
finalDrive = 40/12;
shiftpoint = 14000;
drivetrainLosses = .85;
shift_time = .25;
T_lock = 90/100;
tyreRadius = tire_radius/3.28;
gearTot = gear_ratios(end)*finalDrive*primaryReduction;
VMAX_base = floor(3.28*shiftpoint/(gearTot/tyreRadius*60/(2*pi)));
powertrainpackage = {engineSpeed engineTq primaryReduction gear_ratios finalDrive shiftpoint drivetrainLosses};

W = 660; WDF = 50/100;
Cl = .0418; Cd = .0184; CoP = 48/100;
rg_f=0; rg_r=0; pg=0; WRF=180; WRR=180;
IA_staticf=0; IA_staticr=0; IA_compensationf=10/100; IA_compensationr=20/100;
casterf=0; KPIf=0; casterr=deg2rad(4.1568); KPIr=0;
IA_roll_inducedf = asin(2/(46/12)/12);
IA_roll_inducedr = asin(2/(44/12)/12);
IA_gainf = IA_roll_inducedf*IA_compensationf;
IA_gainr = IA_roll_inducedr*IA_compensationr;

WF = W*WDF; WR = W*(1-WDF);

function t75 = run_accel_sim(cg_ft, l_ft, twf_ft, twr_ft, LLTD_pct, ...
    W, WDF, WF, WR, Cl, Cd, CoP, WRF, WRR, IA_staticf, IA_staticr, ...
    IA_gainf, IA_gainr, pg, sf_x, full_send_x, tire_radius, ...
    powertrainpackage, shiftpoint, shift_time, VMAX)

    global accel grip gear shift_points top_speed

    velocity = 15:5:130;
    g_gear = 1; spcount = 1; shift_points(1) = 0;
    clear A_xr A_Xr fx_pwr gear

    for i = 1:length(velocity)
        gp = g_gear; V = velocity(i);
        DF = Cl*V^2;
        dxf = DF*CoP/2/WRF; dxr = DF*(1-CoP)/2/WRR;
        IA_0f = IA_staticf - dxf*IA_gainf;
        IA_0r = IA_staticr - dxr*IA_gainr;
        wf = (WF+DF*CoP)/2; wr = (WR+DF*(1-CoP))/2;
        Ax = 0; WS = W/2; pitch = -Ax*pg*pi/180;
        wf = wf-Ax*cg_ft*WS/l_ft; wr = wr+Ax*cg_ft*WS/l_ft;
        IA_f = -l_ft*12*sin(pitch)/2*IA_gainf + IA_0f;
        IA_r = l_ft*12*sin(pitch)/2*IA_gainr + IA_0r;
        sl = [0:.01:.11]; fxf_v=[]; fxr_v=[];
        for k = 1:length(sl)
            fxf_v(k) = fnval([sl(k);-wf;rad2deg(-IA_f)],full_send_x)*sf_x;
            fxr_v(k) = fnval([sl(k);-wr;rad2deg(-IA_r)],full_send_x)*sf_x;
        end
        fxf_v(abs(fxf_v)>1000)=[]; fxr_v(abs(fxr_v)>1000)=[];
        FXR = max(fxr_v); FX = abs(2*FXR); AX = FX/W; AX_diff = AX-Ax;
        while AX_diff>0
            Ax=Ax+.01; WS=W/2; pitch=-Ax*pg*pi/180;
            wf=(WF+DF*CoP)/2; wr=(WR+DF*(1-CoP))/2;
            wf=wf-Ax*cg_ft*WS/l_ft/24; wr=wr+Ax*cg_ft*WS/l_ft/24;
            IA_f=-l_ft*12*sin(pitch)/2*IA_gainf+IA_0f;
            IA_r=l_ft*12*sin(pitch)/2*IA_gainr+IA_0r;
            sl=[0:.01:.11]; fxf_v=[]; fxr_v=[];
            for k=1:length(sl)
                fxf_v(k)=fnval([sl(k);-wf;rad2deg(-IA_f)],full_send_x)*sf_x;
                fxr_v(k)=fnval([sl(k);-wr;rad2deg(-IA_r)],full_send_x)*sf_x;
            end
            fxf_v(abs(fxf_v)>1000)=[]; fxr_v(abs(fxr_v)>1000)=[];
            FXR=max(fxr_v); FX=abs(2*FXR); AX=FX/W; AX_diff=AX-Ax;
        end
        A_xr(i) = AX;
        output = Powertrainlapsim(max(7.5,V/3.28));
        FX_pwr=output(1)*.2248; FX_pwr=FX_pwr-Cd*V^2;
        AX_i=min(FX_pwr/W,A_xr(i));
        output = Powertrainlapsim(V/3.28);
        g_gear = output(2); gear(i) = g_gear;
        if g_gear>gp; spcount=spcount+1; shift_points(spcount)=V; end
        A_Xr(i) = max(AX_i, 0);
    end
    accel = csaps(velocity, A_Xr);
    grip = csaps(velocity, A_xr);
    top_speed = VMAX;

    % Run 247ft accel
    count=0; vel_s=0; vmax_s=VMAX; gear_s=1; shifting_s=0; time_shifting_s=0;
    for i_s=1:247
        gears_s=find((shift_points-vel_s)>0);
        if isempty(gears_s); newgear_s=length(shift_points);
        else; newgear_s=gears_s(1)-1; end
        if newgear_s>gear_s; shifting_s=1; end
        AX_s=fnval(accel,max(15,vel_s));
        count=count+1;
        if shifting_s==1 && vel_s<vmax_s
            dt_f_s(count)=1/max(vel_s,0.01);
            time_shifting_s=time_shifting_s+dt_f_s(count);
            v_f_s(count)=vel_s;
        elseif vel_s<vmax_s && vel_s>0
            tt_s=roots([0.5*32.2*AX_s,vel_s,-1]);
            tt_s=tt_s(tt_s>0);
            if isempty(tt_s); dt_f_s(count)=1/max(vel_s,0.01);
            else; dt_f_s(count)=max(tt_s); end
            dv_s=32.2*AX_s*dt_f_s(count);
            dv_s=min(dv_s,vmax_s-vel_s);
            v_f_s(count)=vel_s+dv_s; vel_s=v_f_s(count);
        elseif vel_s==0
            dt_f_s(count)=sqrt(2/(32.2*max(AX_s,0.01)));
            vel_s=32.2*AX_s*dt_f_s(count); v_f_s(count)=vel_s;
        else
            dt_f_s(count)=1/vel_s; v_f_s(count)=vmax_s;
        end
        if time_shifting_s>shift_time
            shifting_s=0; time_shifting_s=0; gear_s=newgear_s;
        end
    end
    t75 = sum(dt_f_s);
end

%% CG Height sweep
cg_values_ft = [150 200 250 300 350 400] / 304.8;  % mm to ft
l_ft = 60.5/12;

fprintf('=== CG Height Sweep (MATLAB CSAPS) ===\n');
fprintf('%10s | %12s | %12s\n', 'h_cg (mm)', 't_75m (s)', 'delta (s)');
fprintf('%s\n', repmat('-', 1, 40));

t_base_cg = 0;
for idx = 1:length(cg_values_ft)
    cg_v = cg_values_ft(idx);
    t = run_accel_sim(cg_v, l_ft, 46/12, 44/12, 51.5, ...
        W, WDF, WF, WR, Cl, Cd, CoP, WRF, WRR, IA_staticf, IA_staticr, ...
        IA_gainf, IA_gainr, pg, sf_x, full_send_x, tire_radius, ...
        powertrainpackage, shiftpoint, shift_time, VMAX_base);
    if idx == 4; t_base_cg = t; end  % 300mm baseline
    fprintf('%10.1f | %12.3f | %+12.3f\n', cg_v*304.8, t, t - t_base_cg);
end

%% Wheelbase sweep
wb_values_ft = [1300 1400 1500 1600 1700 1800] / 304.8;  % mm to ft

fprintf('\n=== Wheelbase Sweep (MATLAB CSAPS) ===\n');
fprintf('%10s | %12s | %12s\n', 'L (mm)', 't_75m (s)', 'delta (s)');
fprintf('%s\n', repmat('-', 1, 40));

t_base_wb = 0;
for idx = 1:length(wb_values_ft)
    l_v = wb_values_ft(idx);
    t = run_accel_sim(13.2/12, l_v, 46/12, 44/12, 51.5, ...
        W, WDF, WF, WR, Cl, Cd, CoP, WRF, WRR, IA_staticf, IA_staticr, ...
        IA_gainf, IA_gainr, pg, sf_x, full_send_x, tire_radius, ...
        powertrainpackage, shiftpoint, shift_time, VMAX_base);
    if idx == 3; t_base_wb = t; end  % 1500mm baseline
    fprintf('%10.1f | %12.3f | %+12.3f\n', l_v*304.8, t, t - t_base_wb);
end

fprintf('\nDone.\n');
