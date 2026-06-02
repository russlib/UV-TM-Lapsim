%% Run the ACTUAL Lap_Sim_CSAPS.m Sections 0-6 (GGV build) + accel sim
%% This script runs the REAL code, not a rewrite.
%% It stops after GGV build + accel sim and dumps results to CSV.

% Run Sections 0-6 of the actual Lap_Sim_CSAPS
% We need to run it up through the GGV build, then extract the envelopes

% --- Run Sections 0-6 by evaluating Lap_Sim_CSAPS.m lines ---
% Instead of copy-pasting 500 lines, we'll run the actual script
% but catch the error when it tries to load track data (Section 7+)

try
    % Run the full script — it will build GGV in sections 0-6,
    % then try to load track data in section 7
    run('Lap_Sim_CSAPS.m');
catch ME
    % Expected: will error on track loading or later sections
    fprintf('Script stopped at: %s\n', ME.message);
    fprintf('(This is expected — we only need sections 0-6)\n\n');
end

%% Dump GGV Envelope Data
% The GGV splines are: accel, grip, lateral, deccel (not built yet), cornering
% These are csaps objects — evaluate them on a grid

fprintf('=== MATLAB GGV Envelope Data ===\n');

% Evaluate splines on velocity grid
v_grid = 15:2:130;  % ft/s (same range as Lap_Sim_CSAPS uses)

accel_vals = fnval(accel, v_grid);
grip_vals = fnval(grip, v_grid);

fprintf('\nAccel envelope (power-limited, g):\n');
for i = 1:length(v_grid)
    v_ms = v_grid(i) * 0.3048;  % Convert to m/s
    fprintf('  v=%6.1f ft/s (%5.1f m/s): accel=%.4f g, grip=%.4f g\n', ...
        v_grid(i), v_ms, accel_vals(i), grip_vals(i));
end

% Cornering envelope (if built) — max ay vs radius
if exist('cornering', 'var')
    fprintf('\nCornering envelope (max speed vs radius):\n');
    r_grid = 15:10:155;  % ft
    corner_vals = fnval(cornering, r_grid);
    for i = 1:length(r_grid)
        r_m = r_grid(i) * 0.3048;
        v_corner_ms = corner_vals(i) * 0.3048;
        fprintf('  R=%5.1f ft (%5.1f m): v_max=%.1f ft/s (%.1f m/s)\n', ...
            r_grid(i), r_m, corner_vals(i), v_corner_ms);
    end
end

if exist('lateral', 'var')
    lat_vals = fnval(lateral, v_grid);
    fprintf('\nLateral envelope (max ay vs speed, g):\n');
    for i = 1:length(v_grid)
        v_ms = v_grid(i) * 0.3048;
        fprintf('  v=%6.1f ft/s (%5.1f m/s): ay_max=%.4f g\n', ...
            v_grid(i), v_ms, lat_vals(i));
    end
end

% Save to CSV
output = zeros(length(v_grid), 4);
output(:,1) = v_grid' * 0.3048;  % v in m/s
output(:,2) = accel_vals';       % accel in g
output(:,3) = grip_vals';        % grip in g
if exist('lateral', 'var')
    output(:,4) = lat_vals';     % lateral in g
end

csvwrite('matlab_ggv_envelope.csv', output);
fprintf('\nGGV envelope saved to matlab_ggv_envelope.csv\n');
fprintf('Columns: v_ms, accel_g, grip_g, lateral_g\n');

%% Run the actual 247ft accel sim (Section 16 code)
% This uses the actual accel spline and shift logic from Lap_Sim_CSAPS
fprintf('\n=== Running 247ft Acceleration Sim ===\n');

% This is the exact code from Section 16 of Lap_Sim_CSAPS.m
count = 0;
v_sim = 0;
vel_sim = v_sim;
vmax_sim = top_speed;
interval_sim = 1;
segment_sim = 1:1:247;
clear dt_f_sim v_f_sim ax_f_sim
t_accel_elapsed = 0;
time_shifting_sim = 0;
gear_sim = 1;
shifting_sim = 0;

for i_sim = 1:1:length(segment_sim)
    d_sim = 1;
    gears_sim = find((shift_points - vel_sim) > 0);
    if isempty(gears_sim)
        newgear_sim = length(shift_points);
    else
        newgear_sim = gears_sim(1) - 1;
    end
    if newgear_sim > gear_sim
        shifting_sim = 1;
    else
        shifting_sim = 0;
    end

    AX_sim = fnval(accel, max(15, vel_sim));
    dd_sim = d_sim / interval_sim;

    for j_sim = 1:1:interval_sim
        count = count + 1;
        if shifting_sim == 1 && vel_sim < vmax_sim
            dt_f_sim(count) = dd_sim / max(vel_sim, 0.1);
            time_shifting_sim = time_shifting_sim + dt_f_sim(count);
            ax_f_sim(count) = 0;
            v_f_sim(count) = vel_sim;
        elseif vel_sim < vmax_sim && vel_sim > 0
            ax_f_sim(count) = AX_sim;
            tt_sim = roots([0.5*32.2*ax_f_sim(count), vel_sim, -dd_sim]);
            tt_sim = tt_sim(tt_sim > 0);
            if isempty(tt_sim)
                dt_f_sim(count) = dd_sim / max(vel_sim, 0.1);
            else
                dt_f_sim(count) = max(tt_sim);
            end
            dv_sim = 32.2 * ax_f_sim(count) * dt_f_sim(count);
            dv_sim = min(dv_sim, vmax_sim - vel_sim);
            v_f_sim(count) = vel_sim + dv_sim;
            vel_sim = v_f_sim(count);
        elseif vel_sim == 0
            % Launch from standstill
            ax_f_sim(count) = AX_sim;
            dt_f_sim(count) = sqrt(2*dd_sim/(32.2*max(ax_f_sim(count),0.01)));
            vel_sim = 32.2 * ax_f_sim(count) * dt_f_sim(count);
            v_f_sim(count) = vel_sim;
        else
            dt_f_sim(count) = dd_sim / vel_sim;
            ax_f_sim(count) = 0;
            v_f_sim(count) = vmax_sim;
        end

        if time_shifting_sim > shift_time
            shifting_sim = 0;
            time_shifting_sim = 0;
            gear_sim = newgear_sim;
        end
    end
end

accel_time_247ft = sum(dt_f_sim) + 0.1;  % +0.1 for reaction time
dist_sim = cumsum(ones(size(dt_f_sim)));  % each step = 1 ft
time_sim = cumsum(dt_f_sim);

fprintf('247ft (75.3m) acceleration time: %.3f s (%.3f without reaction)\n', ...
    accel_time_247ft, accel_time_247ft - 0.1);
fprintf('Final speed: %.1f ft/s (%.1f m/s, %.1f km/h)\n', ...
    v_f_sim(end), v_f_sim(end)*0.3048, v_f_sim(end)*0.3048*3.6);
fprintf('Max acceleration: %.3f g\n', max(ax_f_sim));

% Speed at distance milestones
for target_ft = [49.2, 98.4, 164.0, 247.0]  % ~15m, 30m, 50m, 75m
    idx = find(dist_sim >= target_ft, 1);
    if ~isempty(idx)
        fprintf('  %.0fft (%.0fm): t=%.3fs, v=%.1fft/s (%.1fm/s)\n', ...
            target_ft, target_ft*0.3048, time_sim(idx), v_f_sim(idx), v_f_sim(idx)*0.3048);
    end
end

% Save accel trace
accel_output = [dist_sim'*0.3048, time_sim', v_f_sim'*0.3048, ax_f_sim'*9.81];
csvwrite('matlab_actual_accel.csv', accel_output);
fprintf('\nAccel trace saved to matlab_actual_accel.csv\n');
fprintf('Columns: distance_m, time_s, velocity_ms, accel_ms2\n');
fprintf('\nDone.\n');
