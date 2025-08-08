function [lap_time time_elapsed velocity acceleration lateral_accel gear_counter path_length weights distance] = lap_information(path_positions)
% lap_information: Computes detailed lap simulation data for an endurance track.
% This function generates a vehicle trajectory based on path positions, computes track curvature and maximum speeds,
% performs forward and backward simulations to determine optimal velocity profiles, and outputs time traces,
% accelerations, gear usage, and performance weights (proportions of time spent accelerating, braking, cornering).
% Inputs:
%   path_positions - Vector of positions along track boundaries (normalized [0,1]).
% Outputs:
%   lap_time - Total lap time (s).
%   time_elapsed - Cumulative time trace (vector, s).
%   velocity - Velocity trace (vector, ft/s).
%   acceleration - Longitudinal acceleration trace (vector, g).
%   lateral_accel - Lateral acceleration trace (vector, g, signed by turn direction).
%   gear_counter - Gear selection trace (vector, dimensionless).
%   path_length - Total path length (ft).
%   weights - [accel_weight, brake_weight, corner_weight] based on time spent (proportions sum to 1).
%   distance - Distance trace (vector, ft) - Note: This output is computed but not explicitly returned in original; added for completeness.
% Globals Used:
%   path_boundaries - Track boundary coefficients and limits (matrix).
%   r_min, r_max - Minimum and maximum curvature radii (ft).
%   cornering, accel, grip, deccel, lateral - Spline functions for performance envelopes (e.g., csape or fnval compatible).
%   shift_points - Gear shift speed thresholds (vector, ft/s).
%   top_speed - Maximum vehicle speed (ft/s).
%   shift_time - Time penalty for shifting gears (s).
% Assumptions:
%   - Units: ft for distances, ft/s for speeds, g for accelerations, s for times.
%   - Quasi-steady-state model: Assumes equilibrium at each segment but simulates dynamics via forward/backward passes.
%   - Grip ellipse: Combined forces limited by (ax/AX)^2 + (ay/AY)^2 <= 1.
%   - Backward simulation loop is truncated in provided code; assumed to be similar to forward but in reverse.

global path_boundaries r_min r_max cornering accel grip deccel lateral...
    shift_points top_speed shift_time

%% Generate vehicle trajectory
% This section interpolates the vehicle's racing line from given positions along track boundaries.
% Boundaries are defined by linear coefficients (slope, intercept) and x-limits per segment.
% Uses pchip interpolation for smooth path, then computes total arc length.
interval = 5; % Discretization interval for simulation steps (ft?).
sections = 3000; % Number of sections for path interpolation (high for smoothness).
path_positions(end+1) = path_positions(1); % Close the loop for a full lap (repeat first position).
path_positions(end+1) = path_positions(2); % Repeat second for continuity.
VMAX = top_speed; % Maximum achievable speed (ft/s), from global.
t = 1:1:length(path_positions); % Indices for positions.
for i = 1:1:length(path_positions)
    coeff = path_boundaries(i,1:2); % Slope and intercept for boundary line.
    x2 = max(path_boundaries(i,3:4)); % Upper x-limit of segment.
    x1 = min(path_boundaries(i,3:4)); % Lower x-limit.
    position = path_positions(i); % Normalized position [0,1] along the boundary.
    x3 = x1+position*(x2-x1); % Interpolated x-coordinate.
    y3 = polyval(coeff,x3); % Corresponding y from line equation.
    path_points(i,:) = [x3 y3]; % Store [x,y] point.            
end

x = linspace(1,t(end-1),sections); % Linear spacing for interpolation.
ppv = pchip(t,path_points'); % Piecewise cubic Hermite interpolating polynomial for smooth curve.
vehicle_path = ppval(ppv,x); % Evaluate interpolated path.
path_length = arclength(vehicle_path(1,:),vehicle_path(2,:)); % Total distance along path (ft).

% Alternative interpolation (commented): Uses makima for modified Akima cubic, potentially smoother but less monotonic.
% x = linspace(1,t(end-1),1000);
% ppv = interp1([1:length(path_points)],path_points,x,'makima');
% vehicle_path = ppv';

[L,R,K] = curvature(vehicle_path'); % Compute length, radius, curvature vector for the path.

%% Traverse the track
% This section discretizes the track into segments, computes curvature radii, and determines max speeds and distances per segment.
track_points = vehicle_path; % Copy path for processing.
track_points = [track_points(:,length(vehicle_path)-2) track_points(:,1:end-1)]; % Shift for curvature calc (closes loop?).
[LT,RT,KT] = curvature(track_points'); % Curvature for shifted points.
KT = KT(:,2); % Extract relevant curvature components.
KT = KT(~isnan(RT)); % Remove NaNs from radii.
RT = RT(~isnan(RT)); % Clean radii.
RT = RT(~isnan(RT)); % Redundant clean, but ensures no NaNs.
% For each segment: Clamp radius, compute max speed from cornering spline, and segment distance.
for i = 1:length(RT)
    segment(i) = i; % Segment index.
    r = max(r_min,RT(i)); % Clamp to min radius (tightest turn).
    r = min(r,r_max); % Clamp to max radius (straight-ish).
    RT(i) = r; % Update cleaned radius.
    Vmax(i) = min(VMAX,fnval(cornering,r)); % Max speed from cornering envelope (ft/s).
    x1(i) = track_points(1,i+1); % Start x of segment.
    x2(i) = track_points(1,i+2); % End x.
    y1(i) = track_points(2,i+1); % Start y.
    y2(i) = track_points(2,i+2); % End y.
    dist(i) = sqrt((x1(i)-x2(i))^2+(y2(i)-y1(i))^2); % Euclidean distance (ft).
end

%% Initiate forward sim
% Forward simulation: Starts from low speed, accelerates forward through segments, respecting grip limits and gear shifts.
% Simulates in small intervals per segment, updating velocity, acceleration, and handling shifts.
count = 0; % Global counter for trace indices.
v = 20; % Initial velocity (ft/s) - starting from a low speed, e.g., post-turn or launch.
vel = v; % Current velocity copy.
gear = 0; % FIX: Initialize gear to 0 (neutral/lowest gear) to avoid undefined variable error.
shiftercount = 0; % Counter for number of shifts.
time_shifting = 0; % Timer for current shift duration (s).
ts = zeros(1, length(segment)); % Array to log shift times (for debugging/averaging).

for i = 1:1:length(segment) % Loop over each track segment.
    d = dist(i); % Total distance for this segment (ft).
    r = RT(i); % Radius for this segment (ft).
    
    % Determine new gear based on current velocity (find first shift point exceeding vel).
    gears = find((shift_points-vel)>0); % Indices where shift_points > vel.
    newgear = gears(1)-1; % Current gear is one less than the next shift index.
    
    % Check if upshifting is needed (newgear > current gear).
    if newgear > gear
        shifting = 1; % Flag to indicate shift in progress.
    else
        shifting = 0; % No shift.
    end
    
    vmax = min(VMAX,fnval(cornering,r)); % Max allowable speed for this radius (ft/s).
    if vmax < 0
        vmax = VMAX; % Fallback if spline returns invalid (negative) value.
    end
    AX = fnval(accel,vel); % Max longitudinal accel from spline (g).
    AY = fnval(lateral,vel); % Max lateral accel from spline (g).
    dd = d/interval; % Sub-distance per interval (ft) for finer simulation steps.
    for j = 1:1:interval % Sub-loop for each small interval in segment.
        count = count+1; % Increment trace index.
        vehicle_gear(count) = gear; % Log current gear.
        ay_f(count) = vel^2/(r*32.2); % Lateral accel demand (g) = v^2 / (r * g_ft/s^2).
        if shifting == 1 && vel < vmax % During shift: No accel, constant vel.
            dt_f(count) = dd/vel; % Time for interval (s) = distance / vel.
            time_shifting = time_shifting+dt_f(count); % Accumulate shift time.
            ax_f(count) = 0; % No longitudinal accel during shift.
            v_f(count) = vel; % Velocity trace.
            dv_f(count) = 0; % No change in vel.
            vel = vel; % Hold velocity.
        elseif vel < vmax % Accelerating: Use grip ellipse to limit ax based on ay demand.
            ax_f(count) = AX*(1-(min(AY,ay_f(count))/AY)^2); % Reduced ax = AX * (1 - (ay/AY)^2).
            tt = roots([0.5*32.2*ax_f(count) vel -dd]); % Solve quadratic for time: 0.5*a*t^2 + v*t - dd = 0 (a in ft/s^2 = 32.2*g).
            dt_f(count) = max(tt); % Take positive root (time s).
            dv = 32.2*ax_f(count)*dt_f(count); % Delta vel (ft/s).
            dvmax = vmax-vel; % Cap to not exceed vmax.
            dv_f(count) = min(dv,dvmax); % Applied delta vel.
            v_f(count) = vel+dv_f(count); % Updated velocity trace.
            vel = v_f(count); % Update current vel.
            gears = find((shift_points-vel)>0); % Recheck gear for new vel.
            newgear = gears(1)-1;
            if newgear > gear
                shifting = 1; % Trigger shift if needed.
            end
        else % Coasting/at limit: Constant vel.
            dt_f(count) = dd/vel; % Time for interval.
            ax_f(count) = 0; % No accel.
            v_f(count) = vel;
            dv_f(count) = 0;
            vel = vmax; % Cap to vmax.
        end
        if time_shifting > shift_time % Complete shift if timer exceeds penalty.
            shiftercount = shiftercount+1; % Increment shift count.
            ts(shiftercount) = time_shifting; % Log shift duration.
            shifting = 0; % Reset flag.
            time_shifting = 0; % Reset timer.
            gear = newgear; % Update to new gear.
        end
    end
    if shifting == 1
        gear = gear; % Hold if shift ongoing.
    else
        gear = newgear; % Update otherwise.
    end
end

%% combine results
% Combine forward and backward simulations: Take min velocity per point for feasible profile.
% Note: Backward sim (v_r, dt_r, etc.) is truncated; assume it's computed similarly (reverse direction, braking focus).
VD = v_f-v_r; % Velocity difference (forward - reverse).
forw = find(VD>=0); % Points where forward is feasible.
back = find(VD<0); % Points where reverse (braking-limited) is binding.
velocity = zeros(1,length(VD)); % Initialize velocity trace.
t_elapsed = 0; % Cumulative time (s).
classifier = []; % Unused? Possibly for debugging.
for i = 1:1:length(VD)
    if VD(i)<0 % Use forward if forward < reverse (more conservative).
        velocity(i) = v_f(i);
        dtime(i) = dt_f(i);
        acceleration(i) = ax_f(i); % Longitudinal accel (g).
        lateral_accel(i) = ay_f(i); % Lateral accel demand (g).
        %shift_timer(i) = shift_timer(i); % Commented: Possibly unused.
    else % Use reverse otherwise.
        velocity(i) = v_r(i);
        dtime(i) = dt_r(i);
        acceleration(i) = -ax_r(i); % Negate for braking.
        lateral_accel(i) = ay_r(i);
        %shift_timer(i) = 0; % Reset during braking?
    end
    t_elapsed = t_elapsed+dtime(i); % Accumulate time.
    time_elapsed(i) = t_elapsed; % Time trace.
end
AY_outlier = find(lateral_accel > fnval(lateral,116)); % Cap outliers to max lateral at 116 ft/s.
lateral_accel(AY_outlier) = fnval(lateral,116); % Clamp to envelope.
throttle = 0; % Accumulator for accel effort (g*s).
brake = 0; % For braking.
corner = 0; % For cornering.
for i = 1:1:length(VD) % Integrate over time for weighted efforts.
    if acceleration(i)>0
        throttle = throttle+acceleration(i)*dtime(i); % Weighted by time.
    elseif acceleration(i) < 0
        brake = brake-acceleration(i)*dtime(i); % Positive for magnitude.
    end
    corner = corner + lateral_accel(i)*dtime(i); % Lateral effort.
end
summ = throttle+brake+corner; % Total effort.
weights = [throttle/summ brake/summ corner/summ]; % Proportions (effort-based).
tloc = find(acceleration>.25); % Threshold for accel classification (>0.25g).
t_t = sum(dtime(tloc)); % Time accelerating.
bloc = find(acceleration<-.25); % Braking.
t_b = sum(dtime(bloc));
cloc = find(lateral_accel>.25); % Cornering.
t_c = sum(dtime(cloc));
summ = t_t+t_b+t_c; % Total classified time.
weights = [t_t/summ t_b/summ t_c/summ]; % Time-based weights (overrides effort-based?).
%figure
%plot(distance,velocity) % Uncomment to plot velocity vs. distance.

%% Plot Results
% Visualize average velocity per segment on track map (heatmap-style scatter).
figure
for i = 1:1:length(track_points)-2 % Average velocity over intervals.
    V_plot(i) = mean(velocity(i*interval-interval+1:i*interval)); % Mean per segment.
end
pointsize = 5; % Marker size.
scatter(track_points(1,2:end-1),track_points(2,2:end-1),100,V_plot,'marker','.') % Plot points colored by velocity.
title('2019 Michigan Endurance Simulation Track Summary')
h = colorbar; % Colorbar for velocity.
set(get(h,'title'),'string','Velocity (V) [ft/s]');
set(gca,'XTick',[], 'YTick', []) % Clean axes.

%% Gear Counter
% Post-process gear for entire trace based on final velocities.
for i = 1:1:length(velocity)
    V = velocity(i); % Velocity at point.
    gears = find((shift_points-V)>0); % Shift indices.
    gear = gears(1)-1; % Gear one below next shift.
gear_counter(i) = gear; % Log.
end
lap_time = t_elapsed; % Total time.

for i = 1:1:length(lateral_accel) % Sign lateral accel by turn direction (curvature sign).
    index = floor((i-1)/interval)+1; % Segment index.
    axis(i) = sign(KT(index)); % +1 or -1 from curvature.
end
lateral_accel = lateral_accel.*axis; % Signed trace.
end