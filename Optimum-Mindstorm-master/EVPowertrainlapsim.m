function output = EVPowertrainlapsim(initialV)
% For the love of the FSM, keep everything in metric
% initialV is in m/s **********************************
global powertrainpackage tire_radius
% Input RPM and torque [N-m] curve from dyno. Got rid of xlsread for run time reasons
engineSpeed = powertrainpackage{1};  
engineTq = powertrainpackage{2}; 

tyreRadius = tire_radius/3.28;
%tyreRadius = 0.45974/2;     % metres
% Gear Ratio
finalDrive = powertrainpackage{3};

i=1;
while engineSpeed < ((initialV*60/tyreRadius)/(2*3.1416))/finalDrive
    i=i+1;
end

% Variable declaration
drivetrainLosses = powertrainpackage{4}; 
torque = engineTq(i)*finalDrive*drivetrainLosses;   % torque output at the wheels [N-m]
Fx = torque/tyreRadius;      % force on contact patch from drivetrain [N]

Power = Fx*initialV;
PowerLimit = powertrainpackage{5};

if Power>PowerLimit
    Fx = PowerLimit/initialV;
end

if ((initialV*60/tyreRadius)/(2*3.1416)) > (engineSpeed(end)/finalDrive)
    Fx = 0;
end

output = [Fx]; % [N]