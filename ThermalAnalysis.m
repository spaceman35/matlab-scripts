%% Constants
c = 921.096/2; %J/kgK - Heat capacity for Aluminum/2 as a lumped therm mass assumption
Ap_max = 0.2*0.3; %m^2
Ap_min = 0.1*0.2; %m^2
Ar = 2*(0.1*0.2+0.2*0.3+0.3*0.1); %m^2 Radiative area
S = 1353; %W/m^2 Solar flux
S_ecl = 0; % Solar flux in eclipse
Re = 6378; %km - Earth radius
mu = 398600.4418; % km^3/s^2

%% Variables
m = 10.857; %kg CubeSat mass
albedo = 0.33; % average Earth albedo
Qin_ecl = 17.2975; % Power in eclipse
Qin = 12.6049; % Power in sunlight
Ap_nom = Ar/6; % Nominal absorptive area
h = 600; %km
% surface treatment
% https://www.tiodize.com/emissivity/
% values for Type II TIODIZE PROCESS coatingså
alpha = 0.82; % absorptivity 
eps = 0.51 ; % emissivity

%% Period and Eclipse Calculations
P = round((2*pi*sqrt((Re+h)^3/mu))/60); %Period in minutes
rho = asin(Re/(Re+h));
TE = round(rho/pi*P); % Eclipse time (minutes)
TS = P-TE; % Sunlit time (minutes)

%% Thermal Calculations to determine equilibrium state during nominal ops

To = 293.15; %K - Initial temperature estimate
S_a = S*(1+albedo); % Solar flux with Albedo added. Assumes albedo affects only in sunlight
diff = 1;
step = 0;

while diff > .25
step = step +1; %orbit counter
%Integration over sunlit conditions    
tspan_sun = [0 TS*60];
options = odeset('RelTol',1e-9,'AbsTol',1e-9);
[t_sun,T_sun] = ode45(@(t,T)heat(Qin,Ar,Ap_nom,S_a,T,alpha,eps,m,c),tspan_sun,To,options);
T_sun_val = T_sun(length(T_sun));

% Integration of eclipse conditions
tspan_ecl = [0 TE*60];
[t_ecl,T_ecl] = ode45(@(t,T)heat(Qin_ecl,Ar,Ap_nom,S_ecl,T,alpha,eps,m,c),tspan_ecl,T_sun(length(T_sun)),options);
T_ecl_val = T_ecl(length(T_ecl));

% Difference calculation
diff = abs(T_ecl_val-To);
T_stable = To;
To = T_ecl_val;

% Results plots
t_ecl = t_ecl+TS*60;

sunplot = plot(t_sun/60,T_sun-293,'r','Linewidth',1.5);
grid on, hold on
eclplot = plot(t_ecl/60,T_ecl-293,'b','Linewidth',1.5);
x1 = xline(TS,'k',{'Eclipse Entry'});
x1.LabelVerticalAlignment = 'middle';
x1.LabelHorizontalAlignment = 'center';
xlabel('Time since exiting eclipse, minutes')
ylabel('Equilibrium Temperature, degrees C')
title('Equilibrium Temperature vs Time')
legend({'Sunlit','Eclipse'},'Location','northwest')
hold off
pause(.5)
end

T_ave_in_degC = (sum(T_ecl)+sum(T_sun))/(length(T_ecl)+length(T_sun))-293.15;
T_max_in_degC = T_sun(length(T_sun))-293.15;
T_min_in_degC = T_sun(1)-293.15;

%% Sample Thermal Profile during Sunlit Collision Avoidance Manuever (Peak Power)

%Integration over sunlit conditions    
tspan_sun = [0 (TS-28)*60];
options = odeset('RelTol',1e-9,'AbsTol',1e-9);
[t_sun,T_sun] = ode45(@(t,T)heat(Qin,Ar,Ap_nom,S_a,T,alpha,eps,m,c),tspan_sun,T_sun(1),options);
T_sun_val = T_sun(length(T_sun));

% Integration during propulsion conditions (while sunlit)
Qin_peak = 68.42;
tspan_sunpeak = [0 28*60];
[t_sunpeak,T_sunpeak] = ode45(@(t,T)heat(Qin_peak,Ar,Ap_nom,S_a,T,alpha,eps,m,c),tspan_sunpeak,T_sun(length(T_sun)),options);
T_sun_valpeak = T_sunpeak(length(T_sunpeak));

% Integration of eclipse conditions
tspan_ecl = [0 TE*60];
[t_ecl,T_ecl] = ode45(@(t,T)heat(Qin_ecl,Ar,Ap_nom,S_ecl,T,alpha,eps,m,c),tspan_ecl,T_sunpeak(length(T_sunpeak)),options);
T_ecl_val = T_ecl(length(T_ecl));

% Results plots
t_sunpeak = t_sunpeak+(TS-28)*60;
t_ecl = t_ecl+TS*60;

Propplot = figure();
plot(t_sun/60,T_sun-293,'r','Linewidth',1.5);
grid on, hold on
plot(t_sunpeak/60,T_sunpeak-293,'m','Linewidth',1.5);
plot(t_ecl/60,T_ecl-293,'b','Linewidth',1.5);
x2 = xline(TS,'k',{'Eclipse Entry'});
x2.LabelVerticalAlignment = 'middle';
x2.LabelHorizontalAlignment = 'center';
xlabel('Time since exiting eclipse, minutes')
ylabel('Temperature, degrees C')
title('Temperature Profile with Sunlit Peak Power vs Time')
hold off
legend({'Sunlit','Peak Power','Eclipse'},'Location','northwest')
pause(.5)

T_ave_sunpeak_in_degC = (sum(T_ecl)+sum(T_sun)+sum(T_sunpeak))/(length(T_ecl)+length(T_sun)+length(T_sunpeak))-293.15;
T_max_sunpeak_in_degC = T_sunpeak(length(T_sunpeak))-293.15;
T_min_sunpeak_in_degC = T_sun(1)-293.15;

%% Sample Thermal Profile during Eclipse Collision Avoidance Manuever (Peak Power)

%Integration over sunlit conditions    
tspan_sun = [0 TS*60];
options = odeset('RelTol',1e-9,'AbsTol',1e-9);
[t_sun,T_sun] = ode45(@(t,T)heat(Qin,Ar,Ap_nom,S_a,T,alpha,eps,m,c),tspan_sun,T_sun(1),options);
T_sun_val = T_sun(length(T_sun));

% Integration during propulsion conditions (during eclipse)
Qin_peak = 68.42;
tspan_eclpeak = [0 28*60];
[t_eclpeak,T_eclpeak] = ode45(@(t,T)heat(Qin_peak,Ar,Ap_nom,S_ecl,T,alpha,eps,m,c),tspan_eclpeak,T_sun(length(T_sun)),options);
T_ecl_valpeak = T_eclpeak(length(T_eclpeak));

% Integration of eclipse conditions
tspan_ecl = [0 (TE-28)*60];
[t_ecl,T_ecl] = ode45(@(t,T)heat(Qin_ecl,Ar,Ap_nom,S_ecl,T,alpha,eps,m,c),tspan_ecl,T_eclpeak(length(T_eclpeak)),options);
T_ecl_val = T_ecl(length(T_ecl));

% Results plots
t_eclpeak = t_eclpeak+TS*60;
t_ecl = t_ecl+(TS+28)*60;

Ecl_peak_plot = figure();
plot(t_sun/60,T_sun-293,'r','Linewidth',1.5);
grid on, hold on
plot(t_eclpeak/60,T_eclpeak-293,'m','Linewidth',1.5);
plot(t_ecl/60,T_ecl-293,'b','Linewidth',1.5);
x3 = xline(TS,'k',{'Eclipse Entry'});
x3.LabelVerticalAlignment = 'middle';
x3.LabelHorizontalAlignment = 'center';
xlabel('Time since exiting eclipse, minutes')
ylabel('Temperature, degrees C')
title('Temperature Profile with Eclipse Peak Power vs Time')
hold off
legend({'Sunlit','Peak Power','Eclipse'},'Location','northwest')
pause(.5)

T_ave_eclpeak_in_degC = (sum(T_ecl)+sum(T_sun)+sum(T_eclpeak))/(length(T_ecl)+length(T_sun)+length(T_eclpeak))-293.15;
T_max_eclpeak_in_degC = T_eclpeak(length(T_eclpeak))-293.15;
T_min_eclpeak_in_degC = T_sun(1)-293.15;

