% Collision Avoidance Phasing Maneuver
clear all, close all, clc
set(0,'defaulttextinterpreter','latex')

% Starting at relative position (0,0,0) km with no initial relative
% velocity, rendezvous to position advanced 10km in the in-track direction
% in 24 hours

mu = 3.986004415*10^5; % Earth's gravitational constant, km^3/s^2
Re = 6378;
a = Re+600;
n = sqrt(mu/(a^3));

% Set Initial Conditions
r0 = [0, 0, 0]';
v0 = [0, 0, 0]';

% Desired final position
rf = [0, 10, 0]';

% Define STM submatrices
t = 24*3600; % 24 hours converted to seconds
c = cos(n*t);
s = sin(n*t);
PHI_rr = [4-3*c, 0, 0;6*(s-n*t), 1, 0;0, 0, c];
PHI_rv = [1/n*s, 2/n*(1-c), 0;-2/n*(1-c), 1/n*(4*s-3*n*t), 0;0, 0, s/n];
PHI_vr = [3*n*s, 0, 0;-6*n*(1-c), 0, 0;0, 0, -n*s];
PHI_vv = [c, 2*s, 0;-2*s, 4*c-3, 0;0, 0, c];

% Calculate required Delta v's
Dv1 = inv(PHI_rv)*(rf - PHI_rr*r0) - v0
Dv2 = -(PHI_vr*r0 + PHI_vv*(v0 + Dv1))

% Visualize trajectory
tspan = 0:10:t;
r = NaN(3,length(tspan));
v = NaN(3,length(tspan));
for i = 1:length(tspan)
    tt = tspan(i);
    x0 = [r0; v0+Dv1];
    xt = hcwstm(n,x0,tt);
    r(:,i) = xt(1:3);
    v(:,i) = xt(4:6);
    if i == length(tspan)
        vt = v(:,i) + Dv2;
        v(:,i) = vt;
    end
end


Dv_tot = norm(Dv1)+norm(Dv2)

figure
plot(r(2,:),r(1,:),'Linewidth',2), grid on, hold on
plot(r(2,1),r(1,1),'bo','MarkerFaceColor','b')
xlabel('In-track, km','FontSize',14)
ylabel('Radial, km','FontSize',14)

%% CW Relative Motion State Transition Matrix
function [xfinal] = hcwstm(n_param, x0initial, delta_t)
c = cos(n_param*delta_t);
s = sin(n_param*delta_t);
stm = [4 - 3*c, 0, 0, s/n_param, 2/n_param*(1-c), 0;
       6*(s - n_param*delta_t), 1, 0, -2/n_param*(1-c), (4*s - 3*n_param*delta_t)/n_param, 0;
       0, 0, c, 0, 0, s/n_param;
       3*n_param*s, 0, 0, c, 2*s, 0;
       -6*n_param*(1-c), 0, 0, -2*s, 4*c-3, 0;
       0, 0, -n_param*s, 0, 0, c];
xfinal = stm*x0initial;
end