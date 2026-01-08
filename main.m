clear; close all; clc;
addpath('tr_homogenes/');

% Robot PPR 4R R
% MGD
syms l1 l2 l3 l4 l5 l6 l7 l8 theta3 theta4 theta5 theta6 theta7 theta8 real;
L = [l1 l2 l3 l4 l5 l6 l7 l8];
Theta = [theta3 theta4 theta5 theta6 theta7 theta8];
Pi=sym(pi);
    
% world to robot frame
RwT1 = th_rotz(Pi/2) * th_trans(0, 0, 0) * th_rotx(Pi/2);

% prismatic joints - tracks
R1T2 = th_rotz(Pi/2) * th_trans(0, 0, l1) * th_rotx(Pi/2);
R2T3 = th_rotz(-Pi/2) * th_trans(0, 0, l2) * th_rotx(-Pi/2);

% rotoid on z axis
R3T4 = th_rotz(theta3) * th_trans(0, 0, l3) * th_rotx(-Pi/2);

% 4R on xy plane
R4T5 = th_rotz(theta4) * th_trans(-l4, 0, 0);
R5T6 = th_rotz(theta5) * th_trans(-l5, 0, 0);
R6T7 = th_rotz(theta6) * th_trans(-l6, 0, 0);
R7T8= th_rotz(theta7 + Pi/2) * th_trans(l7, 0, 0) *th_rotx(-Pi/2);

% camera roll rotoid
R8T9 = th_rotz(theta8) * th_trans(0, 0, l8);


% world -> end-effector
RwT9 = RwT1 * R1T2 * R2T3 * R3T4 * R4T5 * R5T6 * R6T7 * R7T8 * R8T9;
RwT9(1:3, 4)

RwT9_S = subs(RwT9, [theta3 theta4 theta5 theta6 theta7 theta8], [0 0 0 0 0 0]);
RwT9_S(1:3, 4)  % expected position of end-effector


RwPR9 = RwT9(1:3, 4);
matlabFunction(RwPR9, 'File', 'clc_mgd.m', 'Vars', {[L Theta]});

%{
% Newton's method MGI
% Newton Raphson method : Solve f(x) = 0
% We want the end-effector to be at a specific location 
% P_target = [xp; yp; zp]

% The direct geometric model is located in the last column of the 
% homogeneous transform R1T5
RwPR9 = RwT9(1:3, 4);

% We want to have the end-effector point at a specific location
% in the world frame
syms xp yp zp real;
P_target = [xp; yp; zp]; % world frame

% This means solving
%     RwPR9 = P_target
% <=> RwPR9- P_target = 0
f = RwPR9 - P_target
J = jacobian(f, [theta3 theta4 theta5 theta6 theta7])
size(J)
% By directly computing the inverse (computationally inefficient)
% J_inv = pinv(J)

% Newton-Raphson scheme
% Xn = Xn-1 - Jf(Xn-1)^-1 

%%
l = [0 0 1 1 1 1 1 0];
f = subs(f, [l1 l2 l3 l4 l5 l6 l7 l8 Pi], [l pi]);
J_inv = subs(J_inv, [l1 l2 l3 l4 l5 l6 l7 l8 Pi], [l pi]);

% using values for MGD above
S_P_target = [1, 0, 3];
f = subs(f, [xp yp zp], S_P_target);
J_inv = subs(J_inv, [xp yp zp], S_P_target);

X = [pi; pi; pi; pi; pi];  % theta angles

disp('start newton-raphson method')
for i = 1:6
    % replace with matlabFunction call (faster)
    S_J_inv = subs(J_inv, [theta3 theta4 theta5 theta6 theta7], X);
    S_f = subs(f, [theta3 theta4 theta5 theta6 theta7], X);
    X = X - S_J_inv * S_f;
end
disp('end');
%}

% Newton-Raphson + Moving inside the null-space of J
% get direct kinematic model for each rotoid joint
RwT4 = RwT1 * R1T2 * R2T3 * R3T4;
RwT5 = RwT1 * R1T2 * R2T3 * R3T4 * R4T5;
RwT6 = RwT1 * R1T2 * R2T3 * R3T4 * R4T5 * R5T6;
RwT7 = RwT1 * R1T2 * R2T3 * R3T4 * R4T5 * R5T6 * R6T7;

% generate function to retrieve position of each rotoid joint on the plane
r_poses = [RwT4(1:3, 4) RwT5(1:3, 4) RwT6(1:3, 4) RwT7(1:3, 4)]; % only position required

% define obstacle to avoid
syms xo yo zo real;
Obs_P = [xo;yo;zo];

% compute distances & their gradient for each joint in 4R part
n = length(r_poses);
dist_q = [];
grad_dist_q = {[]}; % ugh not very pretty
for i = 1:n
    dist_q = [dist_q norm(Obs_P - r_poses(:, i))];
    grad_dist_q{i} = gradient(dist_q(i), Theta);
end
grad_dist_q{1}

matlabFunction(dist_q, 'File', 'clc_dist.m', 'Vars', {[L Theta Obs_P']});
matlabFunction( ...
    grad_dist_q{1}, grad_dist_q{2}, grad_dist_q{3}, grad_dist_q{4}, ...
    'File', 'clc_gradient_dist.m', 'Vars', {[L Theta Obs_P']});


function Z = clc_Z(L, q, Obs_P, Obs_r)
    % L : lengths vector
    % q : joint angles
    % Obs_p : 1x3 row vector of obstacle center
    % Obs_r : Radius of obstacle
    
    % get joint with closest distance to obstacle
    dist_q = clc_dist([L q Obs_P]);
    [v, i] = min(dist_q);
    
    % Check that joint is inside obstacle radius
    if v < 2 * Obs_r
        % we'll use this joint's gradient vector as optimization
        % to move it away from the obstacle
        [gdq1, gdq2, gdq3, gdq4] = clc_gradient_dist([L, q, Obs_P])
        gdq = {gdq1 gdq2 gdq3 gdq4}
        Z = gdq{i};
    else
        Z = zeros([length(q) 1]);
    end
end

% define target position X
P_target = [2.5; 0; 1]; % column vector
P_end_effector = RwPR9;

f = P_end_effector - P_target;
matlabFunction(f, 'File', 'clc_f.m', 'Vars', {[L Theta]});
J_f = jacobian(f, [theta3 theta4 theta5 theta6 theta7 theta8]);
matlabFunction(J_f, 'File', 'clc_jacobian_f.m', 'Vars', {[L Theta]});

%% Newton-Raphson scheme

% only 4R plane robot has lengths of 1
% and joint 3 to put it at z-height of 1
Lv = [0 0 1 1 1 1 1 0];
q0 = [0; 0; 0; 0; 0; 0];
In = eye(length(q0));
q = q0;

% obstacle O
Obs_P = [1.5; 0; 1.5];
Obs_r = 0.4;

max_iter = 900;
for step = 1:max_iter
    % this is messier than my own room
    % i'm open to suggestion to prettify this
    Z = clc_Z(Lv, q', Obs_P', Obs_r);

    f = clc_f([Lv q']);
    J_f = clc_jacobian_f([Lv q']);
    pinv_J_f = clc_pinv(J_f);

    q = q - pinv_J_f * f + (In - pinv_J_f * J_f) * Z;
end

q

%% draw robot
Plist = [];
Tr = {RwT1 R1T2 R2T3 R3T4 R4T5 R5T6 R6T7 R7T8 R8T9};
cur_tr = 1;
for i = 1:length(Tr)
    cur_tr = cur_tr * Tr{i};
    P_tr = cur_tr(1:3, 4);
    P_end_effector = subs(P_tr, [L Theta], [Lv q']);
    Plist = [Plist P_end_effector];
end
double(Plist)

axis equal;
view(3);
hold on;
grid on;

xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('Inverse kinematics of 4R component of robot avoiding obstacle')
plot3(Plist(1, :), Plist(2, :), Plist(3, :));

[A,B,C] = sphere
A = A * Obs_r;
B = B * Obs_r;
C = C * Obs_r;
surf(A + Obs_P(1),B + Obs_P(2), C + Obs_P(3), 20, "FaceColor", "black")