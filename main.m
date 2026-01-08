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
matlabFunction(RwPR9, 'File', 'clc_mgd.m', 'Vars', [L Theta]);

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

%% Newton-Raphson + Moving inside the null-space of J
% get direct kinematic model for each rotoid joint
RwT4 = RwT1 * R1T2 * R2T3 * R3T4;
RwT5 = RwT1 * R1T2 * R2T3 * R3T4 * R4T5;
RwT6 = RwT1 * R1T2 * R2T3 * R3T4 * R4T5 * R5T6;
RwT7 = RwT1 * R1T2 * R2T3 * R3T4 * R4T5 * R5T6 * R6T7;

% generate function to retrieve position of each rotoid joint on the plane
r_poses = simplify([RwT4(1:3, 4) RwT5(1:3, 4) RwT6(1:3, 4) RwT7(1:3, 4)]); % only position required

% define obstacle to avoid
syms xo yo zo Obs_r real;
Obs_P = [xo;yo;zo];

% L2 distance vectorized function
function dist = clc_dist(A, B, r)
    % size(A) = 3 1
    d = sqrt(A - B) .^2;
    dist(d <= r) = d;
    dist(d > r) = 0;
end

% create Phi function for vector Z
Phi_f = 0; % scalar function
for i = 1:length(r_poses)
    d = clc_dist(Obs_P, r_poses(1:3, i), Obs_r);
    if d ~= 0
        Phi_f = Phi_f + 1 / d;
    end
end

% Compute the gradient of the Phi function with respect to the joint positions
grad_Phi_f = gradient(Phi_f, [theta3 theta4 theta5 theta6 theta7 theta8])
Z = grad_Phi_f;

matlabFunction(Z, 'File', 'clc_Z.m', 'Vars', [L Theta Obs_P']);
%%
% define target position X
P_target = [2; 0; 3]; % column vector
P_end_effector = RwPR9;

f = P_end_effector - P_target;
matlabFunction(f, 'File', 'clc_f.m', 'Vars', [L Theta]);
J_f = jacobian(f, [theta3 theta4 theta5 theta6 theta7 theta8]);
matlabFunction(J_f, 'File', 'clc_jacobian_f.m', 'Vars', [L Theta]);
%%
% Newton-Raphson scheme
clear;

% only 4R plane robot has lengths of 1
% and joint 3 to put it at z-height of 1
L = [0 0 1 1 1 1 1 0];
q0 = [0; 0; 0; 0; 0; 0];
In = eye(length(q0));
q = q0;

% obstacle O
Obs_P = [1.5; 0 ;1];


max_iter = 200;
for step = 1:max_iter
    % this is messier than my own room
    % i'm open to suggestion to prettify this
    f = clc_f(L(1), L(2), L(3), L(4), L(5), L(6), L(7), L(8), ...
    q(1), q(2), q(3), q(4), q(5), q(6));
    J_f = clc_jacobian_f(L(1), L(2), L(3), L(4), L(5), L(6), L(7), L(8), ...
    q(1), q(2), q(3), q(4), q(5), q(6));
    pinv_J_f = clc_pinv(J_f);
    Z = clc_Z(L(1), L(2), L(3), L(4), L(5), L(6), L(7), L(8), ...
        q(1), q(2), q(3), q(4), q(5), q(6), Obs_P(1), Obs_P(2), Obs_P(3))
    
    q = q - pinv_J_f * f + (In - pinv_J_f * J_f) * Z;
end
function x = clipZeroOne(x)
x(x < 0) = 0;
x(x > 1) = -1;
end
q
