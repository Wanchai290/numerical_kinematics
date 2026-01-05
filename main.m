clear; close all; clc;
addpath('tr_homogenes/');

% Robot PPRR
% MGD
syms l1 l2 l3 l4 l5 l6 l7 theta3 theta4 theta5 theta6 theta7 real;
Pi=sym(pi);

RwT1 = th_rotz(Pi/2) * th_trans(0, 0, 0) * th_rotx(Pi/2);

R1T2 = th_rotz(Pi/2) * th_trans(0, 0, l1) * th_rotx(Pi/2);
R2T3 = th_rotz(-Pi/2) * th_trans(0, 0, l2) * th_rotx(-Pi/2);
R3T4 = th_rotz(theta3) * th_trans(0, 0, l3) * th_rotx(-Pi/2);
R4T5 = th_rotz(theta4) * th_trans(-l4, 0, 0);
R5T6 = th_rotz(theta5) * th_trans(-l5, 0, 0);
R6T7 = th_rotz(theta6 + Pi/2) * th_trans(l6, 0, 0) *th_rotx(-Pi/2);
R7T8 = th_rotz(theta7);

RwT8 = RwT1 * R1T2 * R2T3 * R3T4 * R4T5 * R5T6 * R6T7 * R7T8;
RwT8(1:3, 4)

RwT8_S = subs(RwT8, [theta3 theta4 theta5 theta6 theta7], [0 0 0 0 0]);
RwT8_S(1:3, 4)  % expected position of end-effector

%% Newton's method MGI
% Newton Raphson method : Solve f(x) = 0
% Suppose the world frame is R1
% We want the end-effector to be at a specific location 
% P_target = [xp; yp; zp]

% The direct geometric model is located in the last column of the 
% homogeneous transform R1T5
RwPR8 = RwT8(1:3, 4);

% We want to have the end-effector point at a specific location
% in the world frame
syms xp yp zp real;
P_target = [xp; yp; zp]; % world frame

% This means solving
%     RwPR8 = P_target
% <=> RwPR8- P_target = 0
f = RwPR8 - P_target
J = jacobian(f, [theta3 theta4 theta5 theta6 theta7])

% By directly computing the inverse (computationally inefficient)
J_inv = pinv(J)

% Newton-Raphson theme
% Xn = Xn-1 - Jf(Xn-1)^-1 

%%
l = [0 0 1 1 1 1 1];
f = subs(f, [l1 l2 l3 l4 l5 l6 l7 Pi], [l pi]);
J_inv = subs(J_inv, [l1 l2 l3 l4 l5 l6 l7 Pi], [l pi]);

% using values for MGD above
S_P_target = [1, 0, 3];
f = subs(f, [xp yp zp], S_P_target);
J_inv = subs(J_inv, [xp yp zp], S_P_target);

%f = subs(f, [xp yp zp], [-1, 1, 1]);
%J_inv = subs(J_inv, [xp yp zp], [-1, 1, 1]);

X = [pi; pi];  % theta angles

disp('start newton-raphson method')
for i = 1:6
    disp(double(X))
    disp(i)

    S_J_inv = subs(J_inv, [theta3; theta4], X);
    S_f = subs(f, [theta3; theta4], X);
    X = X - S_J_inv * S_f;
end
disp('end');