clear; close all; clc;
syms x1 x2 real;

% Define two separate functions
f1 = x1^2 + 5*x1 + 3;
f2 = -x2^2 - 3*x2 - 1;
f = [f1;f2];

% Also define their derivative
df1 = 2 * x1 + 5;
df2 = -2 * x2 - 3;
J = [df1, 0;
     0, df2];

hold on;
fplot(f1, Color="red");
fplot(f2, Color="blue");
fplot(0, Color="black");
legend(Location="best");

%% Newton-Raphson method
X0 = [12; 6];
X = X0;
J_inv = pinv(J);
for i = 1:50
    SJ_inv = subs(J_inv, [x1;x2], X);
    Sf = subs(f, [x1; x2], X);
    X = X - SJ_inv * Sf;
end
double(X)

