% BK70A0600 Computational Methods in Mechanics
% Assigment 4
% Jesse Myller 0503199

close all; clear all;   % To make sure results from previously run code doesn't impact the results here

% Initializing constants
a = 0.1;    % meters
b = 0.2;    % meters
fii = 30/180*pi;    % radians
omega = 1;  % rad/s

% Forming symbolic expressions
syms x theta d

f_expr = [a.*cos(fii)+b.*cos(theta)-d; a.*sin(fii)-b.*sin(theta)];   % symbolic expression for f(x)
jacob_expr = jacobian(f_expr, [theta; d]);      % compute f''(x) symbolically
% Turn f_expr and jacob_expr into plain Matlab functions
f = matlabFunction(f_expr);
jacob = matlabFunction(jacob_expr);

eps = 1e-6;     % Accuracy for approximation
x = [10, 10];   % Initial guesses (theta, d)
i = 0;          % Iteration count

F_value = f(x(2), x(1));
F_norm = norm(F_value);

% To find valid values for theta and d, we use loop to make use of
% Newton-Raphson method in finding with which values the function F(x) is
% equal to 0
while abs(F_norm) > eps
    delta = jacob(x(1))\-F_value;
    x = x + delta;
    F_value = f(x(2), x(1));
    F_norm = norm(F_value);
    i = i + 1;
end

%% Part 2 Changing of fii

close all; clear all;   % To make sure results from previously run code doesn't impact the results here

% Initializing constants
a = 0.1;    % meters
b = 0.2;    % meters
omega = 1;  % rad/s

% Forming symbolic expressions
syms x theta d fii t

fii = pi/6*omega*t;     % radians

f_expr = [a.*cos(fii)+b.*cos(theta)-d; a.*sin(fii)-b.*sin(theta)];   % symbolic expression for f(x)
jacob_expr = jacobian(f_expr, [theta; d; t]);       % compute f''(x) symbolically
% Turn f_expr and jacob_expr into plain Matlab functions
f = matlabFunction(f_expr);
jacob = matlabFunction(jacob_expr);

eps = 1e-6;     % Accuracy for approximation
x = [10, 10, 10];       % Initial guesses (theta, d, t)
i = 0;          % Iteration count

F_value = f(x(2), x(3), x(1));
F_norm = norm(F_value);

% To find valid values for theta and d, we use loop to make use of
% Newton-Raphson method in finding with which values the function F(x) is
% equal to 0
while abs(F_norm) > eps
    delta = jacob(x(3), x(1))\-F_value;
    x = x + delta;
    F_value = f(x(2), x(3), x(1));
    F_norm = norm(F_value);
    i = i + 1;
end

% Plotting images

t = 0:0.1:10;
% Functions and derivatives for displacement and angle were calculated with
% a calculator. Could very well be completely wrong way to approach this,
% but I was unable to find any other method to make it simpler
d = 0.1.*(sqrt((cos(t+pi./6)).^2+3)-cos(t+pi./6));
derd = (-0.1.*sin(t+pi./6).*(sqrt((cos(t+pi./6)).^2+3)+cos(t+pi./6)))./(sqrt((cos(t+pi./6)).^2+3));
theta = pi-asin(0.5.*sin(t+pi./6));
dertheta = (cos(t+pi./6))/(sqrt(sin(t+pi./6)+2).*sqrt(-1.*(sin(t+pi./6)-2)));

% Time and displacement
figure('Name', 'Time (t) versus displacement (d)')
plot(t,d)
grid on
xlabel('time (s)')
ylabel('displacement (meters)')

% Time and angle
figure('Name', 'Time (t) versus angle (theta)')
plot(t,theta)
grid on
xlabel('time (s)')
ylabel('angle (radians)')

% Time and derivative of displacement
figure('Name', 'Time (t) versus derivative of displacement (d)')
plot(t,derd)
grid on
xlabel('time (s)')
ylabel('time derivative of displacement')

% Time and derivative of angle, this one is not working for some reason,
% perhaps a negative in sqrt but since there is no error or warnings
% related to this, I do not know
figure('Name', 'Time (t) versus derivative of angle (theta)')
plot(t,dertheta)
grid on
xlabel('time (s)')
ylabel('time derivative of angle')