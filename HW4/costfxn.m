function J = costfxn(x)
%%
% costfxn
%
% The cost function to solve. This is for solving the HD 209458 planet
% orbit fitting, for AST 381 homework 4.
%
% INPUT
% t - [nz double] time of measurements
% z - [nz double] measurements
% R - [nz double] measurements error
% x - [nx double] parameters (t0, P, M2sini)
%
% OUTPUT
% J - [double] cost function

load('data.mat');
t = data(:,1);
z = data(:,2);
R = data(:,3);

M1 = 1.13 * 1.9855e30; % stellar mass [kg]
t0 = x(1)+2450000;
P = x(2)*86400; 
M2sini = x(3)*1.898e27;
a = (6.67408e-11*M1/(4*pi^2)*P^2)^(1/3);

J = 0;
coeff = 2*pi/P * M2sini/M1 * a ;
for ( i = 1 : length(t) )
    zbar = coeff * sin(2*pi*(t(i)-t0)*86400/P);
    J = J + ((zbar-z(i))/R(i))^2;
end

% J = J / length(t);






end