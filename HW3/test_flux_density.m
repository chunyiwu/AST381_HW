clear; clc; close all hidden;

readQtable()

R_star = 1.842;
T_star = 8590;
% T_star = 5000;
r_orb  = 10;
lambdas = linspace(1e-10,3.5,1001) * 1e-6;
R_grain = 1e-1;

fnum = 1;
[Fvs] = flux_density_mks(R_star,T_star,r_orb,lambdas,fnum);

fnum = 2;
Pin = total_power_in(R_star,T_star,r_orb,R_grain,fnum)

fnum = 3;
Teq = dust_spectrum(Pin,R_star,r_orb,R_grain,fnum)