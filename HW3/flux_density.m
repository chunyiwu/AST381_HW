function [Fvs] = flux_density(R_star,T_star,r_orb,lambdas,fnum)
%%
% flux_density
%
% Calculate and plot the flux-calibrated spectrum of a star from certain
% distance.
%
% The program is coded with cgs unit system.
%
% INPUT
% R_star - [double] body radius of star [R_sun]
% T_star - [double] temperature of star [K]
% r_orb  - [double] orbital radius of secondary object [AU]
% lambdas - [double] wavelengths to calculate [micron]
% fnum   - [int] figure number (0 if not plotting)
%
% OUTPUT
% Fvs    - [double] fluxes at each specified wavelength
%
% AUTHOR
% Chun-Yi Wu

%% Constants/parameters
k = 1.38064852e-16; % Boltzmann constant [erg/K]
h = 6.62607004e-27; % Planck constant [erg-s]
c = 2.99792458e+10; % speed of light in vacuum [cm/s]

R_sun = 6.957e8;   % solar radius [m]
AU = 1.496e+11; % astronomical unit [m]

%% Convert inputs to SI units
R_star = R_star * R_sun;
r_orb = r_orb * AU;
lambdas = lambdas * 1e-4;

%% Find spectral flux at surface of star
vs = c ./ lambdas;

Bvs = (2*h*vs.^3/c^2) ./ ( exp(h*vs/(k*T_star)) - 1);
Fvs = Bvs * (R_star/r_orb)^2 * pi;

%% Plotting (if requested)
if ( fnum > 0 )
    figure(fnum); clf(fnum);
    plot(lambdas*1e4,Fvs,'-');
    grid on; hold on;
    xlabel('wavelength [\mum]');
    ylabel('flux density [Jy]');
end
end