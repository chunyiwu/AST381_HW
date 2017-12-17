function [Fvs] = flux_density_mks(R_star,T_star,r_orb,ls,fnum)
%%
% flux_density
%
% Calculate and plot the flux-calibrated spectrum of a star from certain
% distance.
%
% The program is coded with mks unit system.
%
% INPUT
% R_star - [double] body radius of star [R_sun]
% T_star - [double] temperature of star [K]
% r_orb  - [double] orbital radius of secondary object [AU]
% ls     - [double] wavelengths to calculate [m]
% fnum   - [int] figure number (0 if not plotting)
%
% OUTPUT
% Fvs    - [double] fluxes at each specified wavelength
%
% AUTHOR
% Chun-Yi Wu

%% Constants/parameters
k = 1.38064852e-23; % Boltzmann constant [J/K]
h = 6.62607004e-34; % Planck constant [J-s]
c = 2.99792458e+08; % speed of light in vacuum [m/s]

R_sun = 6.957e8;   % solar radius [m]
AU = 1.496e+11; % astronomical unit [m]

%% Convert inputs to SI units
R_star = R_star * R_sun;
r_orb = r_orb * AU;

%% Find spectral flux at surface of star
vs = c ./ ls;  % frequency

Bvs = (2*h*vs.^3/c^2) ./ ( exp(h*vs/(k*T_star)) - 1 );
% Bvs = 2*h*c^2*ls.^-5 ./ ( exp(h*c./(ls*k*T_star)) - 1 );
Fvs = Bvs * (R_star/r_orb)^2 * pi ;

%% Plotting (if requested)
if ( fnum > 0 )
    figure(fnum); clf(fnum); hold on;
%     plot(ls*1e6,Bvs*1e26*1000,'-');
    plot(ls*1e6,Fvs*1e26,'-');
%     plot(vs,Bvs*1e29,'-');
    grid on; hold on;
    xlabel('wavelength [\mum]');
%     xlabel('freq [Hz]');
    ylabel('flux density [Jy]');
    title('dust ring input');
end
end