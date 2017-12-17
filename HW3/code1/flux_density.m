function [] = flux_density(R_star,T_star,r_orb)
%%
% flux_density
%
% Plot the flux-calibrated spectrum of a star from certain distance.
%
% INPUT
% R_star - [double] body radius of star [R_sun]
% T_star - [double] temperature of star [K]
% r_orb  - [double] orbital radius of secondary object [AU]
%
% AUTHOR
% Chun-Yi Wu

%% Constants/parameters
kB = 1.38064852e-23;    % Boltzmann constant [J/K]
h  = 6.626070040e-34;   % Planck constant [J-s]
c  = 299792458;         % speed of light in vacuum [m/s]

Jy = 1e-26;             % Jansky conversion ratio [W/m^2-Hz]
R_sun = 6.957e8;        % solar radius [m]
AU = 149597870700;      % astronomical unit [m]

%% Convert inputs to SI units
R_star = R_star * R_sun;
r_orb = r_orb * AU;

%% Find spectral flux at surface of star
vs = linspace(1e4,1e16,50);
Bvs = zeros(size(vs));

for ( i = 1 : 50 )
    % black body radiation
    Bvs(i) = 2*h*vs(i)^3/c^2 / ( exp(h*vs(i)/(kB*T_star)) - 1 );
end

%% Find spectral flux at orbital distance
% Bvs = Bvs * ( R_star / r_orb )^2 * pi;


%% Plot result
figure; clf; grid on;
loglog(vs,Bvs/Jy,'-');
xlabel('frequency [Hz]');
ylabel('flux density [Jy]');
end