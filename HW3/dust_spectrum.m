function [Teq,BQ] = dust_spectrum(Pin,R_star,r_orb,R_grain,fnum)
%%
% dust_spectrum
%
% Calculate the spectrum and equilibrium temperature from the dust ring
% given grain size and power input.
%
% INPUT
% Pin - [double] power input [J/s]
% R_star - [double] body radius of star [R_sun]
% r_orb  - [double] orbital radius of secondary object [AU]
% R_grain - [double] dust grain size [micron]
% fnum - [int] figure number (0 if not plotting)
%
% OUTPUT
% Teq - [double] equilibrium temperature [K]
% Fv_out - [double] spectrum flux density [W/m2-Hz]
% 
% AUTHOR
% Chun-Yi Wu
global Rdebs lambdaQs Qabs vQs
R_sun = 6.957e8;   % solar radius [m]
AU = 1.496e+11; % astronomical unit [m]

%% Calculation
% determine with set of Qabs values to use
lambdaQi = 22;
for ( i = 1 : 21 )
    if ( R_grain*1e-6 <= Rdebs(i) )
        lambdaQi = i;
        break;
    end
end

if ( lambdaQi < 22 )
    Qvs = Qabs(:,lambdaQi);
else
    Qvs = ones(size(vQs));
end

% solve for T_dust
f = @(T) power_difference(T,R_grain,1/R_sun,1/AU,Pin,Qvs);
T_guess = 10;
Teq = fzero(f,T_guess);

k = 1.38064852e-23; % Boltzmann constant [J/K]
h = 6.62607004e-34; % Planck constant [J-s]
c = 2.99792458e+08; % speed of light in vacuum [m/s]
% Fv_out = flux_density_mks(1/R_sun,Teq,1/AU,lambdaQs,0) .* Qvs;
Bvs = (2*h*vQs.^3/c^2) ./ ( exp(h*vQs/(k*Teq)) - 1 );

% multiply by absorption ratio
BQ = Bvs .* Qvs;
if ( fnum > 0 )
    figure(fnum); clf(fnum);
    semilogx(lambdaQs*1e6,BQ*1e26,'-');
    grid on; hold on;
    xlabel('wavelength [\mum]');
    ylabel('flux density [Jy]');
    title('dust ring output');
    axis([1,1000,0,max(BQ*1e26)*1.1]);
end
end


function dP = power_difference(T_dust,R_grain,R_star,r_orb,P_in,Qvs)
%%
% power_difference
%
% The power difference between two model. Zeroing this function means the
% power matches.
%
% INPUT
% T_dust - [double] temperature
% R_star - [double] body radius of star [R_sun]
% r_orb  - [double] orbital radius of secondary object [AU]
% R_grain - [double] dust grain size [m]
% P_in - [double] power absorbed
% Qvs - [double] absorption ratio
%
% OUTPUT
% dP - [double] difference in power
global Rdebs lambdaQs Qabs vQs

k = 1.38064852e-23; % Boltzmann constant [J/K]
h = 6.62607004e-34; % Planck constant [J-s]
c = 2.99792458e+08; % speed of light in vacuum [m/s]

%% Calculation
% get flux density
% vs = c ./ ls;  % frequency

Bvs = (2*h*vQs.^3/c^2) ./ ( exp(h*vQs/(k*T_dust)) - 1 );

% multiply by absorption ratio
BQ = Bvs .* Qvs;

% find total power
P_dust = 4 * pi^2 * (R_grain*1e-6)^2 * trapz(vQs,BQ);

% find difference
dP = (P_dust - P_in);
end