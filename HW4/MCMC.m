function [xss] = MCMC(xbar,sigbar,xtyp,fxn,nwalker,itermax)
%%
% MCMC
%
% Perform estimation using Markov Chain Monte Carlo algorithm.
%
% INPUT
% xbar - [double] initial guess
% sigbar - [double] initial guess error
% xtyp - [double] typical value for parameters (used to scale step sizes)
% fxn - [function handle] the function that produces chi-square score
% nwalker - [integer] number of walker present
% itermax - [integer] number of iterations each walker take
%
% OUTPUT
% xss - [double] all points
%
% AUTHOR
% Chun-Yi Wu

%% constants/parameters
nx = length(xbar);
tol = 0.2;
ratio = ones(nx);

njumps = zeros(nx,1);
naccept = zeros(nx,1);
           
load('data.mat');
t = data(:,1);
z = data(:,2);
R = data(:,3);
            figure(5); clf(5); 

%% algorithm
xss_cell = cell(nwalker,1);
xss = zeros(nx,itermax*nwalker);

for ( iwalker = 1 : nwalker )
    % initiate walker
    if ( iwalker == 1 )
        x0 = xbar;
    else
        x0 = xbar + randn(1,nx) .* sigbar;
    end
    chi20 = fxn(x0);
    
    % burn-in
    for ( iter = 1 : itermax/5 )
        if ( mod(iter,100) == 0 )
            fprintf('burn: %d/%d\n',iter,itermax);
            for ( i = 1 : nx )
                fprintf('\t%g %g\n',x0(i),ratio(i)*xtyp(i));
            end
            fprintf('\n');
            racc = naccept ./ njumps;
            
            for ( i = 1 : nx )
                if ( racc(i) < 0.2 )
                    ratio(ind) = ratio(ind) * 0.9;
                elseif ( racc(i) > 0.4 )
                    ratio(ind) = ratio(ind) * 1.1;
                end
            end
            

            plot(mod((t-x0(1)),x0(2))/x0(2),z,'+');
            drawnow;
            njumps = zeros(nx,1);
            naccept = zeros(nx,1);
        end
        new_x_accepted = false;
        
        while ( ~new_x_accepted )
            x1 = x0;
%             x1 = x1 + randn(nx,1) .* xtyp * ratio;
            ind = randi(nx);
            x1(ind) = x1(ind) + randn() * xtyp(ind) * ratio(ind);
            njumps(ind) = njumps(ind) + 1;
            chi21 = fxn(x1);

            if ( rand() < min(1,exp(-(chi21-chi20)/2)) )
                x0 = x1;
                chi20 = chi21;
                new_x_accepted = true;
                xs(:,iter) = x0;
                naccept(ind) = naccept(ind) + 1;
            end
        end
    end

    fprintf('|');
    
    % actual minimum finding
    xs = zeros(nx,itermax);
    for ( iter = 1 : itermax )
        if ( mod(iter,100) == 0 )
            fprintf('%d/%d \n',iter,itermax);
            for ( i = 1 : nx )
                fprintf('\t%g %g\n',x0(i),ratio(i)*xtyp(i));
            end
            racc = naccept ./ njumps;
            
            for ( i = 1 : nx )
                if ( racc(i) < 0.2 )
                    ratio(ind) = ratio(ind) * 0.9;
                elseif ( racc(i) > 0.4 )
                    ratio(ind) = ratio(ind) * 1.1;
                end
            end
            
            plot(mod((t-x0(1)),x0(2))/x0(2),z,'+');
            drawnow;
            njumps = zeros(nx,1);
            naccept = zeros(nx,1);
        end
        new_x_accepted = false;
        
        while ( ~new_x_accepted )
            x1 = x0;
%             x1 = x1 + randn(nx,1) .* xtyp * ratio;
            ind = randi(nx);
            x1(ind) = x1(ind) + randn() * xtyp(ind) * ratio(ind);
            njumps(ind) = njumps(ind) + 1;
            chi21 = fxn(x1);

            if ( rand() < min(1,exp(-(chi21-chi20)/2)) )
                x0 = x1;
                chi20 = chi21;
                new_x_accepted = true;
                xs(:,iter) = x0;
                naccept(ind) = naccept(ind) + 1;
            end
        end
    end
        
        
    xss_cell{iwalker} = xs;
end

for ( iwalker = 1 : nwalker )
    xss(:,(iwalker-1)*itermax+(1:itermax)) = xss_cell{iwalker};
end
end