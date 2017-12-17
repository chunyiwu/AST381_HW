clear; clc; close all;

load('data.mat');
t = data(:,1);
z = data(:,2);
R = data(:,3);

xbar = [2452854.825-2450000;3.5247;0.71*sind(86.1)];
% xbar = [2900;3.52;0.7];
% xbar = [3000;1500;1];

J = costfxn_plot(xbar,2)
figure(5); clf(5); hold on;
plot(mod((t-xbar(1)),xbar(2))/xbar(2),z,'+');
drawnow;

% [xss] = MCMC(xbar,[5;0.5;0.5],[1e-2;1e-4;1e-4],@costfxn,1,5e5);

load('run_data_20171216T212319.mat');

figure(1); clf(1); hold on;
plot3(xss(1,:),xss(2,:),xss(3,:),'.');
xlabel('t_0 - 2450000 [MJD]');
ylabel('P [days]');
zlabel('M_2sini [M_{Jup}]');

xmed = median(xss,2)
Jmed = costfxn_plot(xmed,3)

% xmin = fminsearch(@costfxn,xbar)
% Jmin = costfxn_plot(xmin,3)
% 
% nxss = size(xss,2);
% dJ = zeros(2,nxss);
% dJ(2,:) = 1:nxss;
% 
% for ( i = 1 : nxss )
%     dJ(1,i) = costfxn(xss(:,i)) - Jmed;
% end
% 
% dJ = sortrows(dJ',1);

ps = [0.68, 0.95, 0.997];

for ( i = 1 : 3 ) 
    p = ps(i)

    xss_p = xss(:,dJ(1:floor(p*nxss),2))';
    x_lo = min(xss_p)
    x_hi = max(xss_p)
end

figure(5); clf(5); hold on;
errorbar(mod((t-xbar(1)),xbar(2))/xbar(2),z,R,'+');
drawnow;

% save(sprintf('run_data_%s.mat',datestr(now(),30)));

figure(101); clf; hold on;
plot(xss(1,:),xss(2,:),'.');
plot(xmed(1),xmed(2),'r.');

