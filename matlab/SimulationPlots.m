% Run this script to plot the simulations

% Load the workspace
load("workspace.mat");


% #1 --> Altitude (r)
figure('DefaultAxesFontSize',13);

plot(tout,zout(:,1)-rmoon,'linewidth',4);
title("Altitude",'Interpreter','latex','fontsize',24);
grid;

xlabel("Time [s]",'fontsize',15);
ylabel("Altitude [km]",'fontsize',15);
ylim([0,15]);
    
% #2 --> Horizontal velocity (u)
figure('DefaultAxesFontSize',13);

plot(tout,zout(:,3),'linewidth',4);
title("The horizontal velocity",'Interpreter','latex','fontsize',24);
grid;

xlabel("Time [s]",'fontsize',15);
ylabel("Horizontal velocity [km/s]",'fontsize',15);

% #3 --> Vertical velocity (v)
figure('DefaultAxesFontSize',13);

plot(tout,zout(:,4),'linewidth',4);
title("The vertical velocity",'Interpreter','latex','fontsize',24);
grid;

xlabel("Time [s]",'fontsize',15);
ylabel("Vertical velocity [km/s]",'fontsize',15);

% #4 --> Control input (beta)
figure('DefaultAxesFontSize',13);

plot(tout,rad2deg(atan2(-zout(:,9),-zout(:,8))),'linewidth',4);
title("Control input",'Interpreter','latex','fontsize',24);
grid;

xlabel("Time [s]",'fontsize',15);
ylabel("Control input [deg.]",'fontsize',15);


% #5 --> Mass
figure('DefaultAxesFontSize',13);

plot(tout,zout(:,5),'linewidth',4);
title("Mass",'Interpreter','latex','fontsize',24);
grid;

xlabel("Time [s]",'fontsize',15);
ylabel("Control input [deg.]",'fontsize',15);
