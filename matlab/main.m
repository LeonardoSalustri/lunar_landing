
clear; clc;


% Fixed parameters
rmoon = 1737.1; %km
Isp = 310; %s
M = 300; %kg
V0 = 1.69196926; %km/s
T = 0.440; 
g = 9.81e-3;
omega = 2.7e-6;
miu = 4902.7779;

% Initial guesses
tf = 300; %s
r0 = 15+rmoon; 
phi0 = pi; 
u0 = randn(1,1);
v0 = 0;
m0 = M*exp((u0-V0)/(Isp*g));

lambda_rand = randn(1,4);
lambda_r0 = lambda_rand(1);
lambda_phi0 = lambda_rand(2);
lambda_u0 = lambda_rand(3);
lambda_v0 = lambda_rand(4);
lambda_m0 = -lambda_u0/(M/(Isp*g)*exp((u0-V0)/(Isp*g)));

% Initial guesses
guess = [tf; r0; phi0; u0; v0; m0; lambda_r0; lambda_phi0; lambda_u0; lambda_v0; lambda_m0];


options = odeset('RelTol',1e-8,'AbsTol',[1e-8, 1e-8,1e-8,1e-8,1e-8,1e-8,1e-8,1e-8,1e-8,1e-8]);

% Solve the differential augumented state equations 
[tout,zout] = ode45(@augmented_dynamics,[0 guess(1)],guess(2:end),options);

% Final augumented state's parameters
rf = zout(end,1); 
phif = zout(end,2); 
uf = zout(end,3); 
vf = zout(end,4); 
mf = zout(end,5);

lambda_rf = zout(end,6); 
lambda_phif = zout(end,7); 
lambda_uf = zout(end,8); 
lambda_vf = zout(end,9); 
lambda_mf = zout(end,10);

% Control input
beta_f = atan2(-lambda_vf,(-lambda_uf));

% Computing the final value for the differentiated state, costate and 
% doubled differentiated state 
dstatef = dx(zout(end,1:5),beta_f);
dlambdaf = dcostate(zout(end,6:10),zout(end,1:5),beta_f);
ddstatef = ddx(zout(end,1:5),zout(end,6:10),dstatef,dlambdaf);

% Complete constraints matrix
h = [guess(2)-(15+rmoon);
    guess(3)-pi;
    guess(5);
    guess(6)-M*exp((guess(4)-V0)/(Isp*g));
    guess(9)+guess(11)*M/(Isp*g)*exp((guess(4)-V0)/(Isp*g));
    rf-rmoon;
    uf-rmoon*omega;
    vf;
    lambda_phif;
    lambda_mf+1;
    lambda_rf*dstatef(1) + lambda_phif*dstatef(2) + lambda_uf*dstatef(3) + lambda_vf*dstatef(4) + lambda_mf*dstatef(5)];

% Norma
norma = norm(h);

% Count for the iterations
cont = 1;

% Jacobian of the state 
[df,dhz0,dhzf] = compute_jacobian();

% loading the last workspace
load("workspace.mat");
alpha = 0.5;

% Check loop for the algorithm
while norma > 0.0001 
    cont
    if norma<0.05
        alpha = 0.05;
        
    else
        alpha = 0.2;
    end

    % Derivative of 'h' w.r.t. 'tf'
    dhf = [0;
           0;
           0;
           0;
           0;
           dstatef(1);
           dstatef(3);
           dstatef(4);
           dlambdaf(2);
           dlambdaf(5);
           (dlambdaf(1)*dstatef(1)+lambda_rf*ddstatef(1)+dlambdaf(2)*dstatef(2)+lambda_phif*ddstatef(2)+dlambdaf(3)*dstatef(3)+lambda_uf*ddstatef(3)+dlambdaf(4)*dstatef(4)+lambda_vf*ddstatef(4)+dlambdaf(5)*dstatef(5))];
      


    % Computing the update law elements 
    [DHZ0,DHZF,PHI] = linearization(df,dhz0,dhzf,guess(2:end),zout(end,:)',guess(1));
    
    % Equation #24
    ventiquattro = -alpha*inv([dhf DHZ0+DHZF*PHI])*h;
    
    % Updating the initial guesses
    guess = double(guess+ventiquattro);

    
    options = odeset('RelTol',1e-8,'AbsTol',[1e-8, 1e-8,1e-8,1e-8,1e-8,1e-8,1e-8,1e-8,1e-8,1e-8]);

    % Solve the differential augumented state equations 
    [tout,zout] = ode45(@augmented_dynamics,[0 guess(1)],guess(2:end),options);

    % Final augumented state's parameters
    rf = zout(end,1); 
    phif = zout(end,2); 
    uf = zout(end,3); 
    vf = zout(end,4); 
    mf = zout(end,5);
    
    lambda_rf = zout(end,6); 
    lambda_phif = zout(end,7); 
    lambda_uf = zout(end,8); 
    lambda_vf = zout(end,9); 
    lambda_mf = zout(end,10);
    
    % Control input
    beta_f = atan2(-lambda_vf,(-lambda_uf));

    
    % Computing the final value for the differentiated state, costate and 
    % doubled differentiated state 
    dstatef = dx(zout(end,1:5),beta_f);
    dlambdaf = dcostate(zout(end,6:10),zout(end,1:5),beta_f);
    ddstatef = ddx(zout(end,1:5),zout(end,6:10),dstatef,dlambdaf);
    ddcostatef = ddcostate(zout(end,1:5),zout(end,6:10),dstatef,dlambdaf);
 
    % Updating the complete constraints matrix
    h = [guess(2)-(15+rmoon);
         guess(3)-pi;
         guess(5);
         guess(6)-M*exp((guess(4)-V0)/(Isp*g));
         guess(9)+guess(11)*M/(Isp*g)*exp((guess(4)-V0)/(Isp*g));
         rf-rmoon;
         uf-rmoon*omega;
         vf;
         lambda_phif;
         lambda_mf+1;
         lambda_rf*dstatef(1) + lambda_phif*dstatef(2) + lambda_uf*dstatef(3) + lambda_vf*dstatef(4) + lambda_mf*dstatef(5)];

    % Norma
    norma = norm(h)
    cont=cont+1;
end





%[r,phi,u,v,m]%
    