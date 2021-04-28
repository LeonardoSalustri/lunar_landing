% Compute the Hamiltonian after solving the shooting method 

function [H] = compute_hamiltonian(tout,zout)
    load("workspace.mat");
    H = zeros(size(tout));
    for i=1:1:length(tout)
        beta = atan2(-zout(i,9),-zout(i,8));
        dstate = dx(zout(i,:)',beta);
        h = dstate(1).*zout(i,6)+dstate(2)*zout(i,7)+dstate(3)*zout(i,8)+dstate(4)*zout(i,9)+dstate(5)*zout(i,10);
        H(i) = h;
    end
    
    % Plotting
    figure('DefaultAxesFontSize',13);
    
    plot(tout,H,'linewidth',4);
    grid;
    
    title("Hamiltonian function",'Interpreter','latex','fontsize',24);
    xlabel("Time [s]",'fontsize',15);
    ylabel("H",'fontsize',15);

end