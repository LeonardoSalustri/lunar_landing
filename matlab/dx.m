
% Fuction for the state dynamics
function  dstate = dx(x,beta)
    %x = [r, phi, u, v, m] 
    T = 0.440;
    miu = 4902.7779;
    I = 310;
    g = 9.81e-3;
    
    
    dstate(1) = x(4);
    dstate(2) = x(3)/x(1);
    dstate(3) = -x(3)*x(4)/x(1) + T/x(5)*cos(beta);
    dstate(4) = x(3)^2/x(1) - miu/(x(1)^2) + T/x(5)*sin(beta);
    dstate(5) = -T/(I*g);
    dstate = dstate';

end