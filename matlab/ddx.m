
% Function for the derivative of the state dynamics
function [ddstate] = ddx(state,costate,dstate,dcostate)
    T = 0.440;
    %miu = 2.4691e+10;
    miu = 4.899803828918445e+03;
    I = 310;
    g = 9.81e-3;
    ddstate(1) = dstate(4);
    ddstate(2) = (dstate(3)*state(1)-state(3)*dstate(1))/state(1)^2;
    ddstate(3) = -(dstate(3)*state(4)*state(1)+state(3)*dstate(4)*state(1)-dstate(1)*state(3)*state(4))/(state(1)^2)+(-T*sin(atan2d(-costate(4),(-costate(3))))*(-costate(3)*dcostate(4)+costate(4)*dcostate(3))/(costate(3)^2+costate(4)^2)*state(5)-T*dstate(5)*cos(atan2d(-costate(4),(-costate(3)))))/state(5)^2;
    ddstate(4) = (2*state(3)*dstate(3)*state(1)-state(3)^2*dstate(1))/(state(1)^2)+2*miu*dstate(1)/(state(1)^3)+(T*cos(atan2d(-costate(4),(-costate(3))))*(-costate(3)*dcostate(4)+costate(4)*dcostate(3))/(costate(3)^2+costate(4)^2)*state(5)-T*sin(atan2d(-costate(4),(-costate(3))))*dstate(5))/state(5)^2;
    ddstate(5) = 0;
    ddstate = ddstate';
end



