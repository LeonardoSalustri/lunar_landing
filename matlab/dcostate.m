
% Funcrion for the costate dynamics
function dlambda = dcostate(lambda,x,beta)
    T = 0.440;
    miu = 4902.7779;

    
    
    dlambda(1) = lambda(2)*x(3)/(x(1)^2) - lambda(3)*x(3)*x(4)/(x(1)^2) + lambda(4)*(x(3)^2/(x(1)^2) - 2*miu/(x(1)^3));
    dlambda(2) = 0;
    dlambda(3) = -lambda(2)/x(1) + lambda(3)*x(4)/x(1) - lambda(4)*2*x(3)/x(1);
    dlambda(4) = -lambda(1) + lambda(3)*x(3)/x(1);
    dlambda(5) = T/x(5)^2*(lambda(3)*cos(beta)+lambda(4)*sin(beta));
    
    dlambda=dlambda';
end