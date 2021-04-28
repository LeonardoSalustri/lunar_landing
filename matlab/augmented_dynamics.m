% Function for the augumented state dynamincs
function dz = augmented_dynamics(t,z)
    
    % Control input 
    beta = atan2(-z(9),(-z(8)));
    
    dstate = dx(z(1:5),beta);
    dlambda = dcostate(z(6:10),z(1:5),beta);
    
    dz = [dstate;dlambda];
end

