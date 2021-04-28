% Function for computing the jacobian of the augumented state
function [df,dhz0,dhzf] = compute_jacobian()
    
    rmoon = 1737.1;
    Isp = 310;
    M = 300;
    V0 = 1.69196926;
    g = 9.81e-3;
    omega = 2.7e-6;
    
    
    syms z1_0 z2_0 z3_0 z4_0 z5_0 z6_0 z7_0 z8_0 z9_0 z10_0 z1_f z2_f z3_f z4_f z5_f z6_f z7_f z8_f z9_f z10_f real;
    z0 = [z1_0; z2_0; z3_0; z4_0; z5_0; z6_0; z7_0; z8_0; z9_0; z10_0];
    zf = [z1_f; z2_f; z3_f; z4_f; z5_f; z6_f; z7_f; z8_f; z9_f; z10_f];
 
    dstate = dx([z1_f;z2_f;z3_f;z4_f;z5_f],atan2(-z9_f,(-z8_f)));
    dlambda = dcostate([z6_f;z7_f;z8_f;z9_f;z10_f],[z1_f;z2_f;z3_f;z4_f;z5_f],atan2(-z9_f,(-z8_f)));
    f = [dstate;dlambda];
    df = jacobian(f,[z1_f z2_f z3_f z4_f z5_f z6_f z7_f z8_f z9_f z10_f]);
    
    h = [z1_0-(15+rmoon);
        z2_0-pi;
        z4_0;
        z5_0-M*exp((z3_0-V0)/(Isp*g));
        z8_0+z10_0*M/(Isp*g)*exp((z3_0-V0)/(Isp*g));
        z1_f-rmoon;
        z3_f-rmoon*omega;
        z4_f;
        z7_f;
        z10_f+1;
        z6_f*dstate(1)+z7_f*dstate(2)+z8_f*dstate(3)+z9_f*dstate(4)+z10_f*dstate(5)];

    dhz0 = jacobian(h,z0);
    dhzf = jacobian(h,zf);

end

