function [DHZ0, DHZF, PHI] = linearization(df,dhz0,dhzf,zin,zfin,tf)
    
    syms z1_0 z2_0 z3_0 z4_0 z5_0 z6_0 z7_0 z8_0 z9_0 z10_0 z1_f z2_f z3_f z4_f z5_f z6_f z7_f z8_f z9_f z10_f real;
    z0 = [z1_0; z2_0; z3_0; z4_0; z5_0; z6_0; z7_0; z8_0; z9_0; z10_0];
    zf = [z1_f; z2_f; z3_f; z4_f; z5_f; z6_f; z7_f; z8_f; z9_f; z10_f];
    
    df = double(subs(df,zf,zin));
    df = expm(df*tf);
    PHI = double(vpa(df,4));
   
    dhz0 = subs(dhz0,z0,zin);
    DHZ0 = double(vpa(dhz0,4));
    dhzf = subs(dhzf,zf,zfin);
    DHZF = double(vpa(dhzf,4));
    
end
