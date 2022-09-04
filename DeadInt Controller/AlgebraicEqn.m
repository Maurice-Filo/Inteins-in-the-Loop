function [z_5, z_6, z_7, z_8] = AlgebraicEqn(z_p, z_s, a, d, IG)
    z = lsqnonlin(@(z) MyFun(z, z_p, z_s, a, d), IG, zeros(4,1));
    z_5 = z(1); z_6 = z(2); z_7 = z(3); z_8 = z(4);
end

function Value = MyFun(z, z_p, z_s, a, d)
    a_1 = a(1); a_2 = a(2); a_3 = a(3);
    a_4 = a(4); a_5 = a(5); a_6 = a(6);
    d_1 = d(1); d_2 = d(2); d_3 = d(3);
    d_4 = d(4); d_5 = d(5); d_6 = d(6);
    z_5 = z(1); z_6 = z(2); z_7 = z(3); z_8 = z(4);
    Value = [d_2*z_6 - d_1*z_5 + a_1*(2*z_5 + 2*z_6 + 2*z_7 + z_8 - z_s)^2 + a_2*z_5*(z_6 + 2*z_7 + z_8 - z_p); ...
             d_3*z_7 - d_2*z_6 - d_6*z_6 - a_2*z_5*(z_6 + 2*z_7 + z_8 - z_p) + a_3*z_6*(z_6 + 2*z_7 + z_8 - z_p) - a_6*z_8*(2*z_5 + 2*z_6 + 2*z_7 + z_8 - z_s); ...
             a_5*z_8^2 - d_5*z_7 - a_3*z_6^2 - d_3*z_7 - 2*a_3*z_6*z_7 - a_3*z_6*z_8 + a_3*z_6*z_p; ...
             2*d_5*z_7 - d_4*z_8 + d_6*z_6 - 2*a_5*z_8^2 + a_4*(z_6 + 2*z_7 + z_8 - z_p)*(2*z_5 + 2*z_6 + 2*z_7 + z_8 - z_s) + a_6*z_8*(2*z_5 + 2*z_6 + 2*z_7 + z_8 - z_s)];
    
end
