function Hm = computeHmMatrix(m,Ea,Ta,Ia,a3,a11,a22,a12)
% Compute Hm matrix 

I11 = [zeros(3), transV(Ia*a11); -transV(Ia*a11),zeros(3)];
I22 = [zeros(3), transV(Ia*a22); -transV(Ia*a22),zeros(3)];
I12 = [zeros(3), transV(Ia*a12); -transV(Ia*a12),zeros(3)];
H11 = Ea'*a11*a3'*Ta + Ta'*a3*a11'*Ea + a11'*a3*Ta'*Ea - I11;
H22 = Ea'*a22*a3'*Ta + Ta'*a3*a22'*Ea + a22'*a3*Ta'*Ea - I22;
H12 = Ea'*a12*a3'*Ta + Ta'*a3*a12'*Ea + a12'*a3*Ta'*Ea - I12;
Hm = m(1)*H11 + m(2)*H22 + 2*m(3)*H12;
end

