function dB2 = computeBernstein2ndDers(k,u)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%
%  Compute the second derivatives of non-zero bernstein basis at the paramter points u
%
%  Input:
%    k  - degree, for 1D, k is scalar; for 2D, k=[k1,k2]
%    u  - parameter point,for 1D, u is scalar; for 2D, u=[u,v]
%  Output:
%    dB- bernstein basis 2nd derivatives,dB2= [dB2du; dB2dv; dB2duv]
%
%  ---------------------------------------
%
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%  - 29-NOV-2021
%
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

if 1 == length(k)   % one-dimentional
    B  = computeBernsteinBasis(k-2,u);  % compute basis values of degree k-2
    B1  = [0,0,B];
    B2  = [0,B,0];
    B3  = [B,0,0];
    dB2 = k*(k-1)*(B1-2*B2+B3);
elseif 2 == length(k)  % two-dimentional
    dBdu2    = computeBernstein2ndDers(k(1),u(1));
    dBdv2    = computeBernstein2ndDers(k(2),u(2));
    [Bu,dBdu] = computeBernsteinDers(k(1),u(1));
    [Bv,dBdv] = computeBernsteinDers(k(2),u(2));
    dBdu2    = reshape(dBdu2'*Bv,1,[]);
    dBdv2    = reshape(Bu'*dBdv2,1,[]);
    dBduv2   = reshape(dBdu'*dBdv,1,[]);
    dB2      = [dBdu2; dBdv2; dBduv2];
end

