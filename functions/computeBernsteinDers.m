function [B,dB] = computeBernsteinDers(k,u)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%
%  Construct bernstein basis and its derivatives at the paramter points u
%
%  Input:
%    k  - degree, for 1D, k is scalar; for 2D, k=[k1,k2]
%    u  - parameter point,for 1D, u is scalar; for 2D, u=[u,v]
%  Output:
%    B - bernstein basis, [B1,B2,B,...,Bn];
%    dB- bernstein basis derivatives, [dB1du,dB2du,...,dBndu; dB1dv,dB2dv,...,dBndv]
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
    B = computeBernsteinBasis(k-1,u);  % compute basis values of degree k-1
    B1 = [0,B];
    B2 = [B,0];
    dB = k*(B1-B2);
    B = computeBernsteinBasis(k,u); 
elseif 2 == length(k)  % two-dimentional
    [B1,dBdu] = computeBernsteinDers(k(1),u(1));
    [B2,dBdv] = computeBernsteinDers(k(2),u(2));
    dB(1,:) = reshape(dBdu'*B2,1,[]);
    dB(2,:) = reshape(B1'*dBdv,1,[]);
    B  = computeBernsteinBasis(k,u);
end

