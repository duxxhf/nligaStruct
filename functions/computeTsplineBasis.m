function T = computeTsplineBasis(k,u,Ce,we)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%
%  Construct tspline basis functions and its 1st derivatives from bernstein basis 
%
%  Input:
%    k  - degree, for 1D, k is scalar; for 2D, k=[k1,k2]
%    u  - parameter point,for 1D, u is scalar; for 2D, u=[u,v]
%    Ce - elemental bezier extraction operator
%    we - weights vector for tspline comtrol points
%  Output:
%    T - tspline basis
%
%  ---------------------------------------
%
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%  - 29-NOV-2021
%
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
We = diag(we);
B = computeBernsteinBasis(k,u);
B  = B';
W  = we'*Ce*B;  
T  = We*Ce*B/W;
T  = T';
end

