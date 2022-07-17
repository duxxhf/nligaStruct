function B = computeBernsteinBasis(k,u)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%
%  Construct bernstein basis on [0,1] according to its recurcive
%  properties
%
%  Input:
%    k  - degree, for 1D, k is scalar; for 2D, k=[k1,k2]
%    u  - parameter point,for 1D, u is scalar; for 2D, u=[u,v]
%  Output:
%    B - bernstein basis
%
% please refer to <the nurbs book, 2nd> by piegl et al
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
    if 0 == k
        B(1) = 1.0;
    elseif 0 < k
        B  = zeros(1,k+1);
        B(1) = 1.0;
        u1   = 1.0 - u;
        for j = 2:k+1
            saved = 0.0;
            for i = 1:j-1
                temp = B(i);
                B(i) = saved + u1*temp;
                saved = u*temp;
            end
            B(j) = saved;
        end
    end
elseif 2 == length(k)  % two-dimentional
    B1 = computeBernsteinBasis(k(1),u(1));
    B2 = computeBernsteinBasis(k(2),u(2));
    B = reshape(B1'*B2,1,[]);
end

end


