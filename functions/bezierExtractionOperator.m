function [Ce, nb] = bezierExtractionOperator(U,p,V,q,W,k)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%
%  build local extraction operator for each element
%
%  Input:
%    U,V,W - knot vectors
%    p,q,k - degrees 
%  Output:
%    Ce    - elemental extraction operator
%    nb    - number of elements for knot vector U
%    pleas refer to <Borden,2011,Isogeometric finite element data structures 
%     based on Bezier extraction of NURBS>
%
%  Usage: [Ce, nb] = bezierExtractionOperator(U,p),for 1D
%  Usage: [Ce, nb] = bezierExtractionOperator(U,p,V,q),for 2D
%  Usage: [Ce, nb] = bezierExtractionOperator(U,p,V,q,W,k),for 3D
%
%  ---------------------------------------
%
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%  - 2022
%
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
if nargin == 2  % for 1D
    m   = length(U)-p;
    a   = p+1;
    b   = a+1;
    nb  = 1;
    Ce(:,:,1) = eye(p+1);
    while b < m
        Ce(:,:,nb+1) = eye(p+1);
        ii = b;
        while b < m && U(b+1) == U(b)
            b = b+1;
        end
        mult = b - ii + 1;
        if mult < p    % multiplicity of knots U(b)
            numer = U(b) - U(a);
            for jj = p:-1:mult+1
                alphas(jj-mult) = numer/(U(a+jj)-U(a));
            end
            r = p - mult;
            for jj = 1:r
                saved = r - jj + 1;
                s     = mult + jj;
                for kk = p+1:-1:s+1
                    alpha = alphas(kk-s);
                    Ce(:,kk,nb) = alpha*Ce(:,kk,nb) + (1-alpha)*Ce(:,kk-1,nb);
                end
                if b < m
                    Ce(saved:jj+saved,saved,nb+1) = Ce(p-jj+1:p+1,p+1,nb);
                end
            end
        end  
        nb = nb + 1;
        if b < m
            a = b; 
            b = b+1;
        end 
    end
elseif nargin == 4  % for 2D
    [Ce1, nb1] = bezierExtractionOperator(U,p);
    [Ce2, nb2] = bezierExtractionOperator(V,q);
    Ce = zeros((p+1)*(q+1),(p+1)*(q+1),nb1*nb2);
    nb = nb1*nb2;
    for e2 = 1:nb2
        for e1 = 1:nb1
            eA = nb1*(e2-1) + e1;
            for ii = 1:q+1
                i1 = (p+1)*(ii-1)+1;
                i2 = (p+1)*ii;
                for jj = 1:q+1
                    j1 = (p+1)*(jj-1)+1;
                    j2 = (p+1)*jj;
                    Ce(i1:i2,j1:j2,eA) = Ce2(ii,jj,e2)*Ce1(:,:,e1);
                end
            end
        end
    end
elseif nargin == 6  % for 3D
%    TBA
%     [Ce1, nb1] = bezierExtractionOperator(U,p);
%     [Ce2, nb2] = bezierExtractionOperator(V,q);
%     [Ce3, nb3] = bezierExtractionOperator(W,k);
end

end





















