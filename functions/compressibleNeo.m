function [Sij,Cijkl] = compressibleNeo(G,iC,trC,J,mu,lambda,type)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%
%  Compute stress and constitutive tensor for hyperelastic material
%
%  Input:
%    G     - G^{ij} contravariant metric cofficients
%    iC    - inverse of deformaton tensor C
%    iC    - trace of deformaton tensor C
%    J     - Jacobian determinant
%    mu    - Lame's second parameter
%   lambda - Lame's first parameter
%  Output:
%    Sij   - the second Piola-Kirchhoff stress tensor
%    Cijkl - material tensor
%
%  ---------------------------------------
%
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%  - 2022
%
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
K = lambda + 2*mu/3;
Sij = zeros(3,3);
Cijkl = zeros(3,3,3,3);
if type == 1
    for j=1:3
        for i = 1:3
            Sij(i,j) = mu*J^(-2/3)*(G(i,j) - 1/3*trC*iC(i,j)) + K/2*(J^2-1)*iC(i,j);
        end
    end
    for j=1:3
        for i = 1:3
            for k=1:3
                for l=1:3
                    Cijkl(i,j,k,l) = 1/9*mu*J^(-2/3)*( trC*(2*iC(i,j)*iC(k,l) ...
                        + 3*iC(i,k)*iC(j,l) + 3*iC(i,l)*iC(j,k)) ...
                        - 6*(G(i,j)*iC(k,l) + iC(i,j)*G(k,l))) + K*(J^2*iC(i,j)*iC(k,l)...
                        - 0.5*(J^2-1)*( iC(i,k)*iC(j,l) + iC(i,l)*iC(j,k)));
                end
            end
        end
    end
elseif type == 2
    for j=1:3
        for i = 1:3
            Sij(i,j) = mu*G(i,j) - (mu-lambda/2*(J^2-1))*iC(i,j);
        end
    end
    for j=1:3
        for i = 1:3
            for k=1:3
                for l=1:3
                    Cijkl(i,j,k,l) = lambda*J^2*iC(i,j)*iC(k,l) ...
                        + (mu-lambda/2*(J^2-1))*( iC(i,k)*iC(l,j) + iC(i,l)*iC(k,j)  );
                end
            end
        end
    end
end


end

