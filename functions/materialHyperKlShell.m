function [D0,D1,D2,n,m] = materialHyperKlShell(mat, a, A)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%
%  Compute material constitutive matrix for kl shell
%
%  Input:
%    met   - material properties
%    a     - local base vectors defined on midsurface of deformed
%    geometries
%    A     - local base vectors defined on midsurface of undeformed
%    geometries 
%  Output:
%    D0,D1,D2 - constitutive matrices
%
%  ---------------------------------------
%
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%  - 2022
%
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
global gpz wgz
if 40 == mat(1)  % mat=[40,E,nu,h]
    E   = mat(2);
    nu  = mat(3);
    h   = mat(4);
    a3  = cross(a(1,:),a(2,:))/norm(cross(a(1,:),a(2,:)));
    aij = [a(1,:)*a(1,:)'; a(2,:)*a(2,:)'; a(1,:)*a(2,:)']; 
    bij = a(3:5,:)*a3';
    A3  = cross(A(1,:),A(2,:))/norm(cross(A(1,:),A(2,:)));
    Aij = [A(1,:)*A(1,:)';A(2,:)*A(2,:)';A(1,:)*A(2,:)'];
    Bij = A(3:5,:)*A3';
    D0  = zeros(3);   D1 = zeros(3);   D2 = zeros(3); 
    n   = zeros(3,1); m = zeros(3,1);
    for ipt = 1:length(gpz)
        zeta = 0.5*h*gpz(ipt);
        wt   = 0.5*h*wgz(ipt);
        g    = aij - 2*zeta*bij;
        G    = Aij - 2*zeta*Bij;
        tG   = [G(1),G(3); G(3),G(2)];
        cG   = inv(tG);
        C1111 = cG(1,1)*cG(1,1);
        C1122 = nu*cG(1,1)*cG(2,2)+(1-nu)*cG(1,2)*cG(1,2);
        C1112 = cG(1,1)*cG(1,2);
        C2222 = cG(2,2)*cG(2,2);
        C2212 = cG(2,2)*cG(1,2);
        C1212 = 0.5*(1-nu)*cG(1,1)*cG(2,2)+0.5*(1+nu)*cG(1,2)^2;
        Dk    = E/(1-nu*nu)*[C1111,C1122,C1112; C1122,C2222,C2212; C1112,C2212,C1212 ];
        Eij   = 0.5*(g-G);
        Eij(3)= 2*Eij(3);
        S     = Dk*Eij;
        D0    = D0 + Dk*wt;
        D1    = D1 + Dk*wt*zeta;
        D2    = D2 + Dk*wt*zeta^2;      
        n     = n  + S*wt;
        m     = m  + S*zeta*wt;
    end
    
elseif 41 == mat(1) 
    % mat = [index, A10, h], % Neo-Hookean incompressible material
    A10 = mat(2);
    h   = mat(3);
    a3  = cross(a(1,:),a(2,:))/norm(cross(a(1,:),a(2,:)));
    aij = [a(1,:)*a(1,:)'; a(2,:)*a(2,:)'; a(1,:)*a(2,:)']; 
    bij = a(3:5,:)*a3';
    A3  = cross(A(1,:),A(2,:))/norm(cross(A(1,:),A(2,:)));
    Aij = [A(1,:)*A(1,:)';A(2,:)*A(2,:)';A(1,:)*A(2,:)'];
    Bij = A(3:5,:)*A3';
    D0  = zeros(3);   D1 = zeros(3);   D2 = zeros(3); 
    n   = zeros(3,1); m = zeros(3,1);
    for ipt = 1:length(gpz)
        zeta = 0.5*h*gpz(ipt);
        wt   = 0.5*h*wgz(ipt);
        g    = aij - 2*zeta*bij;
        tg   = [g(1),g(3); g(3),g(2)];
        G    = Aij - 2*zeta*Bij;
        tG   = [G(1),G(3); G(3),G(2)];
        cg   = inv(tg);
        cG   = inv(tG);
        J0   = sqrt(det(tg)/det(tG));
        tS   = 2*A10*( cG - cg/(J0^2) );
        S    = [tS(1,1);tS(2,2);tS(1,2)];
        C1111 = 2*A10*J0^(-2)*( 2*cg(1,1)*cg(1,1)+cg(1,1)*cg(1,1)+cg(1,1)*cg(1,1));
        C1122 = 2*A10*J0^(-2)*( 2*cg(1,1)*cg(2,2)+cg(1,2)*cg(1,2)+cg(1,2)*cg(1,2));
        C1112 = 2*A10*J0^(-2)*( 2*cg(1,1)*cg(1,2)+cg(1,1)*cg(1,2)+cg(1,2)*cg(1,1));
        C2222 = 2*A10*J0^(-2)*( 2*cg(2,2)*cg(2,2)+cg(2,2)*cg(2,2)+cg(2,2)*cg(2,2));
        C2212 = 2*A10*J0^(-2)*( 2*cg(2,2)*cg(1,2)+cg(2,1)*cg(2,2)+cg(2,2)*cg(2,1));
        C1212 = 2*A10*J0^(-2)*( 2*cg(1,2)*cg(1,2)+cg(1,1)*cg(2,2)+cg(1,2)*cg(2,1));
        Dk    = [C1111,C1122,C1112; C1122,C2222,C2212; C1112,C2212,C1212 ];
        D0    = D0 + Dk*wt;
        D1    = D1 + Dk*wt*zeta;
        D2    = D2 + Dk*wt*zeta^2;      
        n     = n  + S*wt;
        m     = m  + S*zeta*wt;
    end
elseif 42 == mat(1)
    % mat = [index, mu, lambda, h], % compressible Neo-Hookean material
    type  = 2; % here we consider two types of compressible Neo-Hookean material
    mu    = mat(2);                 % Lame's second parameter
    lambda = mat(3);                % Lame's first parameter
    K     = lambda + 2*mu/3;        % Bulk modulus 
    h     = mat(4);                 % thickness
    a3    = cross(a(1,:),a(2,:))/norm(cross(a(1,:),a(2,:)));    
    aij   = [a(1,:)*a(1,:)'; a(2,:)*a(2,:)'; a(1,:)*a(2,:)']; 
    bij   = a(3:5,:)*a3';
    A3    = cross(A(1,:),A(2,:))/norm(cross(A(1,:),A(2,:)));
    Aij   = [A(1,:)*A(1,:)';A(2,:)*A(2,:)';A(1,:)*A(2,:)'];
    Bij   = A(3:5,:)*A3';
    D0    = zeros(3);   D1 = zeros(3);   D2 = zeros(3); 
    n     = zeros(3,1); m = zeros(3,1);
    for ipt = 1:length(gpz)
        zeta = 0.5*h*gpz(ipt);
        wt   = 0.5*h*wgz(ipt);
        g    = aij - 2*zeta*bij;
        tg   = [g(1),g(3); g(3),g(2)];
        G    = Aij - 2*zeta*Bij;
        tG   = [G(1),G(3); G(3),G(2)];
        cg   = inv(tg);
        cG   = inv(tG);
        J0   = sqrt(det(tg)/det(tG));
        S33  = 1E3;
        C33  = 1; tol  = 1e-3;
        while abs(S33) > tol
            J     = J0*sqrt(C33);
            if type == 1
                trC   = tg(1,1)*cG(1,1) + tg(1,2)*cG(1,2) + tg(2,1)*cG(2,1) + tg(2,2)*cG(2,2) + C33;
                S33   = mu*J^(-2/3)*(1 - 1/3*trC/C33) + K/2*(J^2-1)/C33;
                C3333 = 1/9*mu*J^(-2/3)*( trC*8/C33^2-12/C33 ) ...
                    + K*(J^2/C33^2-1/2*(J^2-1)*2/C33^2);
            elseif type == 2
                S33   = mu - (mu-lambda/2*(J^2-1))/C33;
                C3333 = (lambda+2*mu)/C33^2;
            end
            dC33  = -2*S33/C3333;
            C33   = C33 + dC33;         
        end
        Cb   = blkdiag(cg,1/C33);
        cG   = blkdiag(cG,1);
        J     = J0*sqrt(C33);
        trC   = tg(1,1)*cG(1,1) + tg(1,2)*cG(1,2) + tg(2,1)*cG(2,1) + tg(2,2)*cG(2,2) + C33;
        [Sij,Cijkl] = compressibleNeo(cG,Cb,trC,J,mu,lambda,type);
        S    = [Sij(1,1);Sij(2,2);Sij(1,2)];
        C1111  = Cijkl(1,1,1,1) - Cijkl(1,1,3,3)*Cijkl(3,3,1,1)/Cijkl(3,3,3,3);
        C1122  = Cijkl(1,1,2,2) - Cijkl(1,1,3,3)*Cijkl(3,3,2,2)/Cijkl(3,3,3,3);
        C1112  = Cijkl(1,1,1,2) - Cijkl(1,1,3,3)*Cijkl(3,3,1,2)/Cijkl(3,3,3,3);
        C2222  = Cijkl(2,2,2,2) - Cijkl(2,2,3,3)*Cijkl(3,3,2,2)/Cijkl(3,3,3,3);
        C2212  = Cijkl(2,2,1,2) - Cijkl(2,2,3,3)*Cijkl(3,3,1,2)/Cijkl(3,3,3,3);
        C1212  = Cijkl(1,2,1,2) - Cijkl(1,2,3,3)*Cijkl(3,3,1,2)/Cijkl(3,3,3,3);
        Dk    = [C1111,C1122,C1112; C1122,C2222,C2212; C1112,C2212,C1212 ];
        D0 = D0 + Dk*wt;
        D1 = D1 + Dk*wt*zeta;
        D2 = D2 + Dk*wt*zeta^2;      
        n  = n  + S*wt;
        m  = m  + S*zeta*wt;
    end
end

end