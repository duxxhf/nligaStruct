function [ Kt, P ] = globalStiffnessHyperKlShell( mesh, mat, u )
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%
%  Assemble global stiffness matrix for nonlinear kl shell
%
%  Input:
%    mesh  - mesh structure
%    mat   - material properties
%    u     - displacement 
%  Output:
%    k - global stiffness
%    P - internal force
%
%  ---------------------------------------
%
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%  - 18-JAN-2022
%
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

global gp wg                            % Gauss points and weights for integration
elDoma = [0,1,0,1];                     % Parametric domain for Bezier elements
dof = 3;  ndofs = dof * mesh.nCpts;     % Total number of nodal dofs and global dofs
Kt  = sparse(ndofs, ndofs);             % Initialize stiffness matrix
P   = zeros(ndofs, 1);                  % Initialize residual vector
for e = 1:mesh.nElems                   % Loop over elements
    sctr   = mesh.elNodeCnt{e,:};       % Global indices of elemental control points
    elCpts = mesh.coords(sctr,1:3);     % Coordinates of elemental control points
    nn     = numel(sctr);               % Number of control points for the element
    nnElem = nn*dof;                    % Dofs for the element
    sctrB = zeros(1, nnElem);           % Scatter nodal dofs
    for i = 1:dof
        sctrB(i:dof:nnElem) = dof*(sctr-1) + i;
    end
    elDisp = reshape(u(sctrB),dof,nn)'; % Elemental displacements
    pu = mesh.elDegree(e,1);            % Degree for the Bezier element
    pv = mesh.elDegree(e,2);            % Degree for the Bezier element
    Ce = mesh.elExtOpe{e,1};            % Extraction operator for the Bezier element
    we = mesh.coords(sctr,4);           % Weights for T-spline elemental control points
    for ipt = 1:size(gp,1)              % Loop over integration points
        pt = gp(ipt,:);                 % Coordinates of the ipt-th integration point
        wt = wg(ipt);                   % Weight of the ipt-th integration point
        % Gauss integration mapping
        gauPts = parameterGaussMapping(elDoma,pt);
        % Jacobian value for gauss mapping, = 0.25
        j1     = jacobianGaussMapping( elDoma );
        % The first and second derivatives of T-spline basis functions
        [~,dR,dR2] = computeTsplineBasis2ndDers([pu,pv],gauPts,Ce,we);
        newCpts    = elCpts + elDisp;   % Deformed control points
        aab = [dR;dR2]*newCpts;         % Derivatives of the deformed configuration
        Aab = [dR;dR2]*elCpts;          % Derivatives of the initial configuration
        % Tangent derivatives and normal vector
        a1  = aab(1,:)'; a2  = aab(2,:)'; a3  = cross(a1,a2)/norm(cross(a1,a2));
        % Second derivatives of the deformed configuration
        a11 = aab(3,:)'; a22 = aab(4,:)'; a12 = aab(5,:)';
        % Tangent derivatives of the initial configurtion
        A1  = Aab(1,:)'; A2  = Aab(2,:)';
        % Jacobian value for parametric mapping
        j2  = norm(cross(A1,A2));
        % Compute material constitutive matrices, force and moment resultants
        [D0,D1,D2,n,m] = materialHyperKlShell(mat,aab,Aab);
        Rm  = [transR(dR(1,:));transR(dR(2,:))]; Ia = eye(3)-a3*a3';        % Eq. (49)
        Ta  = [transV(a2),-transV(a1)]; Ea = (Ia*Ta)./norm(cross(a1,a2));   % Eq. (49)
        Rb  = [Rm; transR(dR2(1,:)); transR(dR2(2,:)); transR(dR2(3,:))];   % Eq. (52)
        Ee  = [a1, zeros(3,1), a2; zeros(3,1), a2, a1];                     % Eq. (52)
        Ek  = -[Ea'*a11,Ea'*a22,2*Ea'*a12; a3,zeros(3,2);                   % Eq. (52)
                zeros(3,1),a3,zeros(3,1); zeros(3,2),2*a3];
        En  = [n(1)*eye(3), n(3)*eye(3); n(3)*eye(3),n(2)*eye(3)];          % Eq. (54)
        Em  = -[zeros(6,6), m(1)*Ea', m(2)*Ea', 2*m(3)*Ea';                 % Eq. (56)
                m(1)*Ea, zeros(3,9); m(2)*Ea, zeros(3,9); 2*m(3)*Ea, zeros(3,9)];
        Hm  = computeHmMatrix(m,Ea,Ta,Ia,a3,a11,a22,a12)/norm(cross(a1,a2));% Eq. (56)
        BL  = [Ee'*Rm; Ek'*Rb];                                             % Eq. (61)
        D   = [D0,D1;D1,D2];                                                % Eq. (61)
        fac = j1*j2*wt;
        % Stiffness matrix and residual vector
        Kt(sctrB,sctrB) = Kt(sctrB,sctrB) + ...                             % Eq. (59)
             (BL'*D*BL + Rm'*En*Rm + Rm'*Hm*Rm + Rb'*Em*Rb)*fac;
        P(sctrB)  = P(sctrB) + (Rm'*Ee*n + Rb'*Ek*m)*fac;                   % Eq. (60)
    end
end

end









