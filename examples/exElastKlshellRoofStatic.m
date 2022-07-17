%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%
%  An example for structural analysis of the benchmark problem 
%  - Scordelis¨CLo roof using KL shell
%
%  ---------------------------------------
%
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%  - 2022
%
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

clc; 
clear;

%****************************************************************%
%             Import and plot geometrical model                  %
%****************************************************************%
fprintf("Import and plot geometrical model \n");
meshType = 'T-SPLINE';  % or 'NURBS'
% meshType = 'NURBS';  % or 'NURBS'
if strcmp(meshType, 'T-SPLINE')
    R = 25; L = 50;
    filename = 'geoQuarterRoofR5.iga';
    geo  = readIgaFile(filename);
    mesh = buildIgaMesh( geo );
elseif strcmp(meshType, 'NURBS')
    R = 25; L = 50; theta = 40; p = 4; q = p; elemX = 10; elemY = elemX; 
    geo  = geoQuarterRoofShellNurbs(L,R,theta,p,q,elemX,elemY);
    mesh = buildIgaMesh( geo );
end

figure(1)
plotIgaMesh(mesh,1);
view(37,37);
%****************************************************************%
%               Define material properties and dof               %
%****************************************************************%
fprintf("Define material properties and dof \n");
tic
% material properties
nu   = 0.0;             % Poisson ratio
E    = 4.32e8;          % Youngs modulus 
g    = 90;            % Gravity loading
h    = 0.25;            % thickness   
dof  = 3;               % degree of freedom
toc

%****************************************************************%
%                   Enforce boundary conditions                  %
%****************************************************************%
fprintf("Enforce boundary conditions \n");
[dbc,sym_dbc,nodeA] = setBcRoof(mesh);
scatdbc = [];
scattbc = [];
if ~isempty(dbc)
    scatdbc = dof * (dbc(:,1)-1) + dbc(:,2);   % scatter dbc
end
toc

%****************************************************************%
%           Initialize and assemble stiffness matrix             %
%****************************************************************%
fprintf("Initialize and assemble stiffness matrix \n");
tic
nDofs = dof * mesh.nCpts;    % total dofs
K = sparse(nDofs,nDofs);     % stiffness matrix 
F = zeros(nDofs,1);          % external force matrix

% use gaussian integration rule
gp_x = max(mesh.elDegree(:,1))+1;           % number of integration points in x-direction
gp_y = max(mesh.elDegree(:,2))+1;           % number of integration points in y-direction
[gp, wgt] = gaussQuadrature(gp_x, gp_y);   % calculate integration points and its weights
elDoma = [0,1,0,1];
for e = 1:mesh.nElems                 % loop over elements
    sctr   = mesh.elNodeCnt{e,:};       % element control points index
    elCpts = mesh.coords(sctr,1:3);     % coordinates of element control points
    nn     = numel(sctr);                 % number of control points for each element
    nnElem = nn*dof;                  % dof for each element
    sctrB = zeros(1, nnElem); 
    for i = 1:dof
        sctrB(i:dof:nnElem) = dof*(sctr-1) + i;  % displacement in i-th direction
    end
    pu = mesh.elDegree(e,1);
    pv = mesh.elDegree(e,2);
    Ce = mesh.elExtOpe{e,1};
    we = mesh.coords(sctr,4); % Tspline control points' weights
    Ke = zeros(nnElem,nnElem);
    for ipt = 1:size(gp,1)   % loop over integration points
        pt   = gp(ipt,1:2);      % reference parametric coordinates for each integration point
        wt   = wgt(ipt);       % weigths for each integration point
        gPts = parameterGaussMapping( elDoma, pt );   % gauss integration mapping
        j1   = jacobianGaussMapping( elDoma );         % jacobian value for gauss mapping   
        [~,dR,dR2] = computeTsplineBasis2ndDers([pu,pv],gPts,Ce,we);
        Aab = [dR;dR2]*elCpts;
        A1  = Aab(1,:)'; A2  = Aab(2,:)'; 
        A11 = Aab(3,:)'; A22 = Aab(4,:)'; A12 = Aab(5,:)';
        j2  = norm(cross(A1,A2));
        A3  = cross(A1,A2)/norm(cross(A1,A2));
        Gab = Aab(1:2,:)*Aab(1:2,:)';
        Gab = inv(Gab);
        Dkl = [Gab(1,1)^2,  nu*Gab(1,1)*Gab(2,2)+(1-nu)*Gab(1,2)^2,  Gab(1,1)*Gab(1,2);
               nu*Gab(1,1)*Gab(2,2)+(1-nu)*Gab(1,2)^2,   Gab(2,2)^2, Gab(2,2)*Gab(1,2);
               Gab(1,1)*Gab(1,2), Gab(2,2)*Gab(1,2), 0.5*(1-nu)*Gab(1,1)*Gab(2,2)+0.5*(1+nu)*Gab(1,2)^2];
        Dkl = Dkl*E/(1-nu^2);
        D0  = h*Dkl; D2 = h^3*Dkl/12;
        Rm  = [transR(dR(1,:)); transR(dR(2,:))];
        Rb  = [Rm; transR(dR2(1,:)); transR(dR2(2,:)); transR(dR2(3,:))];
        Ee  = [A1, zeros(3,1), A2; zeros(3,1), A2, A1];
        Ia  = eye(3)-A3*A3'; Ta = [transV(A2), -transV(A1)];
        Ea  = (Ia*Ta)./norm(cross(A1,A2));
        Ek  = -[Ea'*A11,Ea'*A22,2*Ea'*A12; A3,zeros(3,2); zeros(3,1),A3,zeros(3,1); zeros(3,2),2*A3];
        fac = j1 *j2 * wt;
        Ke  = Ke + (Rm'*Ee*D0*Ee'*Rm + Rb'*Ek*D2*Ek'*Rb)*fac;
    end   
    K(sctrB,sctrB) = K(sctrB,sctrB) +  Ke;
end
toc;

%****************************************************************%
%                   Enforce gravity loading                      %
%****************************************************************%
fprintf("Enforce gravity loading \n");
tic
for e = 1:mesh.nElems                 % loop over elements
    sctr   = mesh.elNodeCnt{e,:};     % element control points index
    elCpts = mesh.coords(sctr,:);     % coordinates of element control points
    strf   = sctr*3;
    pu     = mesh.elDegree(e,1);
    pv     = mesh.elDegree(e,2);
    Ce     = mesh.elExtOpe{e,1};
    we     = mesh.coords(sctr,4); % Tspline control points' weights
    for ipt = 1:size(gp,1)        % loop over integration points
        pt     = gp(ipt,:);      % reference parametric coordinates for each integration point
        wt     = wgt(ipt);       % weigths for each integration point
        gPts   = parameterGaussMapping( elDoma, pt );       % gauss integration mapping
        j1     = jacobianGaussMapping( elDoma );            % jacobian value for gauss mapping   
        [R,dR] = computeTsplineBasisDers([pu,pv],gPts,Ce,we);
        jmatrix  = dR*elCpts(:,1:3); 
        j2       = norm(cross(jmatrix(1,:),jmatrix(2,:)));
        fac      = j1 *j2 * wt;        
        F(strf)  = F(strf) -  R' * g * fac;
    end
end
toc   

%****************************************************************%
%                   Solve stiffness equation                     %
%****************************************************************%
sctrs = [];
sctrs = [sctrs;3 * (sym_dbc-1)+1];
sctrs = [sctrs;3 * (sym_dbc-1)+2];
sctrs = [sctrs;3 * (sym_dbc-1)+3];
alpha  = max(diag(K))*1e5;
pStiff = alpha*[1 -1;-1 1];
for i = 1:size(sctrs,1)  
    sctri = sctrs(i,:);
    K(sctri,sctri) = K(sctri,sctri) + pStiff;
end
fprintf("Solve stiffness equation \n");
tic
ndbc = size(dbc,1);   % number of displacement constrained nodes
K(scatdbc,:) = zeros(ndbc, nDofs);
K(scatdbc,scatdbc) = eye(ndbc);
F(scatdbc,:) = dbc(:,3);  
u = K\F;
toc

%****************************************************************%
%                       Post-processing                          %
%****************************************************************%
fprintf("Post-processing \n");
tic
% scale factor for better visualization
u = reshape(u,3,[])';
vmesh = getVisualMesh(mesh,u);
msize = getMeshSize(mesh.coords);
len = norm([msize(1), msize(3), msize(5)] - [msize(2), msize(4), msize(6)]);
% factor = len/max(abs(u(:,3)))/10;  
factor = 10;
vmesh.node = vmesh.node + factor*vmesh.u;
vMesh{1,1} = vmesh;
vMesh{1,2} = vmesh;
vMesh{1,2}.node(:,2) = L - vMesh{1,2}.node(:,2);
vMesh{1,3} = vmesh;
vMesh{1,3}.node(:,1) = - vMesh{1,3}.node(:,1);
vMesh{1,4} = vMesh{1,2}; 
vMesh{1,4}.node(:,1) = - vMesh{1,4}.node(:,1); 
figure(2)
for i = 1:4
    hold on;
    p = patch('Faces',vMesh{1,i}.element, 'Vertices', vMesh{1,i}.node);         % plot the surface
    cdata = vMesh{1,i}.u(:,3);
    set(p,'FaceColor','interp','FaceVertexCData',cdata,'EdgeColor','none');      
    for j = 1:mesh.nElems
        hold on;
        plot3(vMesh{1,i}.node(vMesh{1,i}.linmesh{j,1},1),...
              vMesh{1,i}.node(vMesh{1,i}.linmesh{j,1},2),...
              vMesh{1,i}.node(vMesh{1,i}.linmesh{j,1},3),'k-');
    end
end

view(37,37);
axis equal;
C = colorbar;
caxis manual
colorMin = C.Limits(1);
colorMax = C.Limits(2);
caxis([colorMin colorMax]);
set(C,'Ticks',linspace(colorMin,colorMax,7));
mycolor = abaqusColorMap(24);
colormap(mycolor);
toc;

fprintf('The maximum vertical displacement is %12.6e.\n',-u(nodeA,3));

