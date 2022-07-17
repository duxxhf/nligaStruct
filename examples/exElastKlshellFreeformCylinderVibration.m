%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%
%  An example for structural vibration using KL shell
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
tic
L    = 10;
geo  = geoCurvedCylinder(L);
mesh = buildIgaMesh( geo );

% plot scaled model
figure(1)
plotIgaMesh(mesh,1);
axis([0,1,0,L,0,1]);
view(10,5);
toc

%****************************************************************%
%               Define material properties and dof               %
%****************************************************************%
fprintf("Define material properties and dof \n");
% material properties
nu   = 0.3;             % Poisson ratio
E    = 210e9;          % Youngs modulus 
h    = 0.004;            % thickness   
dof  = 3;               % degree of freedom
rho  = 7800;            % density 
nModes = 10;            % number of mode shapes


%****************************************************************%
%                   Enforce boundary conditions                  %
%****************************************************************%
bc = 1;  % =1 both ends clamped; = 2 one end clamped, another free
% = 3 one end clamped, another simply supported uy = 0
fprintf("Enforce boundary conditions \n");
tic
% find four boundaries
for i = 1:length(mesh.nodeSets)
    if strcmp(mesh.nodeSets{i,1}.name,'bottom')
        bottomNodes = mesh.nodeSets{i,1}.gloInx;
    elseif strcmp(mesh.nodeSets{i,1}.name,'nextBottom')
        nextBottomNodes = mesh.nodeSets{i,1}.gloInx;
    elseif strcmp(mesh.nodeSets{i,1}.name,'top')
        topNodes = mesh.nodeSets{i,1}.gloInx;
    elseif strcmp(mesh.nodeSets{i,1}.name,'nextTop')
        nextTopNodes = mesh.nodeSets{i,1}.gloInx;
    end
end

% plot constrained control points
hold on;
plot3(mesh.coords(bottomNodes,1),mesh.coords(bottomNodes,2),mesh.coords(bottomNodes,3),'co');
plot3(mesh.coords(nextBottomNodes,1),mesh.coords(nextBottomNodes,2),mesh.coords(nextBottomNodes,3),'c*');
plot3(mesh.coords(topNodes,1),mesh.coords(topNodes,2),mesh.coords(topNodes,3),'ro');
plot3(mesh.coords(nextTopNodes,1),mesh.coords(nextTopNodes,2),mesh.coords(nextTopNodes,3),'r*');

if bc == 1  % both ends clamped
    fixedNodes  = unique([bottomNodes;nextBottomNodes;topNodes;nextTopNodes]);
    fixed_dof   = zeros(length(fixedNodes)*3,1);
    fixed_dof(1:3:end) = (fixedNodes-1)*3+1;  % ux
    fixed_dof(2:3:end) = (fixedNodes-1)*3+2;  % uy
    fixed_dof(3:3:end) = (fixedNodes-1)*3+3;  % uz
elseif bc == 2  % one end clamped, another free
    fixedNodes  = unique([bottomNodes;nextBottomNodes]);
    fixed_dof   = zeros(length(fixedNodes)*3,1);
    fixed_dof(1:3:end) = (fixedNodes-1)*3+1;  % ux
    fixed_dof(2:3:end) = (fixedNodes-1)*3+2;  % uy
    fixed_dof(3:3:end) = (fixedNodes-1)*3+3;  % uz
elseif bc == 3  % one end clamped, another simply supported uy = 0
    fixedNodes1 = unique([bottomNodes;nextBottomNodes]);
    fixedNodes2 = unique(topNodes);
    fixed_dof1  = zeros(length(fixedNodes1)*3,1);
    fixed_dof1(1:3:end) = (fixedNodes1-1)*3+1;  % ux
    fixed_dof1(2:3:end) = (fixedNodes1-1)*3+2;  % uy
    fixed_dof1(3:3:end) = (fixedNodes1-1)*3+3;  % uz
    fixed_dof2  = (fixedNodes2-1)*3+2;  % uy
    fixed_dof   = [fixed_dof1; fixed_dof2];
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
M = sparse(nDofs,nDofs);     % mass matrix

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
    Me = zeros(nnElem,nnElem);
    for ipt = 1:size(gp,1)   % loop over integration points
        pt   = gp(ipt,1:2);      % reference parametric coordinates for each integration point
        wt   = wgt(ipt);       % weigths for each integration point
        gPts = parameterGaussMapping( elDoma, pt );   % gauss integration mapping
        j1   = jacobianGaussMapping( elDoma );        % jacobian value for gauss mapping   
        [R,dR,dR2] = computeTsplineBasis2ndDers([pu,pv],gPts,Ce,we);
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
        N   = zeros(3, nnElem);
        N(1,1:3:end)  =  R;
        N(2,2:3:end)  =  R;
        N(3,3:3:end)  =  R;
        Me  = Me + rho * h * (N' * N) * fac;
    end   
    K(sctrB,sctrB) = K(sctrB,sctrB) +  Ke;
    M(sctrB,sctrB) = M(sctrB,sctrB) +  Me;
end
toc;

%****************************************************************%
%                   Solve stiffness equation                     %
%****************************************************************%
fprintf("Solve stiffness equation \n");
tic
activedof = setdiff((1:nDofs), fixed_dof);
[modeShapes, freq] = eigs(K(activedof,activedof),M(activedof,activedof),nModes,1);
beta = diag(freq);
[sortbeta, index_beta] = sort(beta);
T = sqrt(sortbeta)./2./pi;
toc


%****************************************************************%
%                       Post-processing                          %
%****************************************************************%
fprintf("Post-processing \n");
tic
figure(2);clf;
set( gcf, 'Color', 'White', 'Unit', 'Normalized', 'Position', [0.1,0.1,0.6,0.6] );
title( 'Mode Shapes of Freeform Cylindrical Shell Using IGA', 'FontSize', 14', 'FontWeight', 'Bold') ;
axis off;
nCol = 5;  nRow = ceil( nModes / nCol );
rowH = 0.6 / nRow ;  colW = 0.8 / nCol;
colX = 0.06 + linspace( 0, 0.94, nCol+1);  
colX = colX(1:end-1) ;
rowY = 0.1 + linspace( 0.9, 0, nRow+1 ) ;  
rowY = rowY(2:end) ;

for i = 1:nModes     % Loop over each mode shapes
    U = zeros(nDofs,1);
    U(activedof) = modeShapes(:,index_beta(i));
    U      = reshape(U,3,[])';
    scale = 1;   % scale factor
    mesh2 = mesh;
    mesh2.coords(:,1) = mesh2.coords(:,1) + U(:,1)*scale;
    mesh2.coords(:,2) = mesh2.coords(:,2) + U(:,2)*scale;
    mesh2.coords(:,3) = mesh2.coords(:,3) + U(:,3)*scale;
    vmesh = getVisualMesh(mesh2,U);
    rowId = ceil( i / nCol ) ;
    colId = i - (rowId - 1) * nCol ;
    axes( 'Position', [colX(colId), rowY(rowId), colW, rowH] ) ;
    p = patch('Faces',vmesh.element, 'Vertices', vmesh.node);         % plot the surface
    axis equal
    axis off;
    cdata = sqrt(vmesh.u(:,1).^2+vmesh.u(:,2).^2+vmesh.u(:,3).^2);
    set(p,'FaceColor','interp','FaceVertexCData',cdata,'EdgeColor','none');
    hold on;
%     for j = 1:mesh.nElems
%         hold on;
%         plot3(vmesh.node(vmesh.linmesh{j,1},1),vmesh.node(vmesh.linmesh{j,1},2),vmesh.node(vmesh.linmesh{j,1},3),'k-');
%     end
    view(10,5);
    axis([0,1,0,L,0,1]);
    set(gcf,'renderer','opengl');
    light;
    set(gca,'xtick',[],'ytick',[],'ztick',[]);
    title( [sprintf( 'Modes %d', i), ', \omega = ', sprintf('%3.3f',beta(index_beta(i)))]) ;  
    mycolor = abaqusColorMap(12);
    colormap(mycolor);
    hold off;
end
toc

