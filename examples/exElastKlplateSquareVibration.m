
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%
%  An example for vibration of square plate using tspline model based on
%  Kirchhoff-Love plate assumption
%  
%
%  ---------------------------------------
%
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%  - 29-JAN-2022
%
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

clc; 
clear;
%****************************************************************%
%                        Import model                            %
%****************************************************************%
fprintf("Import geometrical model \n");
tic
meshType = 'T-SPLINE';  % or 'NURBS'
% meshType = 'NURBS';  % or 'NURBS'
a = 1; b = 1;
if strcmp(meshType, 'T-SPLINE')
%     filename = 'geoSquarePlateUnstrucTspline.iga';
    filename = 'geoSquarePlateTspline.iga';
    geo  = readIgaFile(filename);
    mesh = buildIgaMesh( geo );
elseif strcmp(meshType, 'NURBS')
    pts = [0,0]; p = 3; q = p; elemX = 20; elemY = elemX; 
    geo  = geoRectangularPlate( a,b,pts,p,q,elemX,elemY );
    mesh = buildIgaMesh( geo );
end
maxNodeX = max(mesh.coords(:,1));
minNodeX = min(mesh.coords(:,1));
maxNodeY = max(mesh.coords(:,2));
minNodeY = min(mesh.coords(:,2));
mesh.coords(:,1) = a*(mesh.coords(:,1) - minNodeX)/( maxNodeX - minNodeX );
mesh.coords(:,2) = b*(mesh.coords(:,2) - minNodeY)/( maxNodeY - minNodeY );
% plot scaled model
figure(1)
plotIgaMesh(mesh,1);
view(2);
toc

%****************************************************************%
%               Define material properties and dof               %
%****************************************************************%
fprintf("Define material properties and dof \n");
E    = 200*10^9;       % Youngs modulus
nu   = 0.3;            % Poisson ratio
h    = 0.001;          % thickness of the plate
dof  = 1;              % degree if freedom for each control point
rho  = 8000;           % density 
nModes = 15;           % number of mode shapes
D0   = E*h^3/(12*(1-nu*nu));
D    = D0*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];
bcs  = 2; % 1-simply supported, =2 clamped
%****************************************************************%
%                   Enforce boundary conditions                  %
%****************************************************************%
fprintf("Enforce boundary conditions \n");
% find four boundaries
tol = 1e-3;
fprintf("Enforce boundary conditions \n");
tic
% find boundary nodes
for i = 1:length(mesh.nodeSets)
    if strcmp(mesh.nodeSets{i,1}.name,'bottom')
        bottomNodes = mesh.nodeSets{i,1}.gloInx;
    elseif strcmp(mesh.nodeSets{i,1}.name,'nextBottom')
        nextBottomNodes = mesh.nodeSets{i,1}.gloInx;
    elseif strcmp(mesh.nodeSets{i,1}.name,'top')
        topNodes = mesh.nodeSets{i,1}.gloInx;
    elseif strcmp(mesh.nodeSets{i,1}.name,'nextTop')
        nextTopNodes = mesh.nodeSets{i,1}.gloInx;
    elseif strcmp(mesh.nodeSets{i,1}.name,'left')
        leftNodes =mesh.nodeSets{i,1}.gloInx;
    elseif strcmp(mesh.nodeSets{i,1}.name,'nextLeft')
        nextLeftNodes = mesh.nodeSets{i,1}.gloInx;
    elseif strcmp(mesh.nodeSets{i,1}.name,'right')
        rightNodes = mesh.nodeSets{i,1}.gloInx;
    elseif strcmp(mesh.nodeSets{i,1}.name,'nextRight')
        nextRightNodes = mesh.nodeSets{i,1}.gloInx;
    end
end


% plot constrained control points
% plot constrained control points
hold on;
plot3(mesh.coords(leftNodes,1),mesh.coords(leftNodes,2),mesh.coords(leftNodes,3),'mo');
plot3(mesh.coords(nextLeftNodes,1),mesh.coords(nextLeftNodes,2),mesh.coords(nextLeftNodes,3),'m*');
plot3(mesh.coords(rightNodes,1),mesh.coords(rightNodes,2),mesh.coords(rightNodes,3),'bo');
plot3(mesh.coords(nextRightNodes,1),mesh.coords(nextRightNodes,2),mesh.coords(nextRightNodes,3),'b*');
plot3(mesh.coords(bottomNodes,1),mesh.coords(bottomNodes,2),mesh.coords(bottomNodes,3),'co');
plot3(mesh.coords(nextBottomNodes,1),mesh.coords(nextBottomNodes,2),mesh.coords(nextBottomNodes,3),'c*');
plot3(mesh.coords(topNodes,1),mesh.coords(topNodes,2),mesh.coords(topNodes,3),'ro');
plot3(mesh.coords(nextTopNodes,1),mesh.coords(nextTopNodes,2),mesh.coords(nextTopNodes,3),'r*');

% four boundaries are clamped/fixed 
dbc = [];    % dbc = [node index, direction, prescribed displacement]

if bcs == 1
    fixedNodes = unique([bottomNodes;topNodes;leftNodes;rightNodes]);
elseif bcs == 2
    fixedNodes = unique([bottomNodes;nextBottomNodes;topNodes;nextTopNodes; ...
                        leftNodes;nextLeftNodes;rightNodes;nextRightNodes]);
end
fixed_dof  = fixedNodes;
%****************************************************************%
%           Initialize and assemble stiffness matrix             %
%****************************************************************%
fprintf("Initialize and assemble stiffness matrix \n");
tic
% initialize stiffness and force matrices
ndofs = dof * mesh.nCpts;    % total dofs
K = sparse(ndofs,ndofs);     % stiffness matrix 
M = sparse(ndofs,ndofs);     % mass matrix

% use gaussian integration rule
gp_x = max(mesh.elDegree(:,1))+1;           % number of integration points in x-direction
gp_y = max(mesh.elDegree(:,2))+1;           % number of integration points in y-direction
[gp, wgt] = gaussQuadrature(gp_x, gp_y);   % calculate integration points and its weights
elDoma    = [0,1,0,1];

for e = 1:mesh.nElems                 % loop over elements
    sctr   = mesh.elNodeCnt{e,:};       % element control points index
    elCpts = mesh.coords(sctr,:);     % coordinates of element control points
    pu  = mesh.elDegree(e,1);
    pv  = mesh.elDegree(e,2);
    Ce = mesh.elExtOpe{e,1};
    we = mesh.coords(sctr,4); % Tspline control points' weights
    for ipt = 1:size(gp,1)   % loop over integration points
        pt = gp(ipt,:);      % reference parametric coordinates for each integration point
        wt = wgt(ipt);       % weigths for each integration point
        gauPts = parameterGaussMapping( elDoma, pt );   % gauss integration mapping
        j1 = jacobianGaussMapping( elDoma );            % jacobian value for gauss mapping   
        [R,dR,dR2] = computeTsplineBasis2ndDers([pu,pv],gauPts,Ce,we);
        jmat1 = dR*elCpts(:,1:2); 
        jmat2 = dR2*elCpts(:,1:2); 
        j2 = det(jmat1);
        dxdxi  = jmat1(1,1); dydxi  = jmat1(1,2);
        dxdeta = jmat1(2,1); dydeta = jmat1(2,2);
        jmat3  = [     dxdxi^2      dydxi^2              2*dxdxi*dydxi;
                      dxdeta^2      dydeta^2           2*dxdeta*dydeta;
                  dxdxi*dxdeta  dydxi*dydeta  dxdxi*dydeta+dxdeta*dydxi];
        dRdx   = jmat1\dR;       
        dR2dx  = jmat3\(dR2-jmat2*dRdx);
        B          = dR2dx;
        B(3,:)     = B(3,:)*2;   
        fac = j1 *j2 * wt;      
        K(sctr,sctr) = K(sctr,sctr) + B' * D * B * fac;
        N = [R; dR];
        m = rho*[h,0,0; 0,h^3/12,0; 0,0,h^3/12];
		M(sctr,sctr) = M(sctr,sctr) + N' * m * N * fac;
    end   
end
toc

%****************************************************************%
%                   Solve stiffness equation                     %
%****************************************************************%
fprintf("Solve stiffness equation \n");
tic
activedof = setdiff((1:ndofs), fixed_dof);
[modeShapes, freq] = eigs(K(activedof,activedof),M(activedof,activedof),nModes,0);
beta = abs(diag(freq)).^0.5.*a^2/pi^2*(rho*h/D0)^0.5;
[sortbeta, index_beta] = sort(beta);
toc

%****************************************************************%
%                       Post-processing                          %
%****************************************************************%
fprintf("Post-processing \n");
tic
figure(2);clf;
set( gcf, 'Color', 'White', 'Unit', 'Normalized', 'Position', [0.1,0.1,0.6,0.6] );
title( sprintf('Mode Shapes of KL Square Plate Using T-spline based IGA with h = %2.4f',h), 'FontSize', 14', 'FontWeight', 'Bold','FontName','Times New Roman') ;
axis off;
nCol = 5;  nRow = ceil( nModes / nCol );
rowH = 0.6 / nRow ;  colW = 0.8 / nCol;
colX = 0.06 + linspace( 0, 0.94, nCol+1);  
colX = colX(1:end-1) ;
rowY = 0.1 + linspace( 0.9, 0, nRow+1 ) ;  
rowY = rowY(2:end) ;

for i = 1:nModes     % Loop over each mode shapes
    U = zeros(ndofs,1);
    U(activedof) = real(modeShapes(:,index_beta(i)));
    scale = 0;   % scale factor
    mesh2 = mesh;
    mesh2.coords(:,3) = mesh2.coords(:,3) + U*scale;
    vmesh = getVisualMesh(mesh2,U);
    rowId = ceil( i / nCol ) ;
    colId = i - (rowId - 1) * nCol ;
    axes( 'Position', [colX(colId), rowY(rowId), colW, rowH] ) ;
    p = patch('Faces',vmesh.element, 'Vertices', vmesh.node);         % plot the surface
    axis equal
    axis off;
    cdata = vmesh.u;
    set(p,'FaceColor','interp','FaceVertexCData',cdata,'EdgeColor','none');
    hold on;
    for j = 1:mesh.nElems
        hold on;
        plot3(vmesh.node(vmesh.linmesh{j,1},1),vmesh.node(vmesh.linmesh{j,1},2),vmesh.node(vmesh.linmesh{j,1},3),'k-');
    end
    view(2);
    set(gca,'xtick',[],'ytick',[],'ztick',[]);
    title( [sprintf( 'Modes %d', i), ', \omega = ', sprintf('%3.4f',beta(index_beta(i)))], 'FontName','Times New Roman') ;  
    mycolor = abaqusColorMap(12);
    colormap(mycolor);
    hold off;
end
toc

