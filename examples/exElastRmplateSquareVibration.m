%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%
%  An example for vibration of square plate using tspline model
%
%  ---------------------------------------
%
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%  - 29-NOV-2021
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
tic
% MATERIAL PREPERTIES
E      = 200*10^9;   % Youngs modulus
nu     = 0.3;        % Poisson ratio
h      = 0.1;        % thickness of the plate
rho    = 8000;       % density
nModes = 15;         % number of mode shapes
dof    = 3;          % degree of freedom
bcs = [3 3 3 3];     % [left boundary, bottom boundary, right boundary, top bounday]
                     % 1-free, 2-simply supported,  3-clamped/fixed imposed
                     % on corresponding boundaries
                     % example: [2 3 2 3] means that the left and right boundaries are
                     % simply supported, bottom and top boundaries are clamped
                     % [2 3 2 3] is be denoted by 'SCSC'
                     

%****************************************************************%
%                   Enforce boundary conditions                  %
%****************************************************************%
fprintf("Enforce boundary conditions \n");
% find four boundaries
tol = 1e-3;
left_nodes   = find( abs(mesh.coords(:,1)) < tol );
bottom_nodes = find( abs(mesh.coords(:,2)) < tol );
right_nodes  = find( abs(mesh.coords(:,1)-a) < tol );
top_nodes    = find( abs(mesh.coords(:,2)-b) < tol );

% plot constrained control points
hold on;
plot3(mesh.coords(left_nodes,1),mesh.coords(left_nodes,2),mesh.coords(left_nodes,3),'bo');
plot3(mesh.coords(right_nodes,1),mesh.coords(right_nodes,2),mesh.coords(right_nodes,3),'bo');
plot3(mesh.coords(bottom_nodes,1),mesh.coords(bottom_nodes,2),mesh.coords(bottom_nodes,3),'ro');
plot3(mesh.coords(top_nodes,1),mesh.coords(top_nodes,2),mesh.coords(top_nodes,3),'ro');

bc_nodes = cell(1,4);
bc_nodes{1} = left_nodes;
bc_nodes{2} = bottom_nodes;
bc_nodes{3} = right_nodes;
bc_nodes{4} = top_nodes;

fixed_dof = [];  %  constrained dofs
for i = 1:2:3    % left boundary and right boundary
    if bcs(i) == 2
        fixed_dof = [fixed_dof; bc_nodes{i}*3-2; bc_nodes{i}*3-1 ];
    elseif bcs(i) == 3
        fixed_dof = [fixed_dof; bc_nodes{i}*3-2; bc_nodes{i}*3-1; bc_nodes{i}*3 ];
    end
end
for i = 2:2:4    % bottom boundary and top boundary
    if bcs(i) == 2
        fixed_dof = [fixed_dof; bc_nodes{i}*3-2; bc_nodes{i}*3 ];
    elseif bcs(i) == 3
        fixed_dof = [fixed_dof; bc_nodes{i}*3-2; bc_nodes{i}*3-1; bc_nodes{i}*3 ];
    end
end
fixed_dof = unique(fixed_dof);


%****************************************************************%
%           Initialize and assemble stiffness matrix             %
%****************************************************************%
fprintf("Initialize and assemble stiffness matrix \n");
tic
% initialize stiffness and force matrices
nDofs = dof * mesh.nCpts;    % total dofs
K = sparse(nDofs,nDofs);     % stiffness matrix 
F = zeros(nDofs,1);          % external force matrix
M = sparse(nDofs,nDofs);     % mass matrix

% Set the D matrix
lamda = 6/5;
D0 = E*h^3/(12*(1-nu*nu));
Db =D0*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];
Ds =E*h/(2*(1+nu)*lamda)*[1 0;0 1];
D  = {Db, zeros(3,2);zeros(2,3), Ds};
D  = cell2mat(D);

% use gaussian integration rule
gp_x = max(mesh.elDegree(:,1))+1;           % number of integration points in x-direction
gp_y = max(mesh.elDegree(:,2))+1;           % number of integration points in y-direction
[gp, wgt] = gaussQuadrature(gp_x, gp_y);   % calculate integration points and its weights
elDoma = [0,1,0,1];
for e = 1:mesh.nElems                 % loop over elements
    sctr   = mesh.elNodeCnt{e,:};     % element control points index
    elCpts = mesh.coords(sctr,:);     % coordinates of element control points
    nn     = numel(sctr);             % number of control points for each element
    nnElem = nn*dof;                  % dof for each element
    sctrB = zeros(1, nnElem); 
    for i = 1:dof
        sctrB(i:dof:nnElem) = dof*(sctr-1) + i;  % displacement in i-th direction
    end
    sctrw = 3*sctr-2;
    p  = mesh.elDegree(e,1);
    q  = mesh.elDegree(e,2);
    Ce = mesh.elExtOpe{e,1};
    we = mesh.coords(sctr,4); % Tspline control points' weights
    for ipt = 1:size(gp,1)   % loop over integration points
        pt = gp(ipt,:);      % reference parametric coordinates for each integration point
        wt = wgt(ipt);       % weigths for each integration point
        gauPts  = parameterGaussMapping( elDoma, pt );   % gauss integration mapping
        j1      = jacobianGaussMapping( elDoma );            % jacobian value for gauss mapping   
        [R,dR]  = computeTsplineBasisDers([p,q],gauPts,Ce,we);
        jmatrix = dR*elCpts(:,1:2); 
        j2 = det(jmatrix);
        ders =  jmatrix \ dR;  
        
        B  = zeros(5, nnElem);
        N  = zeros(3, nnElem);
        B(1,3:3:nnElem) = ders(1,:);    
        B(2,2:3:nnElem) = -ders(2,:);   
        B(3,2:3:nnElem) = -ders(1,:);   
        B(3,3:3:nnElem) = ders(2,:);
        B(4,1:3:nnElem) = -ders(2,:);  
        B(4,2:3:nnElem) = R;
        B(5,1:3:nnElem) = ders(1,:);   
        B(5,3:3:nnElem) = R;
        
        N(1,1:3:nnElem) = R;
        N(2,2:3:nnElem) = R;
        N(3,3:3:nnElem) = R;
        m = rho*[h  0  0;  0  h^3/12  0;  0  0  h^3/12];
        fac = j1 *j2 * wt;      
        K(sctrB,sctrB) = K(sctrB,sctrB) + B' * D * B * fac;   
        M(sctrB,sctrB) = M(sctrB,sctrB) + N' * m * N * fac;
    end   
end
toc

%****************************************************************%
%                   Solve stiffness equation                     %
%****************************************************************%
fprintf("Solve stiffness equation \n");
tic
activedof = setdiff((1:nDofs), fixed_dof);
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
title( sprintf('Mode Shapes of RM Square Plate Using IGA with h = %2.4f',h), 'FontSize', 14', 'FontWeight', 'Bold','FontName','Times New Roman') ;
axis off;
nCol = 5;  nRow = ceil( nModes / nCol );
rowH = 0.6 / nRow ;  colW = 0.8 / nCol;
colX = 0.06 + linspace( 0, 0.94, nCol+1);  
colX = colX(1:end-1) ;
rowY = 0.1 + linspace( 0.9, 0, nRow+1 ) ;  
rowY = rowY(2:end) ;

for i = 1:nModes     % Loop over each mode shapes
    U = zeros(nDofs,1);
    U(activedof) = real(modeShapes(:,index_beta(i)));
    scale = 1.0/max(abs(U(activedof)));   % scale factor
    mesh2 = mesh;
    mesh2.coords(:,3) = U(1:3:end)*scale;
    vmesh = getVisualMesh(mesh2);
    rowId = ceil( i / nCol ) ;
    colId = i - (rowId - 1) * nCol ;
    axes( 'Position', [colX(colId), rowY(rowId), colW, rowH] ) ;
    p = patch('Faces',vmesh.element, 'Vertices', vmesh.node);         % plot the surface
    box on;
    cdata = vmesh.node(:,3);
    set(p,'FaceColor','interp','FaceVertexCData',cdata,'EdgeColor','none');
    axis equal;
    hold on;
    for j = 1:mesh.nElems
        hold on;
        plot3(vmesh.node(vmesh.linmesh{j,1},1),vmesh.node(vmesh.linmesh{j,1},2),vmesh.node(vmesh.linmesh{j,1},3),'k-');
    end
    view(2);
    set(gca,'xtick',[],'ytick',[],'ztick',[]);
    title( [sprintf( 'Modes %d', i), ', \omega = ', sprintf('%3.3f',beta(index_beta(i)))],'FontName','Times New Roman') ;  
    mycolor = abaqusColorMap(12);
    colormap(mycolor);
    hold off;
end
toc







