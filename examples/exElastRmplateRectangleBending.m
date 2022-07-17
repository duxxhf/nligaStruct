
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%
%  An example for bending of a simply supported rectangle plate
%  
%  Please see <Du, X., Zhao, G., & Wang, W. (2015). Nitsche method for 
%  isogeometric analysis of Reissner¨CMindlin plate with non-conforming 
%  multi-patches. Computer Aided Geometric Design, 35, 121-136.> for
%  problems' definition.
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
a = 1; b = 1.5;
if strcmp(meshType, 'T-SPLINE')
%     filename = 'geoSquarePlateUnstrucTspline.iga';
    filename = 'geoSquarePlateTspline.iga';
    geo  = readIgaFile(filename);
    mesh = buildIgaMesh( geo );
elseif strcmp(meshType, 'NURBS')
    pts = [0,0]; p = 2; q = p; elemX = 6; elemY = elemX; 
    geo  = geoRectangularPlate( a,b,pts,p,q,elemX,elemY );
    mesh = buildIgaMesh( geo );
end
maxNodeX = max(mesh.coords(:,1));
minNodeX = min(mesh.coords(:,1));
maxNodeY = max(mesh.coords(:,2));
minNodeY = min(mesh.coords(:,2));
mesh.coords(:,1) = a*(mesh.coords(:,1) - minNodeX)/( maxNodeX - minNodeX );
mesh.coords(:,2) = b*(mesh.coords(:,2) - minNodeY)/( maxNodeY - minNodeY );
toc

%****************************************************************%
%               Define material properties and dof               %
%****************************************************************%
fprintf("Define material properties and dof \n");
E    = 2*10^8;         % Youngs modulus
nu   = 0.3;            % Poisson ratio
h    = 0.01;           % thickness of the plate
q0   = -10;           % uniform pressure
dof  = 3;

% Set the D matrix
lamda = 6/5;
D0 = E*h^3/(12*(1-nu*nu));
Db = D0*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];
Ds = E*h/(2*(1+nu)*lamda)*[1 0;0 1];
D  = {Db, zeros(3,2);zeros(2,3), Ds};
D  = cell2mat(D);

% plot scaled model
figure(1)
plotIgaMesh(mesh,1);
view(2);

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
bc_nodes = (unique([bottom_nodes;left_nodes;top_nodes;right_nodes]));
% plot constrained control points
hold on;
plot3(mesh.coords(left_nodes,1),mesh.coords(left_nodes,2),mesh.coords(left_nodes,3),'bo');
plot3(mesh.coords(right_nodes,1),mesh.coords(right_nodes,2),mesh.coords(right_nodes,3),'bo');
plot3(mesh.coords(bottom_nodes,1),mesh.coords(bottom_nodes,2),mesh.coords(bottom_nodes,3),'ro');
plot3(mesh.coords(top_nodes,1),mesh.coords(top_nodes,2),mesh.coords(top_nodes,3),'ro');

% four boundaries are simply supported
dbc = [];    % dbc = [node index, direction, prescribed displacement]
dbc = [dbc; bc_nodes,   ones(size(bc_nodes)),   zeros(size(bc_nodes))];
scatdbc = [];
scattbc = [];
if ~isempty(dbc)
    scatdbc = dof * (dbc(:,1)-1) + dbc(:,2);   % scatter dbc
end

%****************************************************************%
%           Initialize and assemble stiffness matrix             %
%****************************************************************%
fprintf("Initialize and assemble stiffness matrix \n");
tic
% initialize stiffness and force matrices
ndofs = dof * mesh.nCpts;    % total dofs
k = sparse(ndofs,ndofs);   % stiffness matrix 
f = zeros(ndofs,1);        % external force matrix

% use gaussian integration rule
gp_x = max(mesh.elDegree(:,1))+1;           % number of integration points in x-direction
gp_y = max(mesh.elDegree(:,2))+1;           % number of integration points in y-direction
[gp, wgt] = gaussQuadrature(gp_x, gp_y);   % calculate integration points and its weights
elDoma = [0,1,0,1];
for e = 1:mesh.nElems                 % loop over elements
    sctr   = mesh.elNodeCnt{e,:};       % element control points index
    elCpts = mesh.coords(sctr,:);     % coordinates of element control points
    nn     = numel(sctr);                 % number of control points for each element
    nnElem = nn*dof;                  % dof for each element
    sctrB = zeros(1, nnElem); 
    for i = 1:dof
        sctrB(i:dof:nnElem) = dof*(sctr-1) + i;  % displacement in i-th direction
    end
    sctrw = 3*sctr-2;
    B  = zeros(5, nnElem);
    N  = zeros(3, nnElem);
    p  = mesh.elDegree(e,1);
    q  = mesh.elDegree(e,2);
    Ce = mesh.elExtOpe{e,1};
    we = mesh.coords(sctr,4); % Tspline control points' weights
    for ipt = 1:size(gp,1)   % loop over integration points
        pt = gp(ipt,:);      % reference parametric coordinates for each integration point
        wt = wgt(ipt);       % weigths for each integration point
        gauPts = parameterGaussMapping( elDoma, pt );   % gauss integration mapping
        j1 = jacobianGaussMapping( elDoma );            % jacobian value for gauss mapping   
        [R,dR] = computeTsplineBasisDers([p,q],gauPts,Ce,we);
        jmatrix = dR*elCpts(:,1:2); 
        j2 = det(jmatrix);
        ders =  jmatrix \ dR;     
        B(1,3:3:nnElem) = ders(1,:);    
        B(2,2:3:nnElem) = -ders(2,:);   
        B(3,2:3:nnElem) = -ders(1,:);   
        B(3,3:3:nnElem) = ders(2,:);
        B(4,1:3:nnElem) = ders(1,:);  
        B(4,3:3:nnElem) = R;
        B(5,1:3:nnElem) = -ders(2,:);   
        B(5,2:3:nnElem) = R;
        
        fac = j1 *j2 * wt;      
        xy  = R*elCpts(:,1:2);
        q2  = q0*sin(pi*xy(1)/a)*sin(pi*xy(2)/b);
        f(sctrw) = f(sctrw) + q2 * R'* fac;
        k(sctrB,sctrB) = k(sctrB,sctrB) + B' * D * B * fac;   
    end   
end
toc

%****************************************************************%
%                   Solve stiffness equation                     %
%****************************************************************%
fprintf("Solve stiffness equation \n");
tic
ndbc = size(dbc,1);   % number of displacement constrained nodes
k(scatdbc,:) = zeros(ndbc, ndofs);
k(scatdbc,scatdbc) = eye(ndbc);
f(scatdbc,:) = 0;
f(scatdbc,:) = dbc(:,3);        
u = k\f;
toc

%****************************************************************%
%                       Post-processing                          %
%****************************************************************%
fprintf("Post-processing \n");
tic
% post-processing
num1 = 5;  num2 = 5;  
vmesh = buildVisualMesh(num1,num2,mesh);
count = 1;
elDoma = [0,1,0,1];
for e = 1:mesh.nElems
    sctr   = mesh.elNodeCnt{e,1};       % element control points index
    elCpts = mesh.coords(sctr,1:3);     % coordinates of element control points
    nn     = numel(sctr);               % number of control points for each element
    nnElem = nn*dof;               % dof for each element
    sctrB  = zeros(1, nnElem); 
    for i  = 1:dof
        sctrB(i:dof:nnElem) = dof*(sctr-1) + i;  % displacement in i-th direction
    end
    sctrB = zeros(1, nnElem); 
    sctrB(1:3:nnElem) = 3*sctr-2;  % deflection   
    sctrB(2:3:nnElem) = 3*sctr-1;  % rotation x   
    sctrB(3:3:nnElem) = 3*sctr;    % rotation y
    
    p      = mesh.elDegree(e,1);  % degree u
    q      = mesh.elDegree(e,2);  % degree v
    we     = mesh.coords(sctr,4); % Tspline control points' weights
    Ce     = mesh.elExtOpe{e,1};  % elemental extration operator
    
    for j = 1:(num1+1)*(num2+1)
        pt  = vmesh.tripts(j,:);  
        [R,dR] = computeTsplineBasisDers([p,q],pt,Ce,we);
        jmatrix = dR*elCpts(:,1:2); 
        ders =  jmatrix \ dR;    
        B    = zeros(5,nnElem);
        B(1,3:3:nnElem) = ders(1,:);    
        B(2,2:3:nnElem) = -ders(2,:);   
        B(3,2:3:nnElem) = -ders(1,:);   
        B(3,3:3:nnElem) = ders(2,:);
        B(4,1:3:nnElem) = ders(1,:);  
        B(4,3:3:nnElem) = R;
        B(5,1:3:nnElem) = -ders(2,:);   
        B(5,2:3:nnElem) = R;
        strain = B * u(sctrB);
        stress = D * strain;
        edsp   = u(sctrB);
        edsp   = reshape(edsp, 3, nn);
        vmesh.displacement(count,:) = R * edsp';
        vmesh.node(count,1:2)       = R*elCpts(:,1:2);
        vmesh.stress(count,:)       = stress';
        vmesh.node(count,3)   = 0;
        count = count+1;
    end    
end

figure(2);
title( 'Static Bending of a Rectangular Plate Using RM based IGA', 'FontSize', 14', 'FontWeight', 'Bold','FontName','Times New Roman') ;
axis off;
colX = linspace( 0, 0.7, 4 );
rowY = [0.55, 0.15];
for k = 1:8
    rowId = ceil(k/4);
    colId = k - (rowId-1)*4;
    axes( 'Position', [colX(colId), rowY(rowId), 0.3, 0.3] ) ;
    p = patch('Faces',vmesh.element, 'Vertices', vmesh.node); 
    if k == 1,  cdata = vmesh.displacement(:,1); title('\itw','FontName','Times New Roman');
    elseif k == 2,  cdata = vmesh.displacement(:,2);  title('\itR_x','FontName','Times New Roman');
    elseif k == 3,  cdata = vmesh.displacement(:,3);  title('\itR_y','FontName','Times New Roman');
    elseif k == 4,  cdata = vmesh.stress(:,1);  title('\itM_x','FontName','Times New Roman');  
    elseif k == 5,  cdata = vmesh.stress(:,2);  title('\itM_y','FontName','Times New Roman');  
    elseif k == 6,  cdata = vmesh.stress(:,3);  title('\itM_{xy}','FontName','Times New Roman');  
    elseif k == 7,  cdata = vmesh.stress(:,4);  title('\itq_x','FontName','Times New Roman');  
    elseif k == 8,  cdata = vmesh.stress(:,5);  title('\itq_y','FontName','Times New Roman');  
    end
    set(p,'FaceColor','interp','FaceVertexCData',cdata);
    set(p,'EdgeColor','none');
    hold on;
    for i = 1:mesh.nElems
        plot3(vmesh.node(vmesh.linmesh{i,1},1),vmesh.node(vmesh.linmesh{i,1},2),vmesh.node(vmesh.linmesh{i,1},3),'k-');
    end
    axis equal; 
    axis off;
    hold off;
    C = colorbar;
    colorMin = C.Limits(1);
    colorMax = C.Limits(2);
    set(C,'Ticks',linspace(colorMin,colorMax,7));
    mycolor = abaqusColorMap(24);
    colormap(mycolor);
    axis off;
end

