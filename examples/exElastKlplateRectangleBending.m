
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%
%  An example for bending of square plate using tspline model based on
%  Kirchhoff-Love plate assumption
%  
%
%  ---------------------------------------
%
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%  - 24-JAN-2022
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
    filename = 'geoSquarePlateTspline.iga';
%     filename = 'geoSquarePlateUnstrucTspline.iga';
    geo  = readIgaFile(filename);
    mesh = buildIgaMesh( geo );
elseif strcmp(meshType, 'NURBS')
    pts = [0,0]; p = 5; q = p; elemX = 6; elemY = elemX; 
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
% plot scaled model
figure(1)
plotIgaMesh(mesh,1);
view(2);

%****************************************************************%
%               Define material properties and dof               %
%****************************************************************%
fprintf("Define material properties and dof \n");
E    = 2*10^8;         % Youngs modulus
nu   = 0.3;            % Poisson ratio
h    = 0.01;           % thickness of the plate
q0   = -10;           % uniform pressure
dof  = 1;
D0   = E*h^3/(12*(1-nu*nu));
D    = D0*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];

%****************************************************************%
%                   Enforce boundary conditions                  %
%****************************************************************%
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
hold on;
plot3(mesh.coords(leftNodes,1),mesh.coords(leftNodes,2),mesh.coords(leftNodes,3),'mo');
plot3(mesh.coords(nextLeftNodes,1),mesh.coords(nextLeftNodes,2),mesh.coords(nextLeftNodes,3),'m*');
plot3(mesh.coords(rightNodes,1),mesh.coords(rightNodes,2),mesh.coords(rightNodes,3),'bo');
plot3(mesh.coords(nextRightNodes,1),mesh.coords(nextRightNodes,2),mesh.coords(nextRightNodes,3),'b*');
plot3(mesh.coords(bottomNodes,1),mesh.coords(bottomNodes,2),mesh.coords(bottomNodes,3),'co');
plot3(mesh.coords(nextBottomNodes,1),mesh.coords(nextBottomNodes,2),mesh.coords(nextBottomNodes,3),'c*');
plot3(mesh.coords(topNodes,1),mesh.coords(topNodes,2),mesh.coords(topNodes,3),'ro');
plot3(mesh.coords(nextTopNodes,1),mesh.coords(nextTopNodes,2),mesh.coords(nextTopNodes,3),'r*');

% four boundaries are simply supported 
dbc = [];    % dbc = [node index, direction, prescribed displacement]
bc_nodes = (unique([bottomNodes;topNodes;leftNodes;rightNodes]));
dbc = [dbc; bc_nodes,     ones(size(bc_nodes)),   zeros(size(bc_nodes))];
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
K = sparse(ndofs,ndofs);   % stiffness matrix 
F = zeros(ndofs,1);        % external force matrix

% use gaussian integration rule
gp_x = max(mesh.elDegree(:,1))+1;           % number of integration points in x-direction
gp_y = max(mesh.elDegree(:,2))+1;           % number of integration points in y-direction
[gp, wgt] = gaussQuadrature(gp_x, gp_y);   % calculate integration points and its weights
elDoma = [0,1,0,1];

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
        j1     = jacobianGaussMapping( elDoma );            % jacobian value for gauss mapping   
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
        xt  = R*elCpts(:,1:2);
        q2  = q0*sin(pi*xt(1)/a)*sin(pi*xt(2)/b);
        F(sctr) = F(sctr) + q2 * R'* fac;
        K(sctr,sctr) = K(sctr,sctr) + B' * D * B * fac;
    end   
end
toc

%****************************************************************%
%                   Solve stiffness equation                     %
%****************************************************************%
fprintf("Solve stiffness equation \n");
tic
ndbc = size(dbc,1);   % number of displacement constrained nodes
K(scatdbc,:) = zeros(ndbc, ndofs);
K(scatdbc,scatdbc) = eye(ndbc);
F(scatdbc,:) = 0;
F(scatdbc,:) = dbc(:,3);        
u = K\F;
toc

%****************************************************************%
%                     Convergence study                          %
%****************************************************************%
fprintf("Convergence study: \n");
dw = 0; dMx = 0; ew = 0; eMx = 0;
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
        j1     = jacobianGaussMapping( elDoma );            % jacobian value for gauss mapping   
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
        
        strain = - B * u(sctr);
        stress = D * strain;
        edsp   = u(sctr);
        w      = R * edsp;
        Mx     = stress(1);
        xy     = R*elCpts(:,1:2); 
        x      = xy(1); y = xy(2);
        tmp     = (1/a^2 + 1/b^2)^2;
        exactw  = q0/pi^4/D0/tmp*sin(pi*x/a)*sin(pi*y/b);
        exactMx = q0/pi^2/tmp*(1/a^2+nu/b^2)*sin(pi*x/a)*sin(pi*y/b);
        fac = j1 *j2 * wt;      
        dw  = dw + (w-exactw)^2*fac;
        ew  = ew + exactw^2*fac;
        dMx = dMx + (Mx-exactMx)^2*fac;
        eMx = eMx + exactMx^2*fac;
    end   
end
fprintf('Deflection error is %12.6e.\n',sqrt(dw/ew));
fprintf('Bending moment Mx error is %12.6e.\n',sqrt(dMx/eMx));


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
    pu     = mesh.elDegree(e,1);  % degree u
    pv     = mesh.elDegree(e,2);  % degree v
    we     = mesh.coords(sctr,4); % Tspline control points' weights
    Ce     = mesh.elExtOpe{e,1};  % elemental extration operator
    for j = 1:(num1+1)*(num2+1)
        pt  = vmesh.tripts(j,:);  
        [R,dR,dR2] = computeTsplineBasis2ndDers([pu,pv],pt,Ce,we);
        jmat1 = dR*elCpts(:,1:2); 
        jmat2 = dR2*elCpts(:,1:2); 
        dxdxi  = jmat1(1,1); dydxi  = jmat1(1,2);
        dxdeta = jmat1(2,1); dydeta = jmat1(2,2);
        jmat3  = [     dxdxi^2      dydxi^2              2*dxdxi*dydxi;
                      dxdeta^2      dydeta^2           2*dxdeta*dydeta;
                  dxdxi*dxdeta  dydxi*dydeta  dxdxi*dydeta+dxdeta*dydxi];
        dRdx   = jmat1\dR;       
        dR2dx  = jmat3\(dR2-jmat2*dRdx);
        B          = dR2dx;
        B(3,:)     = B(3,:)*2;   
        strain = - B * u(sctr);
        stress = D * strain;
        edsp   = u(sctr);
        vmesh.displacement(count,1) = R * edsp;
        vmesh.node(count,1:2)       = R*elCpts(:,1:2);
        vmesh.stress(count,:)       = stress';
        vmesh.node(count,3)   = 0;
        count = count+1;
    end    
end


figure(2);
suptitle( 'Static Bending of a Rectangular Plate Using KL based IGA') ;
for gg = 1:4
    subplot(2,2,gg);
    p = patch('Faces',vmesh.element, 'Vertices', vmesh.node); 
    if gg == 1, cdata = vmesh.displacement(:,1); title('Deflection','FontName','Times New Roman');
    elseif gg == 2, cdata = vmesh.stress(:,1); title('Moment stress Mx','FontName','Times New Roman');
    elseif gg == 3, cdata = vmesh.stress(:,2); title('Moment stress My','FontName','Times New Roman');
    elseif gg == 4, cdata = vmesh.stress(:,3); title('Moment stress Mxy','FontName','Times New Roman');
    end
    set(p,'FaceColor','interp','FaceVertexCData',cdata);
    set(p,'EdgeColor','none');
    hold on;
    for i = 1:mesh.nElems
        plot3(vmesh.node(vmesh.linmesh{i,1},1),vmesh.node(vmesh.linmesh{i,1},2),vmesh.node(vmesh.linmesh{i,1},3),'k-');
    end
    C = colorbar;
    colorMin = C.Limits(1);
    colorMax = C.Limits(2);
    set(C,'Ticks',linspace(colorMin,colorMax,7));
    mycolor = abaqusColorMap(24);
    colormap(mycolor);
    axis equal; 
    axis off;
end

