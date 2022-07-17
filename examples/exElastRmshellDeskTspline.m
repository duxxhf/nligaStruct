%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%
%  An example for structural analysis of the benchmark problem 
%  - Desk model built with T-spline
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
filename = 'geoDesk.iga';
geo  = readIgaFile(filename);
mesh = buildIgaMesh( geo );
% plot desk model
figure(1)
vmesh = getVisualMesh(mesh);
p = patch('Faces',vmesh.element, 'Vertices', vmesh.node);         % plot the surface
set(p,'FaceColor','y','EdgeColor','none');
% alpha(0.5);
axis equal;
for i = 1:mesh.nElems
    hold on;
    plot3(vmesh.node(vmesh.linmesh{i,1},1),vmesh.node(vmesh.linmesh{i,1},2),vmesh.node(vmesh.linmesh{i,1},3),'k-');
end
box on;
xlabel('x'); ylabel('y'); zlabel('z');
view(37,37);
toc

%****************************************************************%
%               Define material properties and dof               %
%****************************************************************%
fprintf("Define material properties and dof \n");
tic
% material properties
nu   = 0.3;             % Poisson ratio
E    = 3.0e11;          % Youngs modulus 
kapa = 5/6;             % Shear correction factor
h    = 0.03;            % thickness   
Dl = E/(1-nu^2)*[  1, nu, 0, 0, 0, 0;
                   nu, 1, 0, 0, 0, 0;
                   0,  0, 0, 0, 0, 0;
                   0, 0, 0, (1-nu)/2, 0, 0;
                   0, 0, 0, 0, kapa*(1-nu)/2, 0;
                   0, 0, 0, 0, 0, kapa*(1-nu)/2  ];  % elasticity matrix
dof = 6;             % degree of freedom
toc


%****************************************************************%
%                   Enforce boundary conditions                  %
%****************************************************************%
fprintf("Enforce boundary conditions \n");
tic
% find four boundaries
tol = 1e-3;
bottomNodes   = find( mesh.coords(:,3) < -63.5 );
% plot constrained control points
hold on;
plot3(mesh.coords(bottomNodes,1),mesh.coords(bottomNodes,2),mesh.coords(bottomNodes,3),'co');

% impose displacement boundary conditions
dbc = [];    % dbc = [node index, direction, prescribed displacement]
% rigid diaphram boundary, ux = uy = uz = theta_x = theta_y = theta_z = 0
dbc = [dbc; bottomNodes,     1*ones(size(bottomNodes)),     zeros(size(bottomNodes))];  
dbc = [dbc; bottomNodes,     2*ones(size(bottomNodes)),     zeros(size(bottomNodes))];  
dbc = [dbc; bottomNodes,     3*ones(size(bottomNodes)),     zeros(size(bottomNodes))];
dbc = [dbc; bottomNodes,     4*ones(size(bottomNodes)),     zeros(size(bottomNodes))];  
dbc = [dbc; bottomNodes,     5*ones(size(bottomNodes)),     zeros(size(bottomNodes))];
dbc = [dbc; bottomNodes,     6*ones(size(bottomNodes)),     zeros(size(bottomNodes))];  

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
gp_z = 2;                  % number of integration points in z-direction
[gp, wgt] = gaussQuadrature(gp_x, gp_y, gp_z);   % calculate integration points and its weights
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
    pu = mesh.elDegree(e,1);
    pv = mesh.elDegree(e,2);
    Ce = mesh.elExtOpe{e,1};
    we = mesh.coords(sctr,4); % Tspline control points' weights
    Ke = zeros(nnElem,nnElem);
    for ipt = 1:size(gp,1)   % loop over integration points
        pt = gp(ipt,1:2);      % reference parametric coordinates for each integration point
        wt = wgt(ipt);       % weigths for each integration point
        zeta    = gp(ipt,3);
        gauPts  = parameterGaussMapping( elDoma, pt );   % gauss integration mapping
        j1      = jacobianGaussMapping( elDoma );            % jacobian value for gauss mapping   
        [R,dR,dR2] = computeTsplineBasis2ndDers([pu,pv],gauPts,Ce,we);
        dxdxi = dR*elCpts(:,1:3); 
        n_top = cross(dxdxi(1,:),dxdxi(2,:));
        n_bot = norm(n_top);
        e3 = n_top/n_bot;   % normal vector 
        ea = dxdxi(1,:) + dxdxi(2,:);
        ea = ea/norm(ea);
        eb = cross(e3,ea);
        eb = eb/norm(eb);
        e1 = sqrt(2)/2 * (ea - eb);  % local basis vector
        e2 = sqrt(2)/2 * (ea + eb);
        q = [e1; e2; e3]';
        Q = localGlobalVoigt( q );
        
        dxdxi2 = dR2*elCpts(:,1:3); 
        dN     = computeNormalDers(dxdxi(1,:), dxdxi(2,:), dxdxi2(1,:), dxdxi2(2,:), dxdxi2(3,:)); % normal derivatives

        jmatrix(1:2,1:3) = dxdxi + h /2 * zeta * dN;
        jmatrix(3,1:3)   = e3 * h /2;   % dxdzeta
        j2        = det(jmatrix);
        j2_inv    = inv(jmatrix); 
        dzetadx   = j2_inv(1:3,3);
        dR(3,:)   = 0;
        dN(3,:)   = 0;
        dRdx      =  jmatrix \ dR;
        dndx      =  jmatrix \ dN;
        
        Ri  = zeta *h/2 * R;
        Rix =  (dzetadx(1)*R + zeta*dRdx(1,:))*h/2;
        Riy =  (dzetadx(2)*R + zeta*dRdx(2,:))*h/2;
        Riz =  (dzetadx(3)*R + zeta*dRdx(3,:))*h/2;      
        
        n = e3;
        B = zeros(6,nnElem);    
        B(1,1:6:nnElem) =    dRdx(1,:);
        B(1,5:6:nnElem) = Rix * n(3) + Ri * dndx(1,3);
        B(1,6:6:nnElem) = -( Rix * n(2) + Ri * dndx(1,2) );
        
        B(2,2:6:nnElem) =    dRdx(2,:);
        B(2,4:6:nnElem) = -( Riy * n(3) + Ri * dndx(2,3) );
        B(2,6:6:nnElem) =    Riy * n(1) + Ri * dndx(2,1);
        
        B(3,3:6:nnElem) =    dRdx(3,:);
        B(3,4:6:nnElem) =    Riz * n(2) + Ri * dndx(3,2);
        B(3,5:6:nnElem) =  -(Riz * n(1) + Ri * dndx(3,1));
        
        B(4,1:6:nnElem) =    dRdx(2,:);
        B(4,2:6:nnElem) =    dRdx(1,:);
        B(4,4:6:nnElem) =  -(Rix * n(3) + Ri * dndx(1,3));
        B(4,5:6:nnElem) =    Riy * n(3) + Ri * dndx(2,3);
        B(4,6:6:nnElem) =    Rix * n(1) + Ri * dndx(1,1) - (Riy * n(2) + Ri * dndx(2,2));            
        
        B(5,2:6:nnElem) =   dRdx(3,:);
        B(5,3:6:nnElem) =   dRdx(2,:);
        B(5,4:6:nnElem) =   Riy * n(2) + Ri * dndx(2,2) - (Riz * n(3) + Ri * dndx(3,3));   
        B(5,5:6:nnElem) = -(Riy * n(1) + Ri * dndx(2,1));
        B(5,6:6:nnElem) =   Riz * n(1) + Ri * dndx(3,1);          
        
        B(6,1:6:nnElem) =   dRdx(3,:);
        B(6,3:6:nnElem) =   dRdx(1,:);
        B(6,4:6:nnElem) =   Rix * n(2) + Ri * dndx(1,2);  
        B(6,5:6:nnElem) =   Riz * n(3) + Ri * dndx(3,3) - (Rix * n(1) + Ri * dndx(1,1));   
        B(6,6:6:nnElem) = -(Riz * n(2) + Ri * dndx(3,2));

        Dg  = Q * Dl * Q';
        fac = j1 *j2 * wt;   
        Ke  = Ke + B' * Dg * B * fac;
    end   
    K(sctrB,sctrB) = K(sctrB,sctrB) +  Ke;
end
toc;

%****************************************************************%
%                   Enforce force loading                        %
%****************************************************************%
% find the elements for enforcing load
loadElem = [];
for e = 1:mesh.nElems                 % loop over elements
    sctr   = mesh.elNodeCnt{e,:};     % element control points index
    elCpts = mesh.coords(sctr,:);     % coordinates of element control points
    Ce     = mesh.elExtOpe{e,1};
    [Pbe, wbe] = computeBezierCtrlPts(elCpts(:,1:3), elCpts(:,4), Ce);
    flg    = 1;
    if all( Pbe(:,3) > 90.8 )
        for i = 1:size(Pbe,1)
            r = sqrt(Pbe(i,1)^2 + Pbe(i,2)^2);
            if r > 144
                flg = -1;
                break;
            end
        end
        if 1 == flg
            loadElem = [loadElem; e];
        end
    end
end

fprintf("Enforce force load \n");
tic
gp_x = max(mesh.elDegree(:,1))+1;           % number of integration points in x-direction
gp_y = max(mesh.elDegree(:,2))+1;           % number of integration points in y-direction
[gp, wgt] = gaussQuadrature(gp_x, gp_y);   % calculate integration points and its weights
elDoma = [0,1,0,1];
g      = 1;
for i = 1:length(loadElem)                 % loop over elements
    e      = loadElem(i);
    sctr   = mesh.elNodeCnt{e,:};     % element control points index
    elCpts = mesh.coords(sctr,:);     % coordinates of element control points
    strf   = sctr*6-3;
    pu = mesh.elDegree(e,1);
    pv = mesh.elDegree(e,2);
    Ce = mesh.elExtOpe{e,1};
    we = mesh.coords(sctr,4); % Tspline control points' weights
    for ipt = 1:size(gp,1)        % loop over integration points
        pt       = gp(ipt,:);      % reference parametric coordinates for each integration point
        wt       = wgt(ipt);       % weigths for each integration point
        gauPts   = parameterGaussMapping( elDoma, pt );       % gauss integration mapping
        j1       = jacobianGaussMapping( elDoma );            % jacobian value for gauss mapping   
        [R,dR] = computeTsplineBasisDers([pu,pv],gauPts,Ce,we);
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
msize = getMeshSize(mesh.coords);
len = norm([msize(1), msize(3), msize(5)] - [msize(2), msize(4), msize(6)]);
% factor = len/max(abs(u(3:6:end)))/10;  
factor = 1.0;
% post-processing
num1 = 5;  num2 = 5;  
vmesh = buildVisualMesh(num1,num2,mesh);
count = 1;
for e = 1:mesh.nElems
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
    pu = mesh.elDegree(e,1);
    pv = mesh.elDegree(e,2);
    Ce = mesh.elExtOpe{e,1};
    we = mesh.coords(sctr,4); % Tspline control points' weights
    Ke = zeros(nnElem,nnElem);
    for j = 1:(num1+1)*(num2+1)
        pt = vmesh.tripts(j,:);  
        zeta    = 0;
        [R,~,~] = computeTsplineBasis2ndDers([pu,pv],pt,Ce,we);
        edsp   = u(sctrB);
        edsp   = reshape(edsp, dof, nn);
        vmesh.displacement(count,:) = R * edsp';
        vmesh.node(count,1:3) = R*elCpts(:,1:3) + vmesh.displacement(count,1:3)*factor;
        count = count+1;
    end       
end

figure(2)
title( 'Vertical displacement', 'FontSize', 14', 'FontWeight', 'Bold') ;
p = patch('Faces',vmesh.element, 'Vertices', vmesh.node); 
cdata = vmesh.displacement(:,3);
set(p,'FaceColor','interp','FaceVertexCData',cdata);
view(37,37);
set(p,'EdgeColor','none');
hold on;
for i = 1:mesh.nElems
    plot3(vmesh.node(vmesh.linmesh{i,1},1),vmesh.node(vmesh.linmesh{i,1},2),vmesh.node(vmesh.linmesh{i,1},3),'k-');
end
axis equal;
colorbar;
mycolor = abaqusColorMap(12);
colormap(mycolor);
% colormap jet;
hold on;
vmesh2 = getVisualMesh(mesh);
for i = 1:mesh.nElems
    plot3(vmesh2.node(vmesh2.linmesh{i,1},1),vmesh2.node(vmesh2.linmesh{i,1},2),vmesh2.node(vmesh2.linmesh{i,1},3),'Color',[0.5,0.5,0.5]);
end
xlabel('x'); ylabel('y'); zlabel('z');
toc;








