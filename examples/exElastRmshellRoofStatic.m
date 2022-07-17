%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%
%  An example for structural analysis of the benchmark problem 
%  - Scordelis¨CLo roof using RM shell
%
%  ---------------------------------------
%
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%  - 7-DEC-2021
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
    filename = 'geoQuarterRoofR3.iga';
    geo  = readIgaFile(filename);
    mesh = buildIgaMesh( geo );
elseif strcmp(meshType, 'NURBS')
    R = 25; L = 50; theta = 40; p = 4; q = p; elemX = 20; elemY = elemX; 
    geo  = geoQuarterRoofShellNurbs(L,R,theta,p,q,elemX,elemY);
    mesh = buildIgaMesh( geo );
end

% plot scaled model
figure(1)
plotIgaMesh(mesh,1);
view(50,30);


%****************************************************************%
%               Define material properties and dof               %
%****************************************************************%
fprintf("Define material properties and dof \n");
tic
% material properties
nu   = 0.0;             % Poisson ratio
E    = 4.32e8;          % Youngs modulus 
g    = 90;              % Gravity loading
h    = 0.25;            % thickness   
kapa = 5/6;             % Shear correction factor
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
bottomNodes   = find( abs(mesh.coords(:,2)) < tol );
topNodes      = find( abs(mesh.coords(:,2)-L/2) < tol );
leftNodes     = find( abs(mesh.coords(:,1)) < tol );
% sort topNodes
topNodes = sortrows([floor(mesh.coords(topNodes,1)*1000)/1000, topNodes],1,'ascend');
topNodes = topNodes(:,2);
nodeA = topNodes(end);
% plot constrained control points
hold on;
plot3(mesh.coords(bottomNodes,1),mesh.coords(bottomNodes,2),mesh.coords(bottomNodes,3),'co');
plot3(mesh.coords(topNodes,1),mesh.coords(topNodes,2),mesh.coords(topNodes,3),'ro');
plot3(mesh.coords(leftNodes,1),mesh.coords(leftNodes,2),mesh.coords(leftNodes,3),'mo');

% impose displacement boundary conditions
dbc = [];    % dbc = [node index, direction, prescribed displacement]
% rigid diaphram boundary, ux = uz = theta_y = 0
dbc = [dbc; bottomNodes,     1*ones(size(bottomNodes)),    zeros(size(bottomNodes))];  
dbc = [dbc; bottomNodes,     3*ones(size(bottomNodes)),    zeros(size(bottomNodes))];
dbc = [dbc; bottomNodes,     5*ones(size(bottomNodes)),    zeros(size(bottomNodes))];
% top symmetric boundary, uy = theta_x = theta_z = 0
dbc = [dbc; topNodes,        2*ones(size(topNodes)),       zeros(size(topNodes))];      
dbc = [dbc; topNodes,        4*ones(size(topNodes)),       zeros(size(topNodes))];
dbc = [dbc; topNodes,        6*ones(size(topNodes)),       zeros(size(topNodes))];
% left symmetric boundary, ux = theta_y = theta_z = 0
dbc = [dbc; leftNodes,       1*ones(size(leftNodes)),      zeros(size(leftNodes))];      
dbc = [dbc; leftNodes,       5*ones(size(leftNodes)),      zeros(size(leftNodes))];
dbc = [dbc; leftNodes,       6*ones(size(leftNodes)),      zeros(size(leftNodes))];
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
gp_z = 3;                  % number of integration points in z-direction
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
        gauPts = parameterGaussMapping( elDoma, pt );   % gauss integration mapping
        j1 = jacobianGaussMapping( elDoma );            % jacobian value for gauss mapping   
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
%                   Enforce gravity loading                      %
%****************************************************************%
fprintf("Enforce gravity loading \n");
tic
gp_x = max(mesh.elDegree(:,1))+1;           % number of integration points in x-direction
gp_y = max(mesh.elDegree(:,2))+1;           % number of integration points in y-direction
[gp, wgt] = gaussQuadrature(gp_x, gp_y);   % calculate integration points and its weights
elDoma = [0,1,0,1];

for e = 1:mesh.nElems                 % loop over elements
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
        [R,dR]   = computeTsplineBasisDers([pu,pv],gauPts,Ce,we);
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
factor = len/max(abs(u(3:6:end)))/10;  
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
        [R,dR,dR2] = computeTsplineBasis2ndDers([pu,pv],pt,Ce,we);
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
        strain = B * u(sctrB);
        stress = Dg * strain;
        edsp   = u(sctrB);
        edsp   = reshape(edsp, dof, nn);
        vmesh.displacement(count,:) = R * edsp';
        vmesh.node(count,1:3) = R*elCpts(:,1:3) + vmesh.displacement(count,1:3)*factor;
        vmesh.stress(count,:) = stress';
        count = count+1;
    end       
end
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
    cdata = vMesh{1,i}.displacement(:,3);
    set(p,'FaceColor','interp','FaceVertexCData',cdata,'EdgeColor','none');      
    for j = 1:mesh.nElems
        hold on;
        plot3(vMesh{1,i}.node(vMesh{1,i}.linmesh{j,1},1),...
              vMesh{1,i}.node(vMesh{1,i}.linmesh{j,1},2),...
              vMesh{1,i}.node(vMesh{1,i}.linmesh{j,1},3),'k-');
    end
end

axis equal;
C = colorbar;
caxis manual
colorMin = C.Limits(1);
colorMax = C.Limits(2);
caxis([colorMin colorMax]);
set(C,'Ticks',linspace(colorMin,colorMax,7));
mycolor = abaqusColorMap(24);
colormap(mycolor);
view(37,37);
toc;

fprintf('The maximum vertical displacement is %12.6e.\n',-min(u(nodeA*6-3)));

