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
meshType = 'T-SPLINE';  % or 'NURBS'
% meshType = 'NURBS';  % or 'NURBS'
R = 25; L = 50;
if strcmp(meshType, 'T-SPLINE')
    filename = 'geoRoof.iga';
    geo  = readIgaFile(filename);
    mesh = buildIgaMesh( geo );
    mesh.coords(:,1:3) = 5/6*mesh.coords(:,1:3) ;
elseif strcmp(meshType, 'NURBS')
    theta = 40; p = 3; q = p; elemX = 20; elemY = elemX; 
    geo  = geoRoofShellNurbs(L,R,theta,p,q,elemX,elemY);
    mesh = buildIgaMesh( geo );
end

% plot scaled model
figure(1)
plotIgaMesh(mesh,1);
view(50,30);
toc


%****************************************************************%
%               Define material properties and dof               %
%****************************************************************%
fprintf("Define material properties and dof \n");
tic
% material properties
nu   = 0.0;             % Poisson ratio
E    = 4.32e8;          % Youngs modulus 
kapa = 5/6;             % Shear correction factor
h    = 0.25;            % thickness   
nModes = 10;            % number of mode shapes
rho  = 8000;            % density 
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
xmax = max(mesh.coords(:,1)); xmin = min(mesh.coords(:,1));
bottomNodes   = find( abs(mesh.coords(:,2)) < tol );
topNodes      = find( abs(mesh.coords(:,2)-L) < tol );
leftNodes     = find( abs(mesh.coords(:,1)-xmin) < tol );
rightNodes    = find( abs(mesh.coords(:,1)-xmax) < tol );
fixedNodes    = unique([bottomNodes;topNodes;leftNodes;rightNodes]);
fixed_dof     = zeros(length(fixedNodes)*6,1);
for i = 1:6
    fixed_dof(i:6:end) = (fixedNodes-1)*6+i;
end

% plot constrained control points
hold on;
plot3(mesh.coords(bottomNodes,1),mesh.coords(bottomNodes,2),mesh.coords(bottomNodes,3),'co');
plot3(mesh.coords(topNodes,1),mesh.coords(topNodes,2),mesh.coords(topNodes,3),'ro');
plot3(mesh.coords(leftNodes,1),mesh.coords(leftNodes,2),mesh.coords(leftNodes,3),'bd');
plot3(mesh.coords(rightNodes,1),mesh.coords(rightNodes,2),mesh.coords(rightNodes,3),'m*');
% 

%****************************************************************%
%           Initialize and assemble stiffness matrix             %
%****************************************************************%
fprintf("Initialize and assemble stiffness matrix \n");
tic
nDofs = dof * mesh.nCpts;    % total dofs
K = sparse(nDofs,nDofs);     % stiffness matrix 
M = sparse(nDofs,nDofs);     % mass matrix

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
    sctrw = 3*sctr-2;
    pu = mesh.elDegree(e,1);
    pv = mesh.elDegree(e,2);
    Ce = mesh.elExtOpe{e,1};
    we = mesh.coords(sctr,4); % Tspline control points' weights
    Ke = zeros(nnElem,nnElem);
    for ipt = 1:size(gp,1)   % loop over integration points
        pt = gp(ipt,1:2);      % reference parametric coordinates for each integration point
        wt = wgt(ipt);       % weigths for each integration point
        zeta   = gp(ipt,3);
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
        K(sctrB,sctrB)  = K(sctrB,sctrB) + B' * Dg * B * fac;
        N  = zeros(3, nnElem);
        N(1,1:6:end) = R; N(1,5:6:end) = Ri*n(3); N(1,6:6:end) = -Ri*n(2);
        N(2,2:6:end) = R; N(2,4:6:end) =-Ri*n(3); N(2,6:6:end) =  Ri*n(1);
        N(3,3:6:end) = R; N(3,4:6:end) = Ri*n(2); N(3,5:6:end) = -Ri*n(1);
        M(sctrB,sctrB) = M(sctrB,sctrB) + rho * (N' * N) * fac;
    end   
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
toc


%****************************************************************%
%                       Post-processing                          %
%****************************************************************%
fprintf("Post-processing \n");
tic
figure(2);clf;
set( gcf, 'Color', 'White', 'Unit', 'Normalized', 'Position', [0.1,0.1,0.6,0.6] );
title( 'Mode Shapes of Roof Shell Using IGA', 'FontSize', 14', 'FontWeight', 'Bold') ;
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
    scale = 0;   % scale factor
    mesh2 = mesh;
    mesh2.coords(:,1) = mesh2.coords(:,1) + U(1:6:end)*scale;
    mesh2.coords(:,2) = mesh2.coords(:,2) + U(2:6:end)*scale;
    mesh2.coords(:,3) = mesh2.coords(:,3) + U(3:6:end)*scale;
    vmesh = getVisualMesh(mesh2,[U(1:6:end),U(2:6:end),U(3:6:end)]);
    rowId = ceil( i / nCol ) ;
    colId = i - (rowId - 1) * nCol ;
    axes( 'Position', [colX(colId), rowY(rowId), colW, rowH] ) ;
    p = patch('Faces',vmesh.element, 'Vertices', vmesh.node);         % plot the surface
    axis equal
    axis off;
    cdata = sqrt(vmesh.u(:,1).^2+vmesh.u(:,2).^2+vmesh.u(:,3).^2);
    set(p,'FaceColor','interp','FaceVertexCData',cdata,'EdgeColor','none');
    hold on;
    for j = 1:mesh.nElems
        hold on;
        plot3(vmesh.node(vmesh.linmesh{j,1},1),vmesh.node(vmesh.linmesh{j,1},2),vmesh.node(vmesh.linmesh{j,1},3),'k-');
    end
    view(2);
    set(gca,'xtick',[],'ytick',[],'ztick',[]);
    title( [sprintf( 'Modes %d', i), ', \omega = ', sprintf('%3.3f',beta(index_beta(i)))]) ;  
    mycolor = abaqusColorMap(12);
    colormap(mycolor);
    hold off;
end
toc

