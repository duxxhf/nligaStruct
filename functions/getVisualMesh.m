function vmesh = getVisualMesh(mesh,u,mat)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%
% Build a discrete mesh for visualization
%  Input:
%    mesh - iga mesh structure, using the function 'build_tiga_mesh'
%    u    - displacements
%    mat  - material paramters
%  Output:
%    vmesh - visualization mesh
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%  - 2022
%
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

num1 = 5;  num2 = 5;  
trinum = mesh.nElems*num1*num2*2;        % total elements
trimesh = zeros(trinum,3);   % mesh connectivity
linmesh = cell(mesh.nElems,1);
count = 1;
for k = 1:mesh.nElems
    for j = 1:num2
        for i = 1:num1
            p1 = (k-1)*(num1+1)*(num2+1)+(j-1)*(num1+1)+i;
            p2 = (k-1)*(num1+1)*(num2+1)+(j-1)*(num1+1)+i+1;
            p3 = (k-1)*(num1+1)*(num2+1)+(j-1)*(num1+1)+i+1+(num1+1);
            p4 = (k-1)*(num1+1)*(num2+1)+(j-1)*(num1+1)+i+(num1+1);
            trimesh(count,:) = [p1 p2 p4];
            trimesh(count+1,:) = [p2 p3 p4];
            count = count+2;
        end
    end
    linmesh{k,1} = zeros(num1*2+num2*2+1,1);
    linmesh{k,1}(1:num1) = (1:num1);
    linmesh{k,1}((num1+1):(num1+num2)) =  (num1+1):num1+1:(num1+1)*num2;
    linmesh{k,1}((num1+num2+1):(num1*2+num2)) = (num1+1)*(num2+1):-1:(num1+1)*num2+2;
    linmesh{k,1}((num1*2+num2+1):num1*2+num2*2) = ((num1+1)*num2+1) : -(num1+1) : (num1+2);
    linmesh{k,1} = linmesh{k,1} + (k-1)*(num1+1)*(num2+1);
    linmesh{k,1}(end) = linmesh{k,1}(1);
end

offset = 0;
uu = linspace(0+offset,1-offset,num1+1);
vv = linspace(0+offset,1-offset,num2+1);
count = 1;
tripts = zeros((num1+1)*(num2+1),2);

for j = 1:num2+1
    for i = 1:num1+1
        tripts(count,:) = [uu(i), vv(j)];   % parametric nodal points
        count = count+1;
    end
end
vmesh.node = zeros(mesh.nElems*(num1+1)*(num2+1),3);
vmesh.element = trimesh;
if nargin == 2
    vmesh.u = zeros(mesh.nElems*(num1+1)*(num2+1),size(u,2));
elseif nargin == 3
    vmesh.u = zeros(mesh.nElems*(num1+1)*(num2+1),3);
    vmesh.n = zeros(mesh.nElems*(num1+1)*(num2+1),3);  % normal forces
    vmesh.m = zeros(mesh.nElems*(num1+1)*(num2+1),3);  % bending moments
    vmesh.vms = zeros(mesh.nElems*(num1+1)*(num2+1),1);  % mises stress
    vmesh.s = zeros(mesh.nElems*(num1+1)*(num2+1),3);  % mises stress
end
count = 1;
for e = 1:mesh.nElems
    sctr   = mesh.elNodeCnt{e,1};       % element control points index
    elCpts = mesh.coords(sctr,1:3);     % coordinates of element control points
    if nargin == 2 || nargin == 3
        eU = u(sctr,:);
    end
    nn     = numel(sctr);               % number of control points for each element
    nnElem = nn*mesh.dim;               % dof for each element
    sctrB  = zeros(1, nnElem); 
    for i  = 1:mesh.dim
        sctrB(i:mesh.dim:nnElem) = mesh.dim*(sctr-1) + i;  % displacement in i-th direction
    end
    pu     = mesh.elDegree(e,1);  % degree u
    pv     = mesh.elDegree(e,2);  % degree v
    we     = mesh.coords(sctr,4); % Tspline control points' weights
    Ce     = mesh.elExtOpe{e,1};  % elemental extration operator
    for j = 1:(num1+1)*(num2+1)
        pt  = tripts(j,:);  
        [R,dR,dR2] = computeTsplineBasis2ndDers([pu,pv],pt,Ce,we);
        xe = R * elCpts; 
        vmesh.node(count,:) = xe; 
        if nargin == 2 
            vmesh.u(count,:) = R*eU;
        elseif nargin == 3
            vmesh.u(count,:) = R*eU;
            % TBA
        end
        count = count+1;
    end    
end

vmesh.linmesh = linmesh;


end

