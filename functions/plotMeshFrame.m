function plotMeshFrame(mesh,u)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%
% Plot iga mesh frame withour color rendering
%  Input:
%    mesh - iga mesh structure, using the function 'build_tiga_mesh'
%    u    - displacement results, [n x 3] or [n x 1], it is optional
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%  - 2022
%
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
if nargin == 1
    vmesh = getVisualMesh(mesh); % discretize the iga mesh
elseif nargin == 2
    vmesh = getVisualMesh(mesh,u); % discretize the iga mesh
    vmesh.node = vmesh.node + vmesh.u(:,1:3);
end
axis equal;
for i = 1:mesh.nElems
    hold on;
    plot3(vmesh.node(vmesh.linmesh{i,1},1),vmesh.node(vmesh.linmesh{i,1},2),vmesh.node(vmesh.linmesh{i,1},3),'-','Color',[102,102,102]./255);
end
view(3);
end

