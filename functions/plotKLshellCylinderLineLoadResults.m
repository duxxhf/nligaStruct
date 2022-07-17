function plotKLshellCylinderLineLoadResults(mesh,u,L)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%
% Plot displacement results of large deformation of the pinched cylinder
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
mesh.coords(:,1:3) = mesh.coords(:,1:3) + u;
vmesh = getVisualMesh(mesh,u);
vMesh{1,1} = vmesh;
vMesh{1,2} = vmesh;
vMesh{1,2}.node(:,2) = L - vMesh{1,2}.node(:,2);
vMesh{1,3} = vmesh; 
vMesh{1,3}.node(:,1) = - vMesh{1,3}.node(:,1); 
vMesh{1,4} = vMesh{1,2}; 
vMesh{1,4}.node(:,1) = - vMesh{1,4}.node(:,1); 

cdata = sqrt(vmesh.u(:,1).^2+vmesh.u(:,2).^2+vmesh.u(:,3).^2);
for i = 1:4
    hold on;
    p = patch('Faces',vMesh{1,i}.element, 'Vertices', vMesh{1,i}.node);         % plot the surface
    set(p,'FaceColor','interp','FaceVertexCData',cdata,'EdgeColor','none');
    for j = 1:mesh.nElems
        hold on;
        plot3(vMesh{1,i}.node(vMesh{1,i}.linmesh{j,1},1),...
              vMesh{1,i}.node(vMesh{1,i}.linmesh{j,1},2),...
              vMesh{1,i}.node(vMesh{1,i}.linmesh{j,1},3),'k-');
    end
end

axis equal;
h = colorbar;
t = get(h,'Limits');
set(h,'Ticks',linspace(t(1),t(2),7));
set(get(h,'title'),'string','Um');
mycolor = abaqusColorMap(12);
colormap(mycolor);
view(3);
end

