function plotIgaMesh(mesh,type)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%
% Plot initial iga mesh model or obtained results
%  Input:
%    mesh - iga mesh structure, using the function 'build_tiga_mesh'
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%  - 2022
%
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
vmesh = getVisualMesh(mesh); % discretize the iga mesh
if 1 == type
    p = patch('Faces',vmesh.element, 'Vertices', vmesh.node); % plot the surface
    set(p,'FaceColor','y','EdgeColor','none');
    set(gcf,'renderer','opengl');
    light;
    axis equal;
    for i = 1:mesh.nElems
        hold on;
        plot3(vmesh.node(vmesh.linmesh{i,1},1),vmesh.node(vmesh.linmesh{i,1},2),vmesh.node(vmesh.linmesh{i,1},3),'k-');
    end
    view(3);
elseif 2 == type
    p = patch('Faces',vmesh.element, 'Vertices', vmesh.node); % plot the surface
    set(p,'FaceColor','y','EdgeColor','none');
    set(gcf,'renderer','opengl');
    light;
    axis equal;
elseif 3 == type
    for i = 1:mesh.nElems
        hold on;
        plot3(vmesh.node(vmesh.linmesh{i,1},1),vmesh.node(vmesh.linmesh{i,1},2),vmesh.node(vmesh.linmesh{i,1},3),'-','Color',[102,102,102]./255);
    end
    axis equal;
end

end

