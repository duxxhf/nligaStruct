%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%  Post-Processing: final configuration
%  Large deformation of a pinched hemispherical shell subjected to point loads
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

filename = 'geoHemisphere.iga';
geo  = readIgaFile(filename); 
mesh = buildIgaMesh( geo );

filename = 'postHemisphereHoleM40E500.msh';
[tsteps,usteps] = readMshFile(filename);

u = usteps{1,end};
vmesh = getVisualMesh(mesh,u);
cdata = sqrt(vmesh.u(:,1).^2+vmesh.u(:,2).^2+vmesh.u(:,3).^2);
colorMin = min(cdata);
colorMax = max(cdata);

gg = 36;
for k = 1:length(gg)
    figure(k)
    u = usteps{1,gg(k)};
    mesh2 = mesh;
    mesh2.coords(:,1:3) = mesh2.coords(:,1:3) + u;
    vmesh = getVisualMesh(mesh2,u);
    vMesh{1,1} = vmesh;
    vMesh{1,2} = vmesh;
    vMesh{1,2}.node(:,2) = - vMesh{1,2}.node(:,2);
    vMesh{1,3} = vmesh;
    vMesh{1,3}.node(:,1) = - vMesh{1,3}.node(:,1);
    vMesh{1,4} = vMesh{1,2}; 
    vMesh{1,4}.node(:,1) = - vMesh{1,4}.node(:,1); 
    cdata = sqrt(vmesh.u(:,1).^2+vmesh.u(:,2).^2+vmesh.u(:,3).^2);
    for i = 1:4
        hold on;
        p = patch('Faces',vMesh{1,i}.element, 'Vertices', vMesh{1,i}.node);         % plot the surface
        set(p,'FaceColor','interp','FaceVertexCData',cdata,'EdgeColor','none');
%         set(p,'FaceColor','w','EdgeColor','none');
        for j = 1:mesh.nElems
            hold on;
            plot3(vMesh{1,i}.node(vMesh{1,i}.linmesh{j,1},1),...
                  vMesh{1,i}.node(vMesh{1,i}.linmesh{j,1},2),...
                  vMesh{1,i}.node(vMesh{1,i}.linmesh{j,1},3),'k-');
        end
    end
    view(3);
    caxis manual
    caxis([colorMin colorMax]);
    axis equal;
    axis off;
    mycolor = abaqusColorMap(12);
    colormap(mycolor);
    C  = colorbar;
    set(C,'Ticks',linspace(colorMin,colorMax,13));
    set(get(C,'title'),'string','Um');
end


