%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%  Post-Processing: deformed edges
%  Large deformation of a pinched cylindrical shell with rigid diaphragms 
%  subjected to point loads
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

filename = 'geoPichchedCylinder.iga';
geo  = readIgaFile(filename); 
mesh = buildIgaMesh( geo );
% scale the model into the range [0,100]x[0,100]x[0,100]
L = 200;
mesh.coords(:,1:3) = L/60*mesh.coords(:,1:3);
% plot scaled model
figure(1)
plotIgaMesh(mesh,3);
view(110,30);

filename = 'postPinchedCylinderM40E2860.msh';
[tsteps,usteps] = readMshFile(filename);

u = usteps{1,end};
vmesh = getVisualMesh(mesh,u);
cdata = sqrt(vmesh.u(:,1).^2+vmesh.u(:,2).^2+vmesh.u(:,3).^2);
colorMin = min(cdata);
colorMax = max(cdata);

gg = 103;
for k = 1:length(gg)
    figure(k+1)
    u = usteps{1,gg(k)};
    mesh2 = mesh;
    mesh2.coords(:,1:3) = mesh2.coords(:,1:3) + u;
    vmesh = getVisualMesh(mesh2,u);
    vMesh{1,1} = vmesh;
    vMesh{1,2} = vmesh;
    vMesh{1,2}.node(:,2) = L - vMesh{1,2}.node(:,2);
    vMesh{1,3} = vmesh;
    vMesh{1,3}.node(:,3) = - vMesh{1,3}.node(:,3);
    vMesh{1,4} = vMesh{1,2}; 
    vMesh{1,4}.node(:,3) = - vMesh{1,4}.node(:,3); 
    vMesh{1,5} = vmesh; 
    vMesh{1,5}.node(:,1) = - vMesh{1,5}.node(:,1); 
    vMesh{1,6} = vMesh{1,2}; 
    vMesh{1,6}.node(:,1) = - vMesh{1,6}.node(:,1); 
    vMesh{1,7} = vMesh{1,5}; 
    vMesh{1,7}.node(:,3) = - vMesh{1,7}.node(:,3); 
    vMesh{1,8} = vMesh{1,6}; 
    vMesh{1,8}.node(:,3) = - vMesh{1,8}.node(:,3); 
    for i = 1:8
        hold on;
        p = patch('Faces',vMesh{1,i}.element, 'Vertices', vMesh{1,i}.node);         % plot the surface
        cdata = sqrt(vMesh{1,i}.u(:,1).^2+vMesh{1,i}.u(:,2).^2+vMesh{1,i}.u(:,3).^2);
        set(p,'FaceColor','interp','FaceVertexCData',cdata,'EdgeColor','none');      
%         for j = 1:mesh.nElems
%             hold on;
%             plot3(vMesh{1,i}.node(vMesh{1,i}.linmesh{j,1},1),...
%                   vMesh{1,i}.node(vMesh{1,i}.linmesh{j,1},2),...
%                   vMesh{1,i}.node(vMesh{1,i}.linmesh{j,1},3),'k-');
%         end
    end
    view(75,25);
    axis equal;
    axis on;
    mycolor = abaqusColorMap(24);
    colormap(mycolor);
    set(gcf,'renderer','opengl');
    light;
    C  = colorbar;
    set(C,'Ticks',linspace(colorMin,colorMax,7));
    set(get(C,'title'),'string','Um');
    caxis manual
    caxis([colorMin colorMax]);
end











