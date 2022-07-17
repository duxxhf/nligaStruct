clc;
clear;

filename = 'geoCylinderLineLoad2.iga';
geo  = readIgaFile(filename); 
mesh = buildIgaMesh( geo );
R = 0.09; L = 0.3;
% scale the model into the range 
mesh.coords(:,1:3) = 1/100*mesh.coords(:,1:3);

filename = 'postCylinderLineLoadM42E810.msh';
[tsteps,usteps] = readMshFile(filename);

u = usteps{1,end};
vmesh = getVisualMesh(mesh,u);
cdata = sqrt(vmesh.u(:,1).^2+vmesh.u(:,2).^2+vmesh.u(:,3).^2);
colorMin = min(cdata);
colorMax = max(cdata);
figure(1)
mesh2 = mesh;
mesh2.coords(:,1:3) = mesh2.coords(:,1:3) + u;
vmesh = getVisualMesh(mesh2,u);
vMesh{1,1} = vmesh;
vMesh{1,2} = vmesh;
vMesh{1,2}.node(:,2) = L - vMesh{1,2}.node(:,2);
vMesh{1,3} = vmesh; 
vMesh{1,3}.node(:,1) = - vMesh{1,3}.node(:,1); 
vMesh{1,4} = vMesh{1,2}; 
vMesh{1,4}.node(:,1) = - vMesh{1,4}.node(:,1); 
for i = 1:4
    hold on;
    p = patch('Faces',vMesh{1,i}.element, 'Vertices', vMesh{1,i}.node);         % plot the surface
    cdata = sqrt(vMesh{1,i}.u(:,1).^2+vMesh{1,i}.u(:,2).^2+vMesh{1,i}.u(:,3).^2);
    set(p,'FaceColor','interp','FaceVertexCData',cdata,'EdgeColor','none');
    for j = 1:mesh.nElems
        hold on;
        plot3(vMesh{1,i}.node(vMesh{1,i}.linmesh{j,1},1),...
              vMesh{1,i}.node(vMesh{1,i}.linmesh{j,1},2),...
              vMesh{1,i}.node(vMesh{1,i}.linmesh{j,1},3),'k-');
    end
end
view(40,35);
caxis manual
caxis([colorMin colorMax]);
C  = colorbar;
set(C,'Ticks',linspace(colorMin,colorMax,13));
mycolor = abaqusColorMap(12);
colormap(mycolor);
axis equal;
axis off;

figure(2);
steps = 1:3:18;
for gg = 1:length(steps)
    u = usteps{1,steps(gg)};
    mesh2 = mesh;
    mesh2.coords(:,1:3) = mesh2.coords(:,1:3) + u;
    vmesh = getVisualMesh(mesh2,u);
    vMesh{1,1} = vmesh;
    vMesh{1,2} = vmesh;
    vMesh{1,2}.node(:,2) = L - vMesh{1,2}.node(:,2);
    vMesh{1,3} = vmesh; 
    vMesh{1,3}.node(:,1) = - vMesh{1,3}.node(:,1); 
    vMesh{1,4} = vMesh{1,2}; 
    vMesh{1,4}.node(:,1) = - vMesh{1,4}.node(:,1); 
    h  = subplot(2,3,gg);
    for i = 1:4
        hold on;
        p = patch('Faces',vMesh{1,i}.element, 'Vertices', vMesh{1,i}.node);         % plot the surface
        cdata = sqrt(vMesh{1,i}.u(:,1).^2+vMesh{1,i}.u(:,2).^2+vMesh{1,i}.u(:,3).^2);
        set(p,'FaceColor','interp','FaceVertexCData',cdata,'EdgeColor','none');
        for j = 1:mesh.nElems
            hold on;
            plot3(vMesh{1,i}.node(vMesh{1,i}.linmesh{j,1},1),...
                  vMesh{1,i}.node(vMesh{1,i}.linmesh{j,1},2),...
                  vMesh{1,i}.node(vMesh{1,i}.linmesh{j,1},3),'k-');
        end
    end
    view(40,35);
    axis([-0.12,0.12,0,0.3,-0.09,0.09]);
    caxis manual
    caxis([colorMin colorMax]);
    box on;
    axis equal;
    title(['P/P_{max}=',num2str(tsteps(steps(gg)))],'Interpreter','tex') ;
end

tt = get(subplot(2,3,6),'Position');
C  = colorbar;
set(C,'Ticks',linspace(colorMin,colorMax,13));
mycolor = abaqusColorMap(12);
colormap(mycolor);
C.Position = [tt(1)+tt(3)+0.02, tt(2), 0.01, 0.8];
title(C,'U magnitude');


