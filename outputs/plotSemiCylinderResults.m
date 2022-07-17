%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%  Post-Processing: final configuration
%  Large deformation of a pinched semi-cylindrical shell with a clamped edge
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
filename = 'geoSemicylinder.iga';
geo  = readIgaFile(filename); 
mesh = buildIgaMesh( geo );
figure(1)
plotIgaMesh(mesh,1);
view(30,30);
h    = 0.03;
nu   = 0.3;             % Poisson ratio
E    = 2.0685e7;
mIndex = 40;
if 40 == mIndex       % linear St. Venant-Kirchhoff type material.
    mat    = [mIndex,E,nu,h];     % mat = [index, E, nu, h],  
    filename = 'postSemiCylinderM40E1664.msh';
elseif 41 == mIndex   % incompressible Neo-Hookean material
    A10    = E/2/(1+nu)/2;
    mat    = [mIndex,0, A10, h];  % mat = [index, K, A10, h],
    filename = 'postSemiCylinderM41E1664.msh';
elseif 42 == mIndex   % compressible Neo-Hookean material
    lamda  = E*nu/(1+nu)/(1-2*nu); 
    mu     = E/2/(1+nu);
    mat    = [mIndex,mu,lamda,h];  % mat = [index, mu, lamda, h], 
    filename = 'postSemiCylinderM42E1664.msh';
end

[tsteps,usteps] = readMshFile(filename);

% gg = [5:5:30,32];
gg = 32;
for k = 1:length(gg)
    figure(k+1)
    u = usteps{1,gg(k)};
    mesh2 = mesh;
    mesh2.coords(:,1:3) = mesh2.coords(:,1:3) + u;
    vmesh = getVisualMesh(mesh2,u,mat);
    vMesh{1,1} = vmesh;
    vMesh{1,2} = vmesh;
    vMesh{1,2}.node(:,1) = - vMesh{1,2}.node(:,1);
    vMesh{1,2}.u(:,1)    = - vMesh{1,2}.u(:,1);
    
    for i = 1:2
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
    view(30,30);
    axis equal;
    axis off;
    mycolor = abaqusColorMap(24);
    colormap(mycolor);
    C  = colorbar;
    caxis manual
    colorMin = C.Limits(1);
    colorMax = C.Limits(2);
    caxis([colorMin colorMax]);
    set(C,'Ticks',linspace(colorMin,colorMax,13));
    set(get(C,'title'),'string','Um');
end

gg = 0;
for k = 1:length(gg)
    u = zeros(size(usteps{1,end}));
    mesh2 = mesh;
    mesh2.coords(:,1:3) = mesh2.coords(:,1:3) + u;
    vmesh = getVisualMesh(mesh2,u,mat);
    vMesh{1,1} = vmesh;
    vMesh{1,2} = vmesh;
    vMesh{1,2}.node(:,1) = - vMesh{1,2}.node(:,1);
    vMesh{1,2}.u(:,1)    = - vMesh{1,2}.u(:,1);
    for i = 1:2
        hold on;
        for j = 1:mesh.nElems
            hold on;
            plot3(vMesh{1,i}.node(vMesh{1,i}.linmesh{j,1},1),...
                  vMesh{1,i}.node(vMesh{1,i}.linmesh{j,1},2),...
                  vMesh{1,i}.node(vMesh{1,i}.linmesh{j,1},3),'Color',[0.5,0.5,0.5]);
        end
    end
end









