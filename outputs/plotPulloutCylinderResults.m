%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%  Post-Processing: final configuration
%  Pullout of cylindrical shell under large deformation
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

filename = 'geoPulloutCylinder1.iga';
geo  = readIgaFile(filename); 
mesh = buildIgaMesh( geo );
figure(1)
plotIgaMesh(mesh,3);

R = 4.953; L  = 10.35;
h    = 0.094;
nu   = 0.3125;             % Poisson ratio
E    = 10.5e6;
A10  = E/2/(1+nu)/2;
mIndex = 40;
if 40 == mIndex       % linear St. Venant-Kirchhoff type material.
    mat    = [mIndex,E,nu,h];     % mat = [index, E, nu, h],  
    filename = 'postPulloutCylinderM40E678.msh';
elseif 41 == mIndex   % incompressible Neo-Hookean material
    A10    = E/2/(1+nu)/2;
    mat    = [mIndex,0, A10, h];  % mat = [index, K, A10, h],
    filename = 'postPulloutCylinderM41E678.msh';
elseif 42 == mIndex   % compressible Neo-Hookean material
    lamda  = E*nu/(1+nu)/(1-2*nu); 
    mu     = E/2/(1+nu);
    mat    = [mIndex,mu,lamda,h];  % mat = [index, mu, lamda, h], 
    filename = 'postPulloutCylinderM42E678.msh';
end

[tsteps,usteps] = readMshFile(filename);
u = usteps{1,end};
vmesh = getVisualMesh(mesh,u);
type = 4;  % 1-ux; 2-uy; 3-uz; 4-magnitude, 
cdata = sqrt(vmesh.u(:,1).^2+vmesh.u(:,2).^2+vmesh.u(:,3).^2);
colorMin = min(cdata);
colorMax = max(cdata);

gg = 5:5:25;
gg = 46;
for k = 1:length(gg)
    figure(k+1)
    u = usteps{1,gg(k)};
    mesh2 = mesh;
    mesh2.coords(:,1:3) = mesh2.coords(:,1:3) + u;
    vmesh = getVisualMesh(mesh2,u);
    vMesh{1,1} = vmesh;
    vMesh{1,2} = vmesh;
    vMesh{1,2}.node(:,2) = L - vMesh{1,2}.node(:,2);
    vMesh{1,2}.u(:,2)    = - vMesh{1,2}.u(:,2);
    vMesh{1,3} = vmesh;
    vMesh{1,3}.node(:,3) = - vMesh{1,3}.node(:,3);
    vMesh{1,3}.u(:,3)    = - vMesh{1,3}.u(:,3);
    vMesh{1,4} = vMesh{1,2}; 
    vMesh{1,4}.node(:,3) = - vMesh{1,4}.node(:,3); 
    vMesh{1,4}.u(:,3)    = - vMesh{1,4}.u(:,3);
    vMesh{1,5} = vmesh; 
    vMesh{1,5}.node(:,1) = - vMesh{1,5}.node(:,1); 
    vMesh{1,5}.u(:,1)    = - vMesh{1,5}.u(:,1);
    vMesh{1,6} = vMesh{1,2}; 
    vMesh{1,6}.node(:,1) = - vMesh{1,6}.node(:,1); 
    vMesh{1,6}.u(:,1)    = - vMesh{1,6}.u(:,1);
    vMesh{1,7} = vMesh{1,5}; 
    vMesh{1,7}.node(:,3) = - vMesh{1,7}.node(:,3); 
    vMesh{1,7}.u(:,3)    = - vMesh{1,7}.u(:,3);
    vMesh{1,8} = vMesh{1,6}; 
    vMesh{1,8}.node(:,3) = - vMesh{1,8}.node(:,3); 
    vMesh{1,8}.u(:,3)    = - vMesh{1,8}.u(:,3);
    
    for i = 1:8
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
    view(45,20);
    axis equal;
    axis on;
    axis([-6,6,-1,12,-10,10]);
    mycolor = abaqusColorMap(12);
    colormap(mycolor);
    C  = colorbar;
    caxis manual
    caxis([colorMin colorMax]);
    set(C,'Ticks',linspace(colorMin,colorMax,13));
    set(get(C,'title'),'string','Um');
end




