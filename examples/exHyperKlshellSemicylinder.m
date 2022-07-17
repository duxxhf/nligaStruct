%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%
%  This example demonstrates the hyperelastic analysis of KL shell
%  - semi cylinder 
%  - Only 1/2 part is analyzed due to symmetric property
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
filename = 'geoSemicylinder.iga';
geo  = readIgaFile(filename); 
mesh = buildIgaMesh( geo );

% plot scaled model
figure(1)
plotIgaMesh(mesh,1);
view(30,15);

%****************************************************************%
%                   Enforce boundary conditions                  %
%****************************************************************%
fprintf("Enforce boundary conditions \n");
pLoad = 1000;  % total load 
[dbc,sym_dbc,tbc,nodeA] = setBcSemiCylinder(mesh,pLoad);
h    = 0.03;
nu   = 0.3;             % Poisson ratio
E    = 2.0685e7;
mIndex = 42;
if 40 == mIndex       % linear St. Venant-Kirchhoff type material.
    mat    = [mIndex,E,nu,h];     % mat = [index, E, nu, h],  
elseif 41 == mIndex   % incompressible Neo-Hookean material
    A10    = E/2/(1+nu)/2;
    mat    = [mIndex, A10, h];    % mat = [index, A10, h],
elseif 42 == mIndex   % compressible Neo-Hookean material
    lamda  = E*nu/(1+nu)/(1-2*nu); 
    mu     = E/2/(1+nu);
    mat    = [mIndex,mu,lamda,h];  % mat = [index, mu, lamda, h], 
end

global gp wg gpz wgz
gp_x = max(mesh.elDegree(:,1))+1;           % number of integration points in x-direction
gp_y = max(mesh.elDegree(:,2))+1;           % number of integration points in y-direction
gp_z = 3;                                   % number of integration points in z-direction
[gp,  wg]  = gaussQuadrature(gp_x, gp_y);   % calculate integration points and its weights
[gpz, wgz] = gaussQuadrature(gp_z);         % calculate integration points and its weights


timeInterval = 0.01; % initial interval between load steps
maxreit      = 6;    % maximum times of load step reduction 
maxit        = 30;   % maximum iterative steps 
minInterval  = 1e-6; % minimum step interval
divgRatio    = 1e3;  % iteration diverged extremely
init  = [timeInterval,maxreit,maxit,minInterval,divgRatio];
tic
[tsteps, usteps] = nligaKlShell( mesh, mat, dbc, sym_dbc, tbc, init );
toc

% plot displacement of the whole model at the last step
figure(2)
u = reshape(usteps{1,end},3,[])';
plotKLshellSemiCylinderResults(mesh,u);
axis([-1.5,1.5,-0.1,3.1,-0.8,1.1])

% plot radial displacement history at the point A and B
figure(3)
uA = zeros(length(usteps),1); 
for i = 1:length(usteps)
    uA(i,1) = usteps{1,i}((nodeA-1)*3+3);
end
plot([0;-uA],[0;tsteps*2000],'r-o');
title('Displacement history');
xlabel('Displacement');
ylabel('Load');


% output results
filename = 'postSemiCylinder';
filename = [filename,'M',num2str(mat(1)),'E',num2str(mesh.nElems)];
fname    = getOutputFileName(filename);
fout     = fopen(fname,'w'); 
writeMshFile(usteps,tsteps,fout);
fclose(fout);


