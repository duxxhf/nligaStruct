%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%
%  This example demonstrates the hyperelastic analysis of KL shell
%  - Pinched cylinder
%  - Only one-eighth part is analyzed due to symmetric property
%
%  ---------------------------------------
%
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%  - 17-JAN-2022
%
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

clc;
clear;
%****************************************************************%
%             Import and plot geometrical model                  %
%****************************************************************%
fprintf("Import and plot geometrical model \n");
filename = 'geoPichchedCylinder.iga';
geo  = readIgaFile(filename); 
mesh = buildIgaMesh( geo );

% scale the model into the range [0,100]x[0,100]x[0,100]
L = 200;
mesh.coords(:,1:3) = L/60*mesh.coords(:,1:3);

% plot scaled model
figure(1)
plotIgaMesh(mesh,1);

%****************************************************************%
%                   Enforce boundary conditions                  %
%****************************************************************%
fprintf("Enforce boundary conditions \n");
pLoad = 12000;
[dbc,sym_dbc,tbc,nodeA,nodeB] = setBcPinchedCylinder(mesh,L,pLoad);
h    = 1.0;
nu   = 0.3;             % Poisson ratio
E    = 3e4;
A10  = E/2/(1+nu)/2;
mIndex = 40;
if 40 == mIndex       % linear St. Venant-Kirchhoff type material.
    mat    = [mIndex,E,nu,h];     % mat = [index, E, nu, h],  
elseif 41 == mIndex   % incompressible Neo-Hookean material
    A10    = E/2/(1+nu)/2;
    mat    = [mIndex,0, A10, h];  % mat = [index, K, A10, h],
end

global gp wg gpz wgz
gp_x = max(mesh.elDegree(:,1))+1;           % number of integration points in x-direction
gp_y = max(mesh.elDegree(:,2))+1;           % number of integration points in y-direction
gp_z = 3;                                   % number of integration points in z-direction
[gp,  wg]  = gaussQuadrature(gp_x, gp_y);   % calculate integration points and its weights
[gpz, wgz] = gaussQuadrature(gp_z);         % calculate integration points and its weights

timeInterval = 0.02; % initial interval between load steps
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
plotKLshellCylinderResults(mesh,u,L);
axis([-150,150,0,200,-100,100]);
view(75,20);

% plot radial displacement history at the point A and B
figure(3)
uA = zeros(length(usteps),1); uB = zeros(length(usteps),1);
for i = 1:length(usteps)
    uA(i,1) = usteps{1,i}((nodeA-1)*3+3);
    uB(i,1) = usteps{1,i}((nodeB-1)*3+1);
end
plot(-[0;uA],[0;tsteps*12000],'r-o');
hold on;
plot([0;uB],[0;tsteps*12000],'b-s');
legend('-u_A','u_B');
title('Radial displacement history of the point A and B');
xlabel('Radial displacement');
ylabel('Load');

% output results
filename = 'postPinchedCylinder';
filename = [filename,'M',num2str(mat(1)),'E',num2str(mesh.nElems)];
fname    = getOutputFileName(filename);
fout     = fopen(fname,'w'); 
writeMshFile(usteps,tsteps,fout);
fclose(fout);
