%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%
%  This example demonstrates the hyperelastic analysis of KL shell
%  - Pulling of the cylinder
%  - Only one-eighth part is analyzed due to symmetric property
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
% filename = 'geo_pulloutCylinder_nurbs.iga';
filename = 'geoPulloutCylinder1.iga';
geo  = readIgaFile(filename); 
mesh = buildIgaMesh( geo );

% scale the model into the range [0,R]x[0,L/2]x[0,R]
R = 4.953; L  = 10.35;

% plot scaled model
figure(1)
plotIgaMesh(mesh,1);
view(60,30);
%****************************************************************%
%                   Enforce boundary conditions                  %
%****************************************************************%
fprintf("Enforce boundary conditions \n");
pLoad = 10000;  % total load 
[dbc,sym_dbc,tbc,nodeA,nodeB,nodeC] = setBcPulloutCylinder(mesh,pLoad);
h    = 0.094;
nu   = 0.3125;             % Poisson ratio
E    = 10.5e6;
A10  = E/2/(1+nu)/2;
mIndex = 40;
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

timeInterval = 0.02; % initial interval between load steps
maxreit      = 6;    % maximum times of load step reduction 
maxit        = 30;   % maximum iterative steps 
minInterval  = 1e-6; % minimum step interval
divgRatio    = 1e5;  % iteration diverged extremely
init  = [timeInterval,maxreit,maxit,minInterval,divgRatio];
tic
[tsteps, usteps] = nligaKlShell( mesh, mat, dbc, sym_dbc, tbc, init );
toc

% plot displacement of the whole model at the last step
figure(2)
u = reshape(usteps{1,end},3,[])';
plotKLshellCylinderResults(mesh,u,L);
view(40,20);
axis([-2,2,-1,11,-8,8]);
axis off

% plot radial displacement history at the point A and B
figure(3)
uA = zeros(length(usteps),1); 
uB = zeros(length(usteps),1);
uC = zeros(length(usteps),1);
for i = 1:length(usteps)
    uA(i,1) = usteps{1,i}((nodeA-1)*3+3);
    uB(i,1) = usteps{1,i}((nodeB-1)*3+1);
    uC(i,1) = usteps{1,i}((nodeC-1)*3+1);
end
plot([0;uA],[0;tsteps*40000],'r-o');
hold on;
plot([0;-uB],[0;tsteps*40000],'b-s');
plot([0;-uC],[0;tsteps*40000],'k-^');
legend('-u_A','u_B','u_C');
title('Radial displacement history of the point A, B and C');
xlabel('Radial displacement');
ylabel('Load');

% output results
filename = 'postPulloutCylinder';
filename = [filename,'M',num2str(mat(1)),'E',num2str(mesh.nElems)];
fname    = getOutputFileName(filename);
fout     = fopen(fname,'w'); 
writeMshFile(usteps,tsteps,fout);
fclose(fout);

