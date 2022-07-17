%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%
%  This example demonstrates the compressible hyperelastic analysis of KL shell
%  - Pinched cylinder with line load
%  - Only quarter part is analyzed due to symmetric property
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
filename = 'geoCylinderLineLoad.iga';
geo  = readIgaFile(filename); 
mesh = buildIgaMesh( geo );
R = 0.09; L = 0.3;
% scale the model into the range 
mesh.coords(:,1:3) = 1/100*mesh.coords(:,1:3);
% plot scaled model
figure(1)
plotIgaMesh(mesh,1);
view(50,25);

%****************************************************************%
%                   Enforce boundary conditions                  %
%****************************************************************%
fprintf("Enforce boundary conditions \n");
tLoad = 36000;  % total load 
[dbc,sym_dbc,tbc,nodeA] = setBcPinchedCylinderLineLoad(mesh, tLoad);
h    = 0.2e-2;
nu   = 0.4;             % Poisson ratio
E    = 1.68e11;
mIndex = 41;
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

tic
timeInterval = 0.02; % initial interval between load steps
maxreit      = 6;    % maximum times of load step reduction 
maxit        = 30;   % maximum iterative steps 
minInterval  = 1e-6; % minimum step interval
divgRatio    = 1e3;  % iteration diverged extremely
init  = [timeInterval,maxreit,maxit,minInterval,divgRatio];
[tsteps, usteps] = nligaKlShell( mesh, mat, dbc, sym_dbc, tbc, init );
toc

% plot displacement of the whole model at the last step
figure(2)
u = reshape(usteps{1,end},3,[])';
plotKLshellCylinderLineLoadResults(mesh,u,L);
view(50,30);

% plot radial displacement history at the point A
uA = zeros(length(usteps),1);
for i = 1:length(usteps)
    uA(i,1) = usteps{1,i}((nodeA-1)*3+3);
end
figure(3)
plot(-[0;uA],[0;tsteps*tLoad],'r-o');
legend('-u_A');
title('Vertical displacement history of the point A');
xlabel('Vertical displacement');
ylabel('Load');
x = -0.16;
y = (tsteps(end-1)-tsteps(end-2))*(x-uA(end-2))/(uA(end-1)-uA(end-2)) + tsteps(end-2);
y = y*tLoad;  % the total load for u = 0.16m

% output results
filename = 'postCylinderLineLoad';
filename = [filename,'M',num2str(mat(1)),'E',num2str(mesh.nElems)];
fname    = getOutputFileName(filename);
fout     = fopen(fname,'w'); 
writeMshFile(usteps,tsteps,fout);
fclose(fout);

