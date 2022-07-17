%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%  Post-Processing: deformed edges
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
figure(1)
plotIgaMesh(mesh,3);
view(3);

h    = 0.04;
nu   = 0.3;             % Poisson ratio
E    = 6.825e7;
A10  = E/2/(1+nu)/2;
mIndex = 40;

% find the elements contacing with the free edge
bottomElem = []; elemDir = []; centZ = []; vPt = [];
gt = [0.5,0; 1,0.5; 0.5,1; 0,0.5];
for e = 1:mesh.nElems
    sctr   = mesh.elNodeCnt{e,:};
    elCpts = mesh.coords(sctr,1:3);    
    pu     = mesh.elDegree(e,1);
    pv     = mesh.elDegree(e,2);
    Ce     = mesh.elExtOpe{e,1};
    we     = mesh.coords(sctr,4); % Tspline control points' weights
    for j = 1:4
        R = computeTsplineBasis([pu,pv],gt(j,:),Ce,we);
        x = R*elCpts;
        if abs(x(3)) < 1e-6
            bottomElem = [bottomElem; e];
            elemDir    = [elemDir;j];
            centZ      = [centZ;x(2)];
            vPt        = [vPt; x];
        end
    end
end

plot3(vPt(:,1),vPt(:,2),vPt(:,3),'*');
centZ = sortrows([centZ,bottomElem,elemDir],1,'ascend');
bottomElem = centZ(:,2);
elemDir    = centZ(:,3);


% read msh file for the cylinder with compressible material
filename = 'postHemisphereHoleM40E500.msh';
[tsteps,usteps] = readMshFile(filename);

%% plot deformed free edge
% steps = 0:1:length(tsteps);
steps = [5:5:25,36];
dispt = 0:0.2:1; % parameter points for discretizing the elements of free edge
profileMsh = zeros(length(steps),length(bottomElem),3*length(dispt));
outPutCrvs = cell(length(steps),1);
for i = 1:length(steps)
    u = usteps{1,steps(i)};
    outPutCrvs{i,1} = zeros(length(bottomElem)*length(dispt),3);
    for j = 1:length(bottomElem)
        e      = bottomElem(j);
        sctr   = mesh.elNodeCnt{e,:};
        elCpts = mesh.coords(sctr,1:3);    
        eU     = u(sctr,1:3);
        pu     = mesh.elDegree(e,1);
        pv     = mesh.elDegree(e,2);
        Ce     = mesh.elExtOpe{e,1};
        we     = mesh.coords(sctr,4); % Tspline control points' weights
        sctrB  = zeros(1, (pv+1)*3); 
        sctr2  = abs(elCpts(:,2))<1e-6;
        for ipt = 1:length(dispt)
            if elemDir(j) == 1, gt = [dispt(ipt),0];
            elseif elemDir(j) == 2, gt = [1,dispt(ipt)];
            elseif elemDir(j) == 3, gt = [dispt(ipt),1];
            elseif elemDir(j) == 3, gt = [0,dispt(ipt)];
            end
            R = computeTsplineBasis([pu,pv],gt,Ce,we);
            x = R*(elCpts+eU);    
            profileMsh(i,j,(ipt-1)*3+1:ipt*3) = x;
            outPutCrvs{i,1}((j-1)*length(dispt)+ipt,:) = x;
        end
    end
end

colorIdx = [27,158,119;
            217,95,2;
            117,112,179;
            231,41,138;
            102,166,30;
            230,171,2;
            116,118,29;
            102,102,102;];
A = [];
for i=1:length(outPutCrvs)
    A1 = outPutCrvs{i,1}(:,[1,2]);
    A2 = flipud([-outPutCrvs{i,1}(:,1),outPutCrvs{i,1}(:,2)]);
    A3 = -outPutCrvs{i,1}(:,[1,2]);
    A4 = flipud([outPutCrvs{i,1}(:,1),-outPutCrvs{i,1}(:,2)]);
    A  = [A, [A1;A2;A3;A4]];
end
% A = [A];
figure(2)
for i=1:length(outPutCrvs)
    hold on; 
    plot(A(:,i*2-1),A(:,i*2),'-','Color',colorIdx(i,:)./255);
end
axis equal;
hold on;
t = linspace(0,2*pi,4*size(outPutCrvs{i,1},1));
x = 10*cos(t); y = 10*sin(t);
plot(x,y,'--','Color',colorIdx(end,:)./255);

filename = 'hemisphereEdges.xlsx';
xlswrite(filename,A);



