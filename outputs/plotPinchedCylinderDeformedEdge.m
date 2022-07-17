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

% find the elements contacing with the free edge
topElem = []; elemDir = []; centZ = []; vPt = [];
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
        if abs(x(2)-100) < 1e-6
            topElem = [topElem; e];
            elemDir    = [elemDir;j];
            centZ      = [centZ;x(3)];
            vPt        = [vPt; x];
        end
    end
end

plot3(vPt(:,1),vPt(:,2),vPt(:,3),'*');
centZ = sortrows([centZ,topElem,elemDir],1,'descend');
topElem = centZ(:,2);
elemDir    = centZ(:,3);


% read msh file for the cylinder with compressible material
filename = 'postPinchedCylinderM40E2860.msh';
[tsteps,usteps] = readMshFile(filename);

%% plot deformed free edge
% steps = 0:1:length(tsteps);
steps = [5,20,50,70,80,90,100,103];
dispt = 0:0.5:1; % parameter points for discretizing the elements of free edge
profileMsh = zeros(length(steps),length(topElem),3*length(dispt));
outPutCrvs = cell(length(steps),1);
for i = 1:length(steps)
    u = usteps{1,steps(i)};
    outPutCrvs{i,1} = zeros(length(topElem)*length(dispt),3);
    for j = 1:length(topElem)
        e      = topElem(j);
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

t = linspace(0,2*pi,4*size(outPutCrvs{i,1},1));
x = 100*cos(t); y = 100*sin(t);

A = [];
for i=1:length(outPutCrvs)
    A1 = outPutCrvs{i,1}(:,[1,3]);
    A2 = flipud([outPutCrvs{i,1}(:,1),-outPutCrvs{i,1}(:,3)]);
    A3 = -outPutCrvs{i,1}(:,[1,3]);
    A4 = flipud([-outPutCrvs{i,1}(:,1),outPutCrvs{i,1}(:,3)]);
    A  = [A, [A1;A2;A3;A4]];
end
A = [A,x',y'];
figure(2)
for i=1:length(outPutCrvs)
    hold on; 
    plot(A(:,i*2-1),A(:,i*2),'-','Color',colorIdx(i,:)./255);
end
plot(x,y,'--','Color',colorIdx(end,:)./255);
axis equal;

filename = 'pinchedCylinderEdges.xlsx';  % output deformed edges
xlswrite(filename,A);






