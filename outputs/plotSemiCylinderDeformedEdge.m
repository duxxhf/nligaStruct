%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%  Post-Processing: deformed edges
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
% find the elements contacing with the free edge
bottomElem = []; elemDir = []; centZ = [];
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
        if abs(x(2)) < 1e-6
            bottomElem = [bottomElem; e];
            elemDir    = [elemDir;j];
            centZ      = [centZ;x(3)];
        end
    end
end

centZ = sortrows([centZ,bottomElem,elemDir],1,'descend');
bottomElem = centZ(:,2);
elemDir    = centZ(:,3);


% read msh file for the cylinder with compressible material
filename = 'postSemiCylinderM40E1664.msh';
[tsteps,usteps] = readMshFile(filename);

%% plot deformed free edge
% steps = 0:1:length(tsteps);
steps = [5:5:25,32];
dispt = 0:0.5:1; % parameter points for discretizing the elements of free edge
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

figure(1)
colorIdx = [27,158,119;
            217,95,2;
            117,112,179;
            231,41,138;
            102,166,30;
            230,171,2;
            116,118,29;
            102,102,102;];
% plot deformed free edge
for i  = 1:length(steps)
    rgb = colorIdx(i,:);
    hold on;
    plot(outPutCrvs{i,1}(:,1),outPutCrvs{i,1}(:,3),'-','Color',rgb./255);
    plot(-outPutCrvs{i,1}(:,1),outPutCrvs{i,1}(:,3),'-','Color',rgb./255);
end

axis equal;
box on;
xaxis = [-1.3,1.3]; 
yaxis = [0,0];
axis([xaxis,-0.9,1.2]);
plot(xaxis,yaxis,'--','Color',colorIdx(end,:)./255);
xlabel('x');
ylabel('z');

t = linspace(pi,0,2*size(outPutCrvs{i,1},1));
x = cos(t); y = sin(t);
plot(x,y,'--','Color',colorIdx(end,:)./255);

A = [];
for i=1:length(outPutCrvs)
    A = [A,[-flipud(outPutCrvs{i,1}(:,1)),flipud(outPutCrvs{i,1}(:,3));outPutCrvs{i,1}(:,[1,3])]];
end
A = [A,x',y'];
filename = 'semiCylinderEdges.xlsx';
xlswrite(filename,A);






