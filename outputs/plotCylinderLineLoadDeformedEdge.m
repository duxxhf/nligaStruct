%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%
%  cylinderLineLoad example
%  - we plot the deformed free edge of the cylinder at each step
%  - left part corresponds to incompressible material
%  - right part corresponds to compressible material
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
% read iga file
filename = 'geoCylinderLineLoad2.iga';
geo  = readIgaFile(filename); 
mesh = buildIgaMesh( geo );
R = 0.09; L = 0.3;
mesh.coords(:,1:3) = 1/100*mesh.coords(:,1:3);

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
filename = 'postCylinderLineLoadM42E810.msh';
[tsteps,usteps] = readMshFile(filename);

%% plot load-deflection curve of the point A
for i = 1:length(mesh.nodeSets)
    if strcmp(mesh.nodeSets{i,1}.name,'left')
        leftNodes = mesh.nodeSets{i,1}.gloInx;
    end
end
leftNodes = sortrows([mesh.coords(leftNodes,1:3), leftNodes],2);
leftNodes = leftNodes(:,4);
% find A nodes
nodeA = leftNodes(1);

figure(1)
uA = zeros(length(usteps),1);
for i = 1:length(usteps)
    uA(i,1) = usteps{1,i}(nodeA,3);
end
tLoad = 36000;  % total load 
plot(-[0;uA],[0;tsteps*tLoad],'r-o');
legend('-u_A');
title('Vertical displacement history of the point A');
xlabel('Vertical displacement');
ylabel('Load');
x = -0.16;
y = (tsteps(end-1)-tsteps(end-2))*(x-uA(end-2))/(uA(end-1)-uA(end-2)) + tsteps(end-2);
y = y*tLoad;  % the total load for u = 0.16m
fprintf("Total load at the state ua = -0.16 is %f\n",y);

%% plot deformed free edge
steps = 0:1:length(tsteps);
dispt = 0:0.2:1; % parameter points for discretizing the elements of free edge
profileMsh = zeros(length(steps),length(bottomElem),3*length(dispt));
for i = 1:length(steps)
    if i == 1, u = zeros(size(usteps{1,1}));
    else,      u = usteps{1,steps(i)};
    end
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
        end
    end
end

figure(2)
% plot deformed free edge
for i  = 1:length(steps)
    for j = 1:length(bottomElem)
        hold on;
        x = reshape(profileMsh(i,j,1:3:end),1,[]);
        z = reshape(profileMsh(i,j,3:3:end),1,[]);
        plot(x,z,'-','Color',[27,158,119]./255);
    end
end
axis equal;

% read msh file for the cylinder with incompressible material
filename = 'postCylinderLineLoadM41E810.msh';
[tsteps,usteps] = readMshFile(filename);
steps = 0:1:length(tsteps);
dispt = 0:0.2:1; % parameter points for discretizing the elements of free edge
profileMsh = zeros(length(steps),length(bottomElem),3*length(dispt));
for i = 1:length(steps)
    if i == 1, u = zeros(size(usteps{1,1}));
    else,      u = usteps{1,steps(i)};
    end
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
        end
    end

end

% plot deformed free edge
for i  = 1:length(steps)
    for j = 1:length(bottomElem)
        hold on;
        x = reshape(profileMsh(i,j,1:3:end),1,[]);
        z = reshape(profileMsh(i,j,3:3:end),1,[]);
        plot(-x,z,'-','Color',[217,95,2]./255);
    end
end
axis equal;
box on;
xaxis = [-0.13,0.13]; 
yaxis = [-0.12,0.1];
axis([xaxis,yaxis]);
plot([sum(xaxis)/2,sum(xaxis)/2],yaxis,'--','Color',[102,102,102]./255);
% xlabel('\itx','FontSize',26,'FontName','Times New Roman','FontWeight','bold');
% ylabel('\itz','FontSize',26,'FontName','Times New Roman','FontWeight','bold');
% set(gca,'FontSize',26,'Fontname', 'Times New Roman');




