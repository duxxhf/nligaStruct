function [dbc,sym_dbc,tbc,nodeA] = setBcPinchedCylinderLineLoad(mesh, tLoad)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%
%  Set boundary conditions for kl shell (pinched cylinder with line load) example
%  Input:
%    mesh - iga mesh structure, using the function 'build_tiga_mesh'
%  Output:
%    dbc      - displacement boundary conditions
%    sym_dbc  - symmetric displacement boundary conditions
%    tbc      - external force boundary conditions
%    nodeA    - index number of control point A
%    nodeB    - index number of control point B
%
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%  - 2022
%
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
for i = 1:length(mesh.nodeSets)
    if strcmp(mesh.nodeSets{i,1}.name,'bottom')
        bottomNodes = mesh.nodeSets{i,1}.gloInx;
    elseif strcmp(mesh.nodeSets{i,1}.name,'nextBottom')
        nextBottomNodes = mesh.nodeSets{i,1}.gloInx;
    elseif strcmp(mesh.nodeSets{i,1}.name,'top')
        topNodes = mesh.nodeSets{i,1}.gloInx;
    elseif strcmp(mesh.nodeSets{i,1}.name,'nextTop')
        nextTopNodes = mesh.nodeSets{i,1}.gloInx;
    elseif strcmp(mesh.nodeSets{i,1}.name,'left')
        leftNodes =mesh.nodeSets{i,1}.gloInx;
    elseif strcmp(mesh.nodeSets{i,1}.name,'nextLeft')
        nextLeftNodes = mesh.nodeSets{i,1}.gloInx;
    elseif strcmp(mesh.nodeSets{i,1}.name,'right')
        rightNodes = mesh.nodeSets{i,1}.gloInx;
    elseif strcmp(mesh.nodeSets{i,1}.name,'nextRight')
        nextRightNodes = mesh.nodeSets{i,1}.gloInx;
    end
end

% sort boundary nodes
bottomNodes1 = bottomNodes(mesh.coords(bottomNodes,3) >= 0);
bottomNodes2 = setdiff(bottomNodes,bottomNodes1);
bottomNodes1 = sortrows([floor(mesh.coords(bottomNodes1,1:3)*1000)/1000, bottomNodes1],[1,3],{'ascend' 'descend'});
bottomNodes1 = bottomNodes1(:,4);
bottomNodes2 = sortrows([floor(mesh.coords(bottomNodes2,1:3)*1000)/1000, bottomNodes2],[1,3],{'descend' 'descend'});
bottomNodes2 = bottomNodes2(:,4);
bottomNodes  = [bottomNodes1; bottomNodes2];

nextBottomNodes1 = nextBottomNodes(mesh.coords(nextBottomNodes,3) >= 0);
nextBottomNodes2 = setdiff(nextBottomNodes,nextBottomNodes1);
nextBottomNodes1 = sortrows([floor(mesh.coords(nextBottomNodes1,1:3)*1000)/1000, nextBottomNodes1],[1,3],{'ascend' 'descend'});
nextBottomNodes1 = nextBottomNodes1(:,4);
nextBottomNodes2 = sortrows([floor(mesh.coords(nextBottomNodes2,1:3)*1000)/1000, nextBottomNodes2],[1,3],{'descend' 'descend'});
nextBottomNodes2 = nextBottomNodes2(:,4);
nextBottomNodes  = [nextBottomNodes1; nextBottomNodes2];

topNodes1 = topNodes(mesh.coords(topNodes,3) >= 0);
topNodes2 = setdiff(topNodes,topNodes1);
topNodes1 = sortrows([floor(mesh.coords(topNodes1,1:3)*1000)/1000, topNodes1],[1,3],{'ascend' 'descend'});
topNodes1 = topNodes1(:,4);
topNodes2 = sortrows([floor(mesh.coords(topNodes2,1:3)*1000)/1000, topNodes2],[1,3],{'descend' 'descend'});
topNodes2 = topNodes2(:,4);
topNodes  = [topNodes1; topNodes2];

nextTopNodes1 = nextTopNodes(mesh.coords(nextTopNodes,3) >= 0);
nextTopNodes2 = setdiff(nextTopNodes,nextTopNodes1);
nextTopNodes1 = sortrows([floor(mesh.coords(nextTopNodes1,1:3)*1000)/1000, nextTopNodes1],[1,3],{'ascend' 'descend'});
nextTopNodes1 = nextTopNodes1(:,4);
nextTopNodes2 = sortrows([floor(mesh.coords(nextTopNodes2,1:3)*1000)/1000, nextTopNodes2],[1,3],{'descend' 'descend'});
nextTopNodes2 = nextTopNodes2(:,4);
nextTopNodes  = [nextTopNodes1; nextTopNodes2];

leftNodes = sortrows([mesh.coords(leftNodes,1:3), leftNodes],2);
leftNodes = leftNodes(:,4);

nextLeftNodes = sortrows([mesh.coords(nextLeftNodes,1:3), nextLeftNodes],2);
nextLeftNodes = nextLeftNodes(:,4);

rightNodes = sortrows([mesh.coords(rightNodes,1:3), rightNodes],2);
rightNodes = rightNodes(:,4);

nextRightNodes = sortrows([mesh.coords(nextRightNodes,1:3), nextRightNodes],2);
nextRightNodes = nextRightNodes(:,4);

% find A nodes
nodeA = leftNodes(1);

% plot constrained control points
hold on;
plot3(mesh.coords(leftNodes,1),mesh.coords(leftNodes,2),mesh.coords(leftNodes,3),'mo');
plot3(mesh.coords(nextLeftNodes,1),mesh.coords(nextLeftNodes,2),mesh.coords(nextLeftNodes,3),'m*');
plot3(mesh.coords(rightNodes,1),mesh.coords(rightNodes,2),mesh.coords(rightNodes,3),'bo');
plot3(mesh.coords(nextRightNodes,1),mesh.coords(nextRightNodes,2),mesh.coords(nextRightNodes,3),'b*');
plot3(mesh.coords(bottomNodes,1),mesh.coords(bottomNodes,2),mesh.coords(bottomNodes,3),'co');
plot3(mesh.coords(nextBottomNodes,1),mesh.coords(nextBottomNodes,2),mesh.coords(nextBottomNodes,3),'c*');
plot3(mesh.coords(topNodes,1),mesh.coords(topNodes,2),mesh.coords(topNodes,3),'ro');
plot3(mesh.coords(nextTopNodes,1),mesh.coords(nextTopNodes,2),mesh.coords(nextTopNodes,3),'r*');
view(137,37);

% impose displacement boundary conditions
dbc = [];    % dbc = [node index, direction, prescribed displacement]
% xz plane symmetric, uy = 0 
dbc = [dbc; topNodes,        2*ones(size(topNodes)),       zeros(size(topNodes))];      
% yz plane symmetric, ux = 0
dbc = [dbc; leftNodes,       1*ones(size(leftNodes)),      zeros(size(leftNodes))];     
% symmetric and supported, ux = uz = 0
dbc = [dbc; rightNodes,      1*ones(size(rightNodes)),     zeros(size(rightNodes))];     
dbc = [dbc; rightNodes,      3*ones(size(rightNodes)),     zeros(size(rightNodes))];    
dbc = [dbc; rightNodes(1),   2,    0];     

sym_dbc = [topNodes,nextTopNodes; leftNodes, nextLeftNodes; rightNodes, nextRightNodes];

forceElem = [];
for i = 1:mesh.nElems
    if ~isempty( intersect(mesh.elNodeCnt{i,1},leftNodes) )
        forceElem = [forceElem; i];
    end
end

F = zeros(3*mesh.nCpts,1);
for i = 1:length(forceElem)
    e      = forceElem(i);
    sctr   = mesh.elNodeCnt{e,:};
    elCpts = mesh.coords(sctr,1:3);    
    pu     = mesh.elDegree(e,1);
    pv     = mesh.elDegree(e,2);
    Ce     = mesh.elExtOpe{e,1};
    we     = mesh.coords(sctr,4); % Tspline control points' weights
    sctrB  = zeros(1, (pv+1)*3); 
    sctr2  = abs(elCpts(:,1))<1e-6;
    sctrB(1:3:end) = 3*(sctr(sctr2)-1)+1;
    sctrB(2:3:end) = 3*(sctr(sctr2)-1)+2;
    sctrB(3:3:end) = 3*(sctr(sctr2)-1)+3;
    [gp, wg] = gaussQuadrature(pv+1);        % calculate integration points and its weights
    for ipt = 1:length(gp)
        pt     = parameterGaussMapping( [0,1], gp(ipt) );   % gauss integration mapping
        [R,dR] = computeTsplineBasisDers([pu,pv],[0,pt],Ce,we);
        R      = R(sctr2);
        dR     = dR(2,sctr2);
        j2     = dR*elCpts(sctr2,1:3); 
        j2     = norm(j2); 
        N      = zeros(3,(pv+1)*3);
        N(1,1:3:end) = R;
        N(2,2:3:end) = R;
        N(3,3:3:end) = R;
        fac = 0.5 * j2 * wg(ipt);
        f        = tLoad/2/0.3;
        F(sctrB) = F(sctrB) + N'* [0; 0; f] * fac;
    end
end
tbc = [];        % tbc = [node index, node dof, prescribed force]
for i = 1:length(F)
    if abs(F(i)) > 1e-12
        rem = mod(i,3);   % direction
        quo = (i-rem)/3;  % node number
        if rem == 0, rem = 3; end
        tbc = [tbc; quo, rem, -F(i)];
    end
end
forceNode = tbc(:,1);
plot3(mesh.coords(forceNode,1),mesh.coords(forceNode,2),mesh.coords(forceNode,3),'kd','MarkerSize',10.0);
end
