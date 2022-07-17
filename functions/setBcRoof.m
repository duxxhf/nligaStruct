function [dbc,sym_dbc,nodeA] = setBcRoof(mesh)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%
%  Set boundary conditions for kl shell (roof) example
%  Input:
%    mesh - iga mesh structure, using the function 'build_tiga_mesh'
%    L    - length of the cylinder
%    pLoad - force
%  Output:
%    dbc      - displacement boundary conditions
%    sym_dbc  - symmetric displacement boundary conditions
%    tbc      - external force boundary conditions
%    nodeA    - index number of control point A
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
% sort boundary nodes
bottomNodes = sortrows([floor(mesh.coords(bottomNodes,1)*1000)/1000, bottomNodes],1,'ascend');
bottomNodes = bottomNodes(:,2);
nextBottomNodes = sortrows([floor(mesh.coords(nextBottomNodes,1)*1000)/1000, nextBottomNodes],1,'ascend');
nextBottomNodes = nextBottomNodes(:,2);
topNodes = sortrows([floor(mesh.coords(topNodes,1)*1000)/1000, topNodes],1,'ascend');
topNodes = topNodes(:,2);
nextTopNodes = sortrows([floor(mesh.coords(nextTopNodes,1)*1000)/1000, nextTopNodes],1,'ascend');
nextTopNodes = nextTopNodes(:,2);
leftNodes = sortrows([mesh.coords(leftNodes,2), leftNodes],1);
leftNodes = leftNodes(:,2);
nextLeftNodes = sortrows([mesh.coords(nextLeftNodes,2), nextLeftNodes],1);
nextLeftNodes = nextLeftNodes(:,2);
rightNodes = sortrows([mesh.coords(rightNodes,2), rightNodes],1);
rightNodes = rightNodes(:,2);
nextRightNodes = sortrows([mesh.coords(nextRightNodes,2), nextRightNodes],2);
nextRightNodes = nextRightNodes(:,2);

% find mid node on the free edge
nodeA = topNodes(end);

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
plot3(mesh.coords(nodeA,1),mesh.coords(nodeA,2),mesh.coords(nodeA,3),'kd','MarkerSize',10.0);
view(37,37);

% impose displacement boundary conditions
dbc = [];    % dbc = [node index, direction, prescribed displacement]
% rigid diaphram boundary, ux = uz  = 0
dbc = [dbc; bottomNodes,     1*ones(size(bottomNodes)),     zeros(size(bottomNodes))];  
dbc = [dbc; bottomNodes,     3*ones(size(bottomNodes)),     zeros(size(bottomNodes))];
% top symmetric boundary, uy = 0
dbc = [dbc; topNodes,        2*ones(size(topNodes)),      zeros(size(topNodes))];      
% left symmetric boundary, ux = 0
dbc = [dbc; leftNodes,       1*ones(size(leftNodes)),      zeros(size(leftNodes))];      
sym_dbc = [topNodes,nextTopNodes; leftNodes, nextLeftNodes];
clampedNode = bottomNodes(1);
dbc     = [dbc; clampedNode, 2, 0];

end

