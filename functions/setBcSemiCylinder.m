function [dbc,sym_dbc,tbc,nodeA] = setBcSemiCylinder(mesh,pLoad)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%
%  Set boundary conditions for kl shell (semiCylinder) example
%  Input:
%    mesh - iga mesh structure, using the function 'build_tiga_mesh'
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
bottomNodes = sortrows([floor(mesh.coords(bottomNodes,1:3)*1000)/1000, bottomNodes],[1,3],{'ascend' 'descend'});
bottomNodes = bottomNodes(:,4);
nextBottomNodes = sortrows([floor(mesh.coords(nextBottomNodes,1:3)*1000)/1000, nextBottomNodes],[1,3],{'ascend' 'descend'});
nextBottomNodes = nextBottomNodes(:,4);
topNodes = sortrows([floor(mesh.coords(topNodes,1:3)*1000)/1000, topNodes],[1,3],{'ascend' 'descend'});
topNodes = topNodes(:,4);
nextTopNodes = sortrows([floor(mesh.coords(nextTopNodes,1:3)*1000)/1000, nextTopNodes],[1,3],{'ascend' 'descend'});
nextTopNodes = nextTopNodes(:,4);
leftNodes = sortrows([mesh.coords(leftNodes,1:3), leftNodes],2);
leftNodes = leftNodes(:,4);
nextLeftNodes = sortrows([mesh.coords(nextLeftNodes,1:3), nextLeftNodes],2);
nextLeftNodes = nextLeftNodes(:,4);
rightNodes = sortrows([mesh.coords(rightNodes,1:3), rightNodes],2);
rightNodes = rightNodes(:,4);
nextRightNodes = sortrows([mesh.coords(nextRightNodes,1:3), nextRightNodes],2);
nextRightNodes = nextRightNodes(:,4);

% find forced node
tol = 1e-3;
forceNodeA   = leftNodes( abs(mesh.coords(leftNodes,2)-0) < tol );

% find serial number of A node
nodeA = forceNodeA;

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
plot3(mesh.coords(forceNodeA,1),mesh.coords(forceNodeA,2),mesh.coords(forceNodeA,3),'kd','MarkerSize',10.0);
view(137,37);

% impose displacement boundary conditions
dbc = [];    % dbc = [node index, direction, prescribed displacement]
% top side is clamped, ux = uy = uz = 0
dbc = [dbc; topNodes,     1*ones(size(topNodes)),     zeros(size(topNodes))];
dbc = [dbc; topNodes,     2*ones(size(topNodes)),     zeros(size(topNodes))];
dbc = [dbc; topNodes,     3*ones(size(topNodes)),     zeros(size(topNodes))];
dbc = [dbc; nextTopNodes, 1*ones(size(nextTopNodes)), zeros(size(nextTopNodes))];
dbc = [dbc; nextTopNodes, 2*ones(size(nextTopNodes)), zeros(size(nextTopNodes))];
dbc = [dbc; nextTopNodes, 3*ones(size(nextTopNodes)), zeros(size(nextTopNodes))];
 % left side is symmetric, ux = 0
dbc = [dbc; leftNodes,    1*ones(size(leftNodes)),    zeros(size(leftNodes))];     
 % right side is simply supported, uz = 0
dbc = [dbc; rightNodes,   3*ones(size(rightNodes)),   zeros(size(rightNodes))];     

sym_dbc = [leftNodes, nextLeftNodes; rightNodes, nextRightNodes];
% sym_dbc = sym_dbc(1:end-2,:);

tbc = [];        % tbc = [node index, node dof, prescribed force]
force = 1000;
tbc = [tbc; forceNodeA,      3*ones(size(forceNodeA)),    -pLoad*ones(size(forceNodeA))];

end

