function [dbc,sym_dbc,tbc,nodeA,nodeB] = setBcHemisphere(mesh,pLoad)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%
%  Set boundary conditions for kl shell (hemisphere with a hole) example
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
bottomNodes = []; nextBottomNodes = []; topNodes = []; nextTopNodes = [];
leftNodes = []; nextLeftNodes = []; rightNodes = []; nextRightNodes = [];
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
if ~isempty(bottomNodes)
    bottomNodes = sortrows([floor(mesh.coords(bottomNodes,1:3)*1000)/1000, bottomNodes],[1,2],{'descend' 'ascend'});
    bottomNodes = bottomNodes(:,4);
end
if ~isempty(nextBottomNodes)
    nextBottomNodes = sortrows([floor(mesh.coords(nextBottomNodes,1:3)*1000)/1000, nextBottomNodes],[1,2],{'descend' 'ascend'});
    nextBottomNodes = nextBottomNodes(:,4);
end
if ~isempty(topNodes)
    topNodes = sortrows([floor(mesh.coords(topNodes,1:3)*1000)/1000, topNodes],[1,2],{'descend' 'ascend'});
    topNodes = topNodes(:,4);
end
if ~isempty(topNodes)
    nextTopNodes = sortrows([floor(mesh.coords(nextTopNodes,1:3)*1000)/1000, nextTopNodes],[1,2],{'descend' 'ascend'});
    nextTopNodes = nextTopNodes(:,4);
end
leftNodes = sortrows([mesh.coords(leftNodes,1:3), leftNodes],3);
leftNodes = leftNodes(:,4);
nextLeftNodes = sortrows([mesh.coords(nextLeftNodes,1:3), nextLeftNodes],3);
nextLeftNodes = nextLeftNodes(:,4);
rightNodes = sortrows([mesh.coords(rightNodes,1:3), rightNodes],3);
rightNodes = rightNodes(:,4);
nextRightNodes = sortrows([mesh.coords(nextRightNodes,1:3), nextRightNodes],3);
nextRightNodes = nextRightNodes(:,4);

% find forced node
tol = 1e-3;
forceNodeA   = leftNodes( abs(mesh.coords(leftNodes,3)-0) < tol );
forceNodeB   = rightNodes( abs(mesh.coords(rightNodes,3)-0) < tol );

% find serial number of A, B and C nodes
nodeA = forceNodeA;
nodeB = forceNodeB;

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
plot3(mesh.coords(forceNodeB,1),mesh.coords(forceNodeB,2),mesh.coords(forceNodeB,3),'kd','MarkerSize',10.0);
view(137,37);

% impose displacement boundary conditions
dbc = [];    % dbc = [node index, direction, prescribed displacement]
% forceNode A, uz = 0
dbc = [dbc; forceNodeA,      3*ones(size(forceNodeA)),    zeros(size(forceNodeA))];      
 % xz plane symmetric, uy = 0
dbc = [dbc; leftNodes,       2*ones(size(leftNodes)),     zeros(size(leftNodes))];     
 % yz plane symmetric, ux = 0
dbc = [dbc; rightNodes,      1*ones(size(rightNodes)),    zeros(size(rightNodes))];     

sym_dbc = [leftNodes, nextLeftNodes; rightNodes, nextRightNodes];

tbc = [];        % tbc = [node index, node dof, prescribed force]
tbc = [tbc; forceNodeA,      1*ones(size(forceNodeA)),    0.5*pLoad*ones(size(forceNodeA))];
tbc = [tbc; forceNodeB,      2*ones(size(forceNodeB)),   -0.5*pLoad*ones(size(forceNodeB))];

end

