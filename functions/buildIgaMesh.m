function mesh = buildIgaMesh( geo )
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%
%  Construct mesh structure for bezier geometry obtained using bezier
%  extraction operation
%
%  Input:
%    geo - NURBS geometry structure
%  Output:
%    mesh - mesh structure
%           mesh.dim  - dimension. 1- one dimensional(curve); 2- two dimensional(surface);
%                                  3- three dimensional(volume)
%           mesh.nCpts - total number of control points          
%           mesh.coords - coordinates of control points, [x1 y1 z1 w1; x2 y2 z2 w2 ...]                       
%           mesh.nElems - total number of elements
%           mesh.nElemCpts -  total number of control points in each element
%           mesh.elNodeCnt - element node connectivity
%
%  ---------------------------------------
%
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%  - 2021
%
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
if strcmp(geo.form,'IGA')  % .iga form files
    if strcmp(geo.type,'plane') || strcmp(geo.type,'surface')
        mesh.dim = 2;
    end
    mesh.nCpts  = geo.nodeN;
    mesh.coords = geo.nodes;
    mesh.nElems = geo.elemN;
    mesh.nElemCpts = zeros(geo.elemN,1);
    mesh.elNodeCnt = cell(geo.elemN,1);
    mesh.elExtOpe  = cell(geo.elemN,1);
    mesh.elDegree  = zeros(geo.elemN,mesh.dim);
    for i = 1:geo.elemN
        mesh.nElemCpts(i) = geo.elems{i,1}.n;
        mesh.elNodeCnt{i,1} = geo.elems{i,1}.gloInx +1;
        mesh.elExtOpe{i,1}  = geo.elems{i,1}.extOpe;
        mesh.elDegree(i,:)  = geo.elems{i,1}.degrees;
    end
    mesh.nodeSets = geo.nodeSets;
    mesh.sideSets = geo.sideSets;
    mesh.elemSets = geo.elemSets;
    if ~isempty(mesh.nodeSets)
        for i = 1:length(mesh.nodeSets)
            mesh.nodeSets{i,1}.gloInx = mesh.nodeSets{i,1}.gloInx + 1;
        end
    end
    if ~isempty(mesh.sideSets)
        for i = 1:length(mesh.sideSets)
            mesh.sideSets{i,1}.gloInx = mesh.sideSets{i,1}.gloInx + 1;
        end
    end
    if ~isempty(mesh.elemSets)
        for i = 1:length(mesh.elemSets)
            mesh.elemSets{i,1}.gloInx = mesh.elemSets{i,1}.gloInx + 1;
        end
    end
elseif strcmp(geo.form,'B-NURBS')  % nurbs structure
    mesh.dim = 2;
    U = geo.knots{1};   V = geo.knots{2};
    p = geo.order(1)-1; q = geo.order(2)-1;
    [Ce, nb] = bezierExtractionOperator(U,p,V,q);
    mesh.nCpts  = size(geo.coefs,2)*size(geo.coefs,3);
    mesh.coords = zeros(mesh.nCpts,4);
    for j = 1:size(geo.coefs,3)
        for i = 1:size(geo.coefs,2)
            geo.coefs(1:3,i,j) = geo.coefs(1:3,i,j)./geo.coefs(4,i,j);
            idx = (j-1)*size(geo.coefs,2) + i;
            mesh.coords(idx,:) = geo.coefs(:,i,j);
        end
    end
    mesh.nElems = nb;
    mesh.nElemCpts = zeros(nb,1);
    mesh.elNodeCnt = cell(nb,1);
    mesh.elExtOpe  = cell(nb,1);
    mesh.elDegree  = zeros(nb,mesh.dim);
    for i = 1:nb
        mesh.nElemCpts(i)   = geo.order(1)*geo.order(2);
        mesh.elExtOpe{i,1}  = Ce(:,:,i);
        mesh.elDegree(i,:)  = geo.order-1;
    end

    [elNodeCntU, ~] = build_knot_connectivity( geo.knots{1} );
    [elNodeCntV, ~] = build_knot_connectivity( geo.knots{2} );     
    count = 0;
    for j = 1:size(elNodeCntV,1)
        for i = 1:size(elNodeCntU,1)
            count = count + 1;
            for hh = 1:geo.order(2)
                for gg = 1:geo.order(1)
                    qq = (hh-1)*geo.order(1) + gg;
                    % build element-node connectivity
                    mesh.elNodeCnt{count,1}(qq) = (elNodeCntV(j,hh)-1)*geo.number(1) + elNodeCntU(i,gg);   
                end
            end
        end
    end    
     
    mesh.nodeSets = cell(8,1);
    mesh.nodeSets{1,1}.setType = 'node';
    mesh.nodeSets{1,1}.name    = 'bottom';
    mesh.nodeSets{1,1}.nodeNum = geo.number(1);
    mesh.nodeSets{1,1}.gloInx  = 1:geo.number(1);
    mesh.nodeSets{1,1}.gloInx  = mesh.nodeSets{1,1}.gloInx';
    
    mesh.nodeSets{2,1}.setType = 'node';
    mesh.nodeSets{2,1}.name    = 'nextBottom';
    mesh.nodeSets{2,1}.nodeNum = geo.number(1);
    mesh.nodeSets{2,1}.gloInx  = geo.number(1)+1:geo.number(1)*2;
    mesh.nodeSets{2,1}.gloInx  = mesh.nodeSets{2,1}.gloInx';
    
    mesh.nodeSets{3,1}.setType = 'node';
    mesh.nodeSets{3,1}.name    = 'top';
    mesh.nodeSets{3,1}.nodeNum = geo.number(1);
    mesh.nodeSets{3,1}.gloInx  = (geo.number(2)-1)*geo.number(1)+1:geo.number(1)*geo.number(2);
    mesh.nodeSets{3,1}.gloInx  = mesh.nodeSets{3,1}.gloInx';
   
    mesh.nodeSets{4,1}.setType = 'node';
    mesh.nodeSets{4,1}.name    = 'nextTop';
    mesh.nodeSets{4,1}.nodeNum = geo.number(1);
    mesh.nodeSets{4,1}.gloInx  = (geo.number(2)-2)*geo.number(1)+1:(geo.number(2)-1)*geo.number(1);
    mesh.nodeSets{4,1}.gloInx  = mesh.nodeSets{4,1}.gloInx';
    
    mesh.nodeSets{5,1}.setType = 'node';
    mesh.nodeSets{5,1}.name    = 'left';
    mesh.nodeSets{5,1}.nodeNum = geo.number(2);
    mesh.nodeSets{5,1}.gloInx  = 1:geo.number(1):geo.number(1)*geo.number(2);
    mesh.nodeSets{5,1}.gloInx  = mesh.nodeSets{5,1}.gloInx';
    
    mesh.nodeSets{6,1}.setType = 'node';
    mesh.nodeSets{6,1}.name    = 'nextLeft';
    mesh.nodeSets{6,1}.nodeNum = geo.number(2);
    mesh.nodeSets{6,1}.gloInx  = 2:geo.number(1):geo.number(1)*geo.number(2);
    mesh.nodeSets{6,1}.gloInx  = mesh.nodeSets{6,1}.gloInx';
    
    mesh.nodeSets{7,1}.setType = 'node';
    mesh.nodeSets{7,1}.name    = 'right';
    mesh.nodeSets{7,1}.nodeNum = geo.number(2);
    mesh.nodeSets{7,1}.gloInx  = geo.number(1):geo.number(1):geo.number(1)*geo.number(2);
    mesh.nodeSets{7,1}.gloInx  = mesh.nodeSets{7,1}.gloInx';
    
    mesh.nodeSets{8,1}.setType = 'node';
    mesh.nodeSets{8,1}.name    = 'nextRight';
    mesh.nodeSets{8,1}.nodeNum = geo.number(2);
    mesh.nodeSets{8,1}.gloInx  = geo.number(1)-1:geo.number(1):geo.number(1)*geo.number(2);
    mesh.nodeSets{8,1}.gloInx  = mesh.nodeSets{8,1}.gloInx';
    
    mesh.sideSets = [];
    mesh.elemSets = [];
end

end







