function mesh = buildIgaMeshFromNurbs(nurbs)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%
%  Bild iga mesh from nurbs based model
%  Here we only consider surface model at this stage!
%
%  ---------------------------------------
%
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%  - 2022
%
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
[Ce, nb] = bezierExtractionOperator(nurbs.knots{1},nurbs.knots{1},nurbs.knots{2},q,W,k);
mesh.nCpts  = geo.nodeN;
mesh.coords = geo.nodes;
mesh.nElems = geo.elemN;
mesh.nElemCpts = zeros(geo.elemN,1);
mesh.elNodeCnt = cell(geo.elemN,1);
mesh.elExtOpe  = cell(geo.elemN,1);
mesh.elDegree  = zeros(geo.elemN,mesh.dim);





end

