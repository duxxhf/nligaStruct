function vmesh = buildVisualMesh(num1,num2,mesh)
% Construct visual mesh
trinum = mesh.nElems*num1*num2*2;        % total elements
trimesh = zeros(trinum,3);   % mesh connectivity
linmesh = cell(mesh.nElems,1);
count = 1;
for k = 1:mesh.nElems
    for j = 1:num2
        for i = 1:num1
            p1 = (k-1)*(num1+1)*(num2+1)+(j-1)*(num1+1)+i;
            p2 = (k-1)*(num1+1)*(num2+1)+(j-1)*(num1+1)+i+1;
            p3 = (k-1)*(num1+1)*(num2+1)+(j-1)*(num1+1)+i+1+(num1+1);
            p4 = (k-1)*(num1+1)*(num2+1)+(j-1)*(num1+1)+i+(num1+1);
            trimesh(count,:) = [p1 p2 p4];
            trimesh(count+1,:) = [p2 p3 p4];
            count = count+2;
        end
    end
    linmesh{k,1} = zeros(num1*2+num2*2+1,1);
    linmesh{k,1}(1:num1) = (1:num1);
    linmesh{k,1}((num1+1):(num1+num2)) =  (num1+1):num1+1:(num1+1)*num2;
    linmesh{k,1}((num1+num2+1):(num1*2+num2)) = (num1+1)*(num2+1):-1:(num1+1)*num2+2;
    linmesh{k,1}((num1*2+num2+1):num1*2+num2*2) = ((num1+1)*num2+1) : -(num1+1) : (num1+2);
    linmesh{k,1} = linmesh{k,1} + (k-1)*(num1+1)*(num2+1);
    linmesh{k,1}(end) = linmesh{k,1}(1);
end
vmesh.linmesh = linmesh;
offset = 0;
uu = linspace(0+offset,1-offset,num1+1);
vv = linspace(0+offset,1-offset,num2+1);
count = 1;
tripts = zeros((num1+1)*(num2+1),2);

for j = 1:num2+1
    for i = 1:num1+1
        tripts(count,:) = [uu(i), vv(j)];   % parametric nodal points
        count = count+1;
    end
end
vmesh.node = zeros(mesh.nElems*(num1+1)*(num2+1),3);
vmesh.element = trimesh;
vmesh.tripts  = tripts;

end

