%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%
%  Five examples demonstrating large deformation of Kirchhoff-Love shell
%  1 - Pinched cylindrical shell with free ends subjected to a line load
%  2 - Pinched cylindrical shell with rigid diaphragms subjected to point
%  loads
%  3 - Pinched hemispherical shell subjected to point loads
%  4 - Pullout of cylindrical shell
%  5 - Pinched semi-cylindrical shell with a clamped edge
%
%  ---------------------------------------
%
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%  - 2022
%
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

function demo
% plot results
prompt = 'Please input a positive integer between 1 and 5:\n';
fprintf('1: Pinched cylindrical shell with free ends subjected to a line load;\n');
fprintf('2: Pinched cylindrical shell with rigid diaphragms subjected to point loads;\n');
fprintf('3: Pinched hemispherical shell subjected to point loads;\n');
fprintf('4: Pullout of cylindrical shell;\n');
fprintf('5: Pinched semi-cylindrical shell with a clamped edge.\n');
num = input(prompt);

% if 0 == nargin
%     num = 1;
% end
if 1 == num
    fprintf('Reading the mesh file of the cylinderLineLoad example\n');
    filename = 'geoCylinderLineLoad2.iga';
    geo  = readIgaFile(filename); 
    mesh = buildIgaMesh( geo );
    L = 0.3;
    % scale the model into the range 
    mesh.coords(:,1:3) = 1/100*mesh.coords(:,1:3);
    
    fprintf('Reading the numerical results of the cylinderLineLoad example\n');
    fprintf('Compressible Neo-Hookean material and 810 bi-cubic Bezier elements are used\n');
    filename = 'postCylinderLineLoadM42E810.msh';
    [~,usteps] = readMshFile(filename);
    u = usteps{1,end};
    
    fprintf('Visualizing the final configuration\n');
    fprintf('Please wait a second...\n');
    plotKLshellCylinderLineLoadResults(mesh,u,L);
    fprintf('Completed!\n');
elseif 2 == num
    fprintf('Reading the mesh file of the pinchedCylinder example\n');
    filename = 'geoPichchedCylinder.iga';
    geo  = readIgaFile(filename); 
    mesh = buildIgaMesh( geo );
    L = 200;
    mesh.coords(:,1:3) = L/60*mesh.coords(:,1:3);
    
    fprintf('Reading the numerical results of the pinchedCylinder example\n');
    fprintf('St. Venant-Kirchhoff material and 2860 bi-cubic Bezier elements are used\n');
    filename = 'postPinchedCylinderM40E2860.msh';
    [~,usteps] = readMshFile(filename);
    
    fprintf('Visualizing the final configuration\n');
    fprintf('Please wait a second...\n');
    u = usteps{1,end};
    plotKLshellCylinderResults(mesh,u,L);
    set(gcf,'renderer','opengl');
    light;
    fprintf('Completed!\n');
elseif 3 == num
    fprintf('Reading the mesh file of the pinchedHemisphere example\n');
    filename = 'geoHemisphere.iga';
    geo  = readIgaFile(filename); 
    mesh = buildIgaMesh( geo );
    
    fprintf('Reading the numerical results of the pinchedHemisphere example\n');
    fprintf('St. Venant-Kirchhoff material and 500 bi-cubic Bezier elements are used\n');
    filename = 'postHemisphereHoleM40E500.msh';
    [~,usteps] = readMshFile(filename);
    
    fprintf('Visualizing the final configuration\n');
    fprintf('Please wait a second...\n');
    u = usteps{1,end};
    plotKLshellHemisphereResults(mesh,u);
    axis([-15,15,-8,8,-5,8]);
    fprintf('Completed!\n');
elseif 4 == num        
    fprintf('Reading the mesh file of the pulloutCylinder example\n');
    filename = 'geoPulloutCylinder1.iga';
    geo  = readIgaFile(filename); 
    mesh = buildIgaMesh( geo );
    
    fprintf('Reading the numerical results of the pulloutCylinder example\n');
    fprintf('Incompressible Neo-Hookean material and 678 bi-cubic Bezier elements are used\n');
    filename = 'postPulloutCylinderM41E678.msh';
    [~,usteps] = readMshFile(filename);
    
    fprintf('Visualizing the final configuration\n');
    fprintf('Please wait a second...\n');
    u = usteps{1,end};
    L  = 10.35;
    plotKLshellCylinderResults(mesh,u,L);
    set(gcf,'renderer','opengl');
    light;
    view(40,20);
    axis([-2,2,-1,11,-8,8]);
    fprintf('Completed!\n');
elseif 5 == num
    fprintf('Reading the mesh file of the semiCylinder example\n');
    filename = 'geoSemicylinder.iga';
    geo  = readIgaFile(filename); 
    mesh = buildIgaMesh( geo );
    
    fprintf('Reading the numerical results of the semiCylinder example\n');
    fprintf('Compressible Neo-Hookean material and 1664 bi-cubic Bezier elements are used\n');
    filename = 'postSemiCylinderM42E1664.msh';
    [~,usteps] = readMshFile(filename);
    
    fprintf('Visualizing the final configuration\n');
    fprintf('Please wait a second...\n');
    u = usteps{1,end};
    plotKLshellSemiCylinderResults(mesh,u);
    axis([-1.5,1.5,-0.1,3.1,-0.8,1.1])
    fprintf('Completed!\n');
end

end

