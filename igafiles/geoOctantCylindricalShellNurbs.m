function pinchedCylinder = geoOctantCylindricalShellNurbs(R,L,p,q,elemX,elemY)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Build a pinched cylindrical SHELL model
% Note that we only build one-eight of the classical pinched cylindrical model
% Input:
%   R - radius 
%   L - length
%   p,q - degrees
%   elemX,elemY - number of elements 
%  ---------------------------------------
%  This function should be used with the toolbox NLIGA:
%  <Du, Xiaoxiao, Gang Zhao, Wei Wang, Mayi Guo, Ran Zhang, and Jiaming Yang. 
%  NLIGA: A MATLAB framework for nonlinear isogeometric analysis. 
%  Computer Aided Geometric Design, 2020, 80:101869.>
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

if nargin ~= 6    % default parameters
    R = 300;
    L = 600;
    p = 2; 
    q = 2;
    elemX = 20;
    elemY = 20;
end

theta = 45;
w = cos(theta/180*pi);
% define control points and knots vector
coefs = zeros(4,3,2);
coefs(:,:,1) = [0, 0, R, 1; R*w, 0, R*w, w; R, 0, 0, 1; ]';
coefs(:,:,2) = [0, L/2, R, 1; R*w, L/2*w, R*w, w; R, L/2, 0, 1;]';
knots{1} = [0 0 0 1 1 1];
knots{2} = [0 0 1 1];

% build nurbs solid by using control points and knots vector
pinchedCylinder = nrbmak(coefs,knots);

% degeree elevate
if p < 2 || q < 1
    error('Please input appropriate degrees.');
end
pinchedCylinder = nrbdegelev(pinchedCylinder,[p-2,q-1]);

% insert knots
RefinementX = elemX-1;    % the number of knots inseted in u direction 
RefinementY = elemY-1;    % the number of knots inseted in v direction
iuknots = 1/(RefinementX+1):1/(RefinementX+1):RefinementX/(RefinementX+1);
ivknots = 1/(RefinementY+1):1/(RefinementY+1):RefinementY/(RefinementY+1);
pinchedCylinder = nrbkntins(pinchedCylinder, {iuknots,ivknots});
end

