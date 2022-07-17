function roof = geoRoofShellNurbs(L,R,theta,p,q,elemX,elemY)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Build Scordelis-Lo roof SHELL model
% Input:
%   R - radius 
%   L - length
%   theta - angle
%   p,q - degrees
%   elemX,elemY - number of elements 
% 
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

if nargin ~= 7    % default parameters
    R = 25;
    L = 50;
    theta = 40;
    p = 3; 
    q = 3;
    elemX = 20;
    elemY = 20;
end
rad = theta/180*pi;
w = cos(rad);

% define control points and knots vector
coefs = zeros(4,3,2);
coefs(:,:,1) = [-R*sin(rad), 0, R*cos(rad), 1;
                0, 0, R/cos(rad)*w, w;
                R*sin(rad), 0, R*cos(rad), 1;]';
coefs(:,:,2) = [-R*sin(rad), L, R*cos(rad), 1;
                0, L*w, R/cos(rad)*w, w;
                R*sin(rad), L, R*cos(rad), 1;]';
knots{1} = [0 0 0 1 1 1];
knots{2} = [0 0 1 1];

% build nurbs solid by using control points and knots vector
roof = nrbmak(coefs,knots);

% degeree elevate
roof = nrbdegelev(roof,[p-2,q-1]);

% insert knots
RefinementX = elemX-1;    % the number of knots inseted in u direction 
RefinementY = elemY-1;    % the number of knots inseted in v direction
iuknots = 1/(RefinementX+1):1/(RefinementX+1):RefinementX/(RefinementX+1);
ivknots = 1/(RefinementY+1):1/(RefinementY+1):RefinementY/(RefinementY+1);
roof = nrbkntins(roof, {iuknots,ivknots});

end

