function hemisphere = geoHemisphericalShellNurbs( R,phi,p,q,elemX,elemY )
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Build a quadrant of hemispherical SHELL model.
% model.
% Input:
%   R - radius 
%   phi - angle
%   p,q - degrees
%   elemX,elemY - number of elements 
% 
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

if nargin ~= 6    % default parameters
    R = 10;
    phi = 18;
    p = 2; 
    q = 2;
    elemX = 20;
    elemY = 20;
end

rad1 = phi/180*pi;
rad2 = (90-phi)/2/180*pi;
w1   = cos(45/180*pi);
w2   = cos((90-phi)/2/180*pi);

% define control points and knots vector
coefs = zeros(4,3,3);
coefs(:,:,1) = [R, 0, 0, 1; R*w1, R*w1, 0, w1; 0, R, 0, 1; ]';
coefs(:,:,2) = [R*w2, 0, R*tan(rad2)*w2, w2; R*w1*w2, R*w1*w2, R*tan(rad2)*w1*w2, w1*w2; 0, R*w2, R*tan(rad2)*w2, w2;]';
coefs(:,:,3) = [R*sin(rad1), 0, R*cos(rad1), 1; 
                R*sin(rad1)*w1, R*sin(rad1)*w1, R*cos(rad1)*w1, w1; 
                0, R*sin(rad1), R*cos(rad1), 1;]';
knots{1} = [0 0 0 1 1 1];
knots{2} = [0 0 0 1 1 1];

% build nurbs solid by using control points and knots vector
hemisphere = nrbmak(coefs,knots);

% degeree elevate
hemisphere = nrbdegelev(hemisphere,[p-2,q-2]);

% insert knots
RefinementX = elemX-1;    % the number of knots inseted in u direction 
RefinementY = elemY-1;    % the number of knots inseted in v direction
iuknots = 1/(RefinementX+1):1/(RefinementX+1):RefinementX/(RefinementX+1);
ivknots = 1/(RefinementY+1):1/(RefinementY+1):RefinementY/(RefinementY+1);
hemisphere = nrbkntins(hemisphere, {iuknots,ivknots});


end

