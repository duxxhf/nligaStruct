function rplate = geoRectangularPlate( a,b,pts,p,q,elemX,elemY )
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Build a rectangular plate.
% model.
% Input:
%   a,b - length and width 
%   pts - coordinates of the bottom-left corner vertex, [x y]
%   p,q - degrees
%   elemX,elemY - number of elements 
% 
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

if nargin ~= 7    % default parameters
    a = 1;
    b = 1;
    pts = [0,0];
    p = 2; 
    q = 2;
    elemX = 20;
    elemY = 20;
end
x = pts(1);
y = pts(2);
% define control points and knots vector
coefs = zeros(4,2,2);
coefs(:,:,1) = [x, y,   0, 1;   x+a, y,   0, 1;]';
coefs(:,:,2) = [x, y+b, 0, 1;   x+a, y+b, 0, 1;]';
knots{1} = [0 0 1 1];
knots{2} = [0 0 1 1];

% build nurbs solid by using control points and knots vector
rplate = nrbmak(coefs,knots);

% degeree elevate
rplate = nrbdegelev(rplate,[p-1,q-1]);

% insert knots
RefinementX = elemX-1;    % the number of knots inseted in u direction 
RefinementY = elemY-1;    % the number of knots inseted in v direction
iuknots = 1/(RefinementX+1):1/(RefinementX+1):RefinementX/(RefinementX+1);
ivknots = 1/(RefinementY+1):1/(RefinementY+1):RefinementY/(RefinementY+1);
rplate = nrbkntins(rplate, {iuknots,ivknots});


end

