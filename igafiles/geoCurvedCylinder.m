function cylinder = geoCurvedCylinder(L)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Build a curved cylinder SHELL model
% 
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
if nargin ~= 1
    L = 10;
end
% define control points and knots vector
coefs = zeros(4,12,2);
coefs(:,:,1) = [0.5, 0, 0.2, 1; 0.6, 0, 0.2, 1;   0.6, 0,   0, 1;
                  0, 0,   0, 1;   0, 0, 0.3, 1;   0.3, 0, 0.3, 1;
                0.3, 0, 0.5, 1;   0, 0, 0.5, 1;     0, 0, 0.8, 1;
                0.6, 0, 0.8, 1; 0.6, 0, 0.6, 1;   0.5, 0, 0.6, 1;]';
coefs(:,:,2) = [0.5, L, 0.2, 1; 0.6, L, 0.2, 1;   0.6, L,   0, 1;
                  0, L,   0, 1;   0, L, 0.3, 1;   0.3, L, 0.3, 1;
                0.3, L, 0.5, 1;   0, L, 0.5, 1;     0, L, 0.8, 1;
                0.6, L, 0.8, 1; 0.6, L, 0.6, 1;   0.5, L, 0.6, 1;]';
knots{1} = [0 0 0 0 0 1/8:1/8:1 1 1 1 1];
knots{2} = [0 0 1 1];

% build nurbs solid by using control points and knots vector
cylinder = nrbmak(coefs,knots);

% degeree elevate
cylinder = nrbdegelev(cylinder,[0,3]);

% insert knots
RefinementX = 63;    % the number of knots inseted in u direction 
RefinementY = 63;    % the number of knots inseted in v direction
iuknots = 1/(RefinementX+1):1/(RefinementX+1):RefinementX/(RefinementX+1);
ivknots = 1/(RefinementY+1):1/(RefinementY+1):RefinementY/(RefinementY+1);
iuknots = setdiff(iuknots,unique(knots{1}));
ivknots = setdiff(ivknots,unique(knots{2}));
cylinder = nrbkntins(cylinder, {iuknots,ivknots});

% for k = 1:2
%     hold on;
%     plot3(coefs(1,:,k),coefs(2,:,k),coefs(3,:,k),'r-o');
% end
% view(0,0);
% axis([0,1,0,L,0,1]);
% axis equal;

end

