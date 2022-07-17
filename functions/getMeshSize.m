function msize = getMeshSize(pts)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%
%  Get the bounding box of a mesh structure
%
%  Input:
%    pts  - mesh nodes' coordinates
%  Output:
%    msize - [xmin, xmax, ymin, ymax, zmin, zmax]
%
%  ---------------------------------------
%
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%  - 6-DEC-2021
%
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

xmin = min(pts(:,1));
xmax = max(pts(:,1));
ymin = min(pts(:,2));
ymax = max(pts(:,2));
zmin = min(pts(:,3));
zmax = max(pts(:,3));
msize = [xmin, xmax, ymin, ymax, zmin, zmax];
end

