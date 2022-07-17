function geo = readIgaFile(filename)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%
% Read *.iga file, 
% more detials about the definition of iga file can be
% found in 'M.A. Scott, et. al, An integrated approach to engineering design 
% and analysis using the autodesk t-spline plugin for rhino3d'
%
%  ---------------------------------------
%
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%  - 2022
%
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

fid = fopen(filename, 'r'); % open file
if fid == -1
    error ('Cannot open file %s',filename);
end
geo.form = 'IGA';
tline = fgetl(fid); % read the first line
ln = sscanf(tline(6:end),'%s',1); %
if strcmp(ln,'plane')  % 2d
    geo.type = 'plane';
elseif strcmp(ln,'surface') % 3d
    geo.type = 'surface';
else
    error('Wrong definition of type');
end
tline = fgetl(fid); % read the second line
ln = sscanf(tline,'%s',1); %
if strcmp(ln,'nodeN')  % 2d
    geo.nodeN = sscanf(tline(7:end),'%f');
else
    error('Wrong definition of nodeN');
end
tline = fgetl(fid); % read the third line
ln = sscanf(tline,'%s',1); %
if strcmp(ln,'elemN')  % 2d
    geo.elemN = sscanf(tline(7:end),'%f');
else
    error('Wrong definition of elemN');
end

geo.nodes = zeros(geo.nodeN,4);  % initialize the nodes matrix for saveing nodes' coordinates
for i = 1:geo.nodeN
    tline = fgetl(fid); % read the nodes' line
    ln = sscanf(tline,'%s',1); %
    if strcmp(ln,'node')  % 2d
        geo.nodes(i,:) = sscanf(tline(6:end),'%f')';
    else
        error('Wrong definition of node coordinates');
    end
end

geo.elems = cell(geo.elemN,1);  % initialize the elements cell for saving element structure
for i = 1:geo.elemN
    tline = fgetl(fid); % read the belem line
    ln = sscanf(tline,'%s',1); %
    if strcmp(ln,'belem')  % 2d
        tmp = sscanf(tline(6:end),'%f')';
        geo.elems{i,1}.n = tmp(1);
        geo.elems{i,1}.degrees = tmp(2:3);
    else
        error('Wrong definition of degrees of element');
    end
    tline = fgetl(fid); % read the next line
    geo.elems{i,1}.gloInx = zeros(1,geo.elems{i,1}.n);
    geo.elems{i,1}.gloInx = sscanf(tline,'%f')';  % global indices of each non-zero t-spline basis functions
    geo.elems{i,1}.extOpe = zeros(geo.elems{i,1}.n, (geo.elems{i,1}.degrees(1)+1)*(geo.elems{i,1}.degrees(2)+1)); % element's extraction operator
    for j = 1:geo.elems{i,1}.n
        tline = fgetl(fid); % read the next line
        geo.elems{i,1}.extOpe(j,:) = sscanf(tline,'%f')';
    end
end
geo.nodeSets = [];
geo.sideSets = [];
geo.elemSets = [];
tline = fgetl(fid); % read the set line
if tline == -1
    return;
end
ln = sscanf(tline(1:3),'%s',1); %
nodeSetNum = 0; sideSetNum = 0; elemSetNum = 0;
while strcmp(ln,'set') 
    idx = find(isspace(tline)); % find all spaces and their positions
    if strcmp(tline(idx(2)+1:idx(3)-1),'node')  % save node sets
        nodeSetNum = nodeSetNum + 1;
        geo.nodeSets{nodeSetNum,1}.setType = 'node';
        geo.nodeSets{nodeSetNum,1}.nodeNum = sscanf(tline(idx(1)+1:idx(2)-1),'%f');
        geo.nodeSets{nodeSetNum,1}.name = tline(idx(3)+1:idx(4)-1);
        geo.nodeSets{nodeSetNum,1}.gloInx = sscanf(tline(idx(4)+1:end),'%f');
    elseif strcmp(tline(idx(2)+1:idx(3)-1),'side') % save side sets
        sideSetNum = sideSetNum + 1;
        geo.sideSets{sideSetNum,1}.setType = 'side';
        geo.sideSets{sideSetNum,1}.sideNum = sscanf(tline(idx(1)+1:idx(2)-1),'%f');
        geo.sideSets{sideSetNum,1}.name = tline(idx(3)+1:idx(4)-1);
        geo.sideSets{sideSetNum,1}.gloInx = zeros(geo.sideSets{sideSetNum,1}.sideNum,1);
        geo.sideSets{sideSetNum,1}.gloInxDir = cell(geo.sideSets{sideSetNum,1}.sideNum,1);
        for j = 1:geo.sideSets{sideSetNum,1}.sideNum-1
            geo.sideSets{sideSetNum,1}.gloInx(j) = sscanf(tline(idx(j*2+2)+1:idx(j*2+3)-1),'%f');
            geo.sideSets{sideSetNum,1}.gloInxDir{j,1} = tline(idx(j*2+3)+1:idx(j*2+4)-1);
        end
        j = j+1;
        geo.sideSets{sideSetNum,1}.gloInx(j) = sscanf(tline(idx(j*2+2)+1:idx(j*2+3)-1),'%f');
        geo.sideSets{sideSetNum,1}.gloInxDir{j,1} = tline(idx(j*2+3)+1:end);
    elseif strcmp(tline(idx(2)+1:idx(3)-1),'elem') % save element sets
        elemSetNum = elemSetNum +1;
        geo.elemSets{elemSetNum,1}.setType = 'elem';
        geo.elemSets{elemSetNum,1}.elemNum = sscanf(tline(idx(1)+1:idx(2)-1),'%f');
        geo.elemSets{elemSetNum,1}.name = tline(idx(3)+1:idx(4)-1);
        geo.elemSets{elemSetNum,1}.gloInx = sscanf(tline(idx(4)+1:end),'%f');
    end
    tline = fgetl(fid); % read the set line
    if tline == -1
        return;
    end
    ln = sscanf(tline(1:3),'%s',1); %
end

end

