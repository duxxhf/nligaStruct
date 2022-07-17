function [tsteps,usteps] = readMshFile(filename)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%
% Read *.msh file, which is used to save displacements obtained at 
% each converged load step.
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
tline = fgetl(fid); 
nodeN = sscanf(tline(10:end),'%f');
tline = fgetl(fid);  
nsteps = sscanf(tline(7:end),'%f');
usteps = cell(1,nsteps);
tsteps = zeros(nsteps,1);
for i = 1:nsteps
    usteps{1,i} = zeros(nodeN,3);
    tline = fgetl(fid);  
    idx = find(isspace(tline)); 
    tsteps(i) = sscanf(tline(idx(1)+6:end),'%f')';
    for j = 1:nodeN
        tline = fgetl(fid);  
        usteps{1,i}(j,:) = sscanf(tline,'%f')';
    end
end

end

