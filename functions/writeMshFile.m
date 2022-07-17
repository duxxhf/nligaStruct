function writeMshFile(usteps,tsteps,fout)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%
%  output displacements of kl shell examples and save it into msh file
%
%  Input:
%    usteps  - displacements at each step
%    tsteps  - time at each step
%    fout    - file for saving  
%
%  ---------------------------------------
nctr = length(usteps{1,1})/3;
fprintf(fout,'nCtrlPts %d\n',nctr);  % total control points
fprintf(fout,'nsteps %d\n',length(tsteps));  % total steps
for i = 1:length(tsteps)
    fprintf(fout,'istep=%d, time=%e\n', i, tsteps(i));
    u = reshape(usteps{1,i},3,[])';
    for j = 1:size(u,1)
        fprintf(fout,'%f %f %f\n', u(j,:));
    end
end

end

