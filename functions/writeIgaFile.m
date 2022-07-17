function writeIgaFile(coefs,knots,p,q,filename)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%
% write nurbs model into *.iga file, 
% more detials about the definition of iga file can be
% found in 'M.A. Scott, et. al, An integrated approach to engineering design 
% and analysis using the autodesk t-spline plugin for rhino3d'
%  Input:
%    coefs   - control points,[x,y,z,w]
%    knots   - knot vectors
%    p,q     - degrees
%    fname   - name of iga file for exporting
%
%  ---------------------------------------
%
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%  - 2022
%
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
fname     = fileparts(mfilename('fullpath')); % get current path
index_dir = strfind(fname,'\');
str_temp  = fname(1:index_dir(end));
fname     = [str_temp,'igafiles\', filename, '.iga'];
fout = fopen(fname,'w'); 
fprintf(fout,'type surface\n');
nodeN = size(coefs,3)*size(coefs,2);
fprintf(fout,'nodeN %d\n',nodeN);
[Ce, elemN] = bezierExtractionOperator(knots{1},p,knots{2},q);
fprintf(fout,'elemN %d\n',elemN);
for jj = 1:size(coefs,3)
    for ii = 1:size(coefs,2)
        fprintf(fout,'node %f %f %f %f\n',coefs(:,ii,jj)');
    end
end

[elNodeCntU, ~] = build_knot_connectivity( knots{1} );
[elNodeCntV, ~] = build_knot_connectivity( knots{2} );     
elNodeCnt = zeros( elemN,(p+1)*(q+1) );  % element node connectivity
count = 0;
for j = 1:size(elNodeCntV,1)
    for i = 1:size(elNodeCntU,1)
        count = count + 1;
        for hh = 1:q+1
            for gg = 1:p+1
                qq = (hh-1)*(p+1) + gg;
                elNodeCnt(count,qq) = (elNodeCntV(j,hh)-1)*size(coefs,2) + elNodeCntU(i,gg);   
            end
        end
    end
end 
elNodeCnt = elNodeCnt - 1;

for ii = 1:elemN
    fprintf(fout,'belem %d %d %d\n',(p+1)*(q+1),p,q);
    fprintf(fout,'%d ',elNodeCnt(ii,:));
    fprintf(fout,'\n');
    for jj = 1:(p+1)*(q+1)
        fprintf(fout,'%d ',Ce(jj,:,ii));
        fprintf(fout,'\n');
    end
end


% fprintf(fout,'nCtrlPts %d\n',nctr);  % total control points












end

