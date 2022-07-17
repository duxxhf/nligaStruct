function tR = transR(R)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%
%  Transform the derivatives matrix
%
%  Input:
%    R  - [R1,R2,R3,...Rn]; 1xn
%  Output:
%    tR - [R1,0,0,R2,0,0,...,Rn,0,0;
%          0,R1,0,0,R2,0,...,0,Rn,0;
%          0,0,R1,0,0,R2,...,0,0,Rn];  3x3n
%
%  ---------------------------------------
%
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%  - 18-JAN-2022
%
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
tR = zeros(3,3*length(R));
tR(1,1:3:end) = R;
tR(2,2:3:end) = R;
tR(3,3:3:end) = R;
end

