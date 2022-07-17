function b = transV(a)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%
%  Transform 3d vector
%
%  Input:
%    a  - [a1,a2,a3]; 1x3
%  Output:
%    b - [0,a(3),-a(2); -a(3),0,a(1);a(2),-a(1),0;] 3x3
%
%  ---------------------------------------
%
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%  - 18-JAN-2022
%
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

b = [   0,  a(3),  -a(2); 
    -a(3),     0,   a(1);
     a(2), -a(1),     0];

end

