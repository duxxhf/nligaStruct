function mycolor = abaqusColorMap(n)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%
%  This function provides the colors approximately mapping from commercial software ABAQUS
%
%  Input:
%    n - color steps
%  Output:
%    mycolor - color rgb
%
%  ---------------------------------------
%
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%  - 17-JAN-2022
%
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

if n == 12  
    mycolor = [0, 0, 225;
               0, 93, 255;
               0, 185, 255;
               0, 255, 232;
               0, 255, 131;
               0, 255, 46;
               46, 255, 0;
               139, 255, 0;
               232, 255, 0;
               255, 185, 0;
               255, 93, 0;
               225, 0, 0;];              
    mycolor = mycolor./255;

elseif n == 24
    mycolor = [  0,   0, 225;
                 0,  44, 255;
                 0,  88, 255;
                 0, 133, 255;
                 0, 187, 216;
                 0, 221, 255 ;
                 0, 255, 243;
                 0, 255, 199;
                 0, 255, 155;
                 0, 255, 110;
                 0, 255,  66;
                 0, 255,  22;
                22, 255,   0;
                66, 255,   0;
               110, 255,   0;
               155, 255,   0;
               199, 255,   0;
               243, 255,   0;
               255, 221,   0;
               255, 177,   0;
               255, 133,   0;
               255,  88,   0;
               255,  44,   0;
               255,   0,   0;];              
    mycolor = mycolor./255;
end


end

