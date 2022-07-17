function [Pbe, wbe] = computeBezierCtrlPts(Pte, wte, Ce)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%
%  Construct elemental bezier control points and weights 
%
%  Input:
%    Pte    - control points of T-spline element
%    wte    - weights of T-spline control points for the element
%    Ce     - elemental bezier extraction matrix
%  Output:
%    Pbe    - control points of bezier element
%    wbe    - weights of bezier control points for the element
%
%  ---------------------------------------
%
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%  - 6-DEC-2021
%
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

wbe = Ce'* wte;
Pbe = diag(wbe)\Ce'*diag(wte)*Pte;

end

