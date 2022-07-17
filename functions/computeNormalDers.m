function dN = computeNormalDers(dSdu, dSdv, dSdu2, dSdv2, dSduv2)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%
%  Construct normal derivatives at the parameter u
%
%  Input:
%    dSdu  - derivative along u direction
%    dSdv  - derivative along v direction
%  Output:
%    dN- 1st normal derivatives, dN = [dNdu; dNdv]
%
%  ---------------------------------------
%
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%  - 6-DEC-2021
%
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

n_top   = cross(dSdu,dSdv);  %  numerator 
n_bot   = norm(n_top);       %  denominator
n_vec   = n_top/n_bot;       %  unit normal vector
n_top_u = cross(dSdu2,dSdv) + cross(dSdu,dSduv2); % derivative of numerator along u
n_top_v = cross(dSduv2,dSdv) + cross(dSdu,dSdv2); % derivative of numerator along v

n_bot_u = dot(n_top_u,n_vec); % derivative of denominator along u
n_bot_v = dot(n_top_v,n_vec); % derivative of denominator along v

dNdu    = 1/n_bot * (n_top_u - n_bot_u * n_vec ); % derivative of normal along u
dNdv    = 1/n_bot * (n_top_v - n_bot_v * n_vec ); % derivative of normal along v

dN      = [dNdu; dNdv];
end

