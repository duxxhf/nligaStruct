function [T,dT,dT2] = computeTsplineBasis2ndDers(k,u,Ce,we)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%
%  Construct tspline basis functions and its 2nd derivatives from bernstein basis 
%
%  Input:
%    k  - degree, for 1D, k is scalar; for 2D, k=[k1,k2]
%    u  - parameter point,for 1D, u is scalar; for 2D, u=[u,v]
%    Ce - elemental bezier extraction operator
%    we - weights vector for tspline comtrol points
%  Output:
%    T - tspline basis
%    dT- tspline basis derivatives
%
%  ---------------------------------------
%
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%  - 3-DEC-2021
%
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
We = diag(we);
[B,dB] = computeBernsteinDers(k,u);
B      = B';
dBdu   = dB(1,:)';  % vector
dBdv   = dB(2,:)';  % vector
dB2    = computeBernstein2ndDers(k,u);
dBdu2  = dB2(1,:)'; % vector
dBdv2  = dB2(2,:)'; % vector
dBduv2 = dB2(3,:)'; % vector
W      = we'*Ce*B;  % scalar
dWdu   = we'*Ce*dBdu;    % scalar
dWdv   = we'*Ce*dBdv;    % scalar
dWdu2  = we'*Ce*dBdu2;   % scalar
dWdv2  = we'*Ce*dBdv2;   % scalar
dWduv2 = we'*Ce*dBduv2;  % scalar

dTdu   = We*Ce*( dBdu/W - B*dWdu/W/W );
dTdv   = We*Ce*( dBdv/W - B*dWdv/W/W );
dT     = [dTdu, dTdv]';

dTdu2  = We*Ce*( dBdu2/W  -  2*dWdu*dBdu/(W^2) -  dWdu2*B/(W^2)  +  2*dWdu*dWdu*B/(W^3) );
dTdv2  = We*Ce*( dBdv2/W  -  2*dWdv*dBdv/(W^2) -  dWdv2*B/(W^2)  +  2*dWdv*dWdv*B/(W^3) );
dTduv2 = We*Ce*( dBduv2/W -  dWdv*dBdu/(W^2)   - dWdu*dBdv/(W^2) -  dWduv2*B/(W^2)  + 2*dWdu*dWdv*B/(W^3) );
dT2    = [dTdu2,dTdv2,dTduv2]';
T      = We*Ce*B/W;
T      = T';
end

