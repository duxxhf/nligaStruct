function [tsteps, usteps] = nligaKlShell( mesh, mat, dbc, sym_dbc, tbc, init )
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%
% Main frame for nonlinear kirchhoff-love isogeometric analysis
%  Input:
%    mesh - iga mesh structure
%    mat - material definition
%    dbc - displacements boundary conditions
%    tbc - tractions boundary conditions 
%    fout - output figure handle for visualization
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

% initialization
dof = 3; % degree of freedom
ndofs = dof * mesh.nCpts;   % total dofs
scatdbc = [];
scattbc = [];
if ~isempty(dbc)
    scatdbc = dof * (dbc(:,1)-1) + dbc(:,2);   % scatter dbc
end
if ~isempty(tbc) 
    scattbc = dof * (tbc(:,1)-1) + tbc(:,2);   % scatter tbc
end

alldof = 1:ndofs;
freedof = setdiff(alldof, scatdbc);    % nodes without displacement constraint

sctrs = [];
for i = 1:dof
    sctrs = [sctrs;dof * (sym_dbc-1)+i];
end
nn = size(sctrs, 1);
freedof = setdiff(freedof, sctrs);    % nodes without displacement constraint

tol = 1e-6;               % tolerance of convergence
reit = 0;                 % reduction index
ndbc = size(dbc,1);       % number of displacement constrained nodes
ntbc = size(tbc,1);       % number of displacement constrained nodes
u  = zeros(ndofs,1);      % nodal displacements
cu = zeros(ndofs,1);      % converged nodal displacements
step = 0;                 % load/displacement step index
curtime = 0;              % current time
timeInterval = init(1);   % initial interval between load steps
maxreit      = init(2);   % maximum times of load step reduction 
maxit        = init(3);   % maximum iterative steps 
minInterval  = init(4);   % minimum step interval
divgRatio    = init(5);   % iteration diverged extremely
cnit   = [];              % record the iterative steps
tsteps = [];              % record time at each step
usteps = [];              % record displacements at each step

while curtime ~= 1                      % Achieve the maximum time(load)?
    curtime = curtime + timeInterval;   % Update current time step
    if curtime > 1                      % If current time exceeds 1, amend it
        timeInterval = 1 - curtime + timeInterval;  % Update time interval 
        curtime = 1;                    % Let current time be 1
    end
    err = 1e6;                          % Initialize iterative error
    perr = err;                         % Iterative error at the last iterative step
    nit = 0;                            % Initialize iterative steps
    fprintf(1,'\n \t time   time step   iter \t  residual \n');
    % While the iterative error is larger than the termination tolerance and the current iterative steps is smaller than allowed steps, enter the next loop
    while (err > tol) && (nit <= maxit) 
        nit    = nit+1;                 % Increase iterative steps steps by 1
        % Compute stiffness matrix and internal force
        [Kt,P]  = globalStiffnessHyperKlShell(mesh,mat,u);  
        alpha  = max(diag(Kt))*1e5;      % Penalty factor
        pStiff = alpha*[1 -1;-1 1];      % Penalty stiffness
        for i  = 1:size(sctrs,1)         % Iterate each pair of symmetry conditions
            sctri = sctrs(i,:);          % Extract the i-th pair 
            % Add penalty stiffness to global stifness
            Kt(sctri,sctri) = Kt(sctri,sctri) + pStiff; 
        end
        f      = zeros(ndofs,1);        % Initialize external force
        if ntbc~=0, f(scattbc) = tbc(:,3); end         % Enforce traction forces
        if ndbc~=0                                     % Enforce displacement conditions
            Kt(scatdbc,:)       = zeros(ndbc,ndofs);   % Amend stiffness matrix 
            Kt(scatdbc,scatdbc) = eye(ndbc);           % Amend stiffness matrix 
            f(scatdbc,:)       = 0;                    % Amend external forces
            if nit == 1, f(scatdbc,:) = dbc(:,3); end  % Amend external forces
        end
        r = curtime*f - P;              % Residual vector
        if ndbc~=0, r(scatdbc) = curtime*dbc(:,3) - u(scatdbc); end
        du = Kt\r;                       % Solve equation
        u = u + du;                     % Update displacement
        if nit > 1                              % Compute iterative error
            num = r(freedof)' * r(freedof);     % Norm of residual vector
            denom = 1+f(freedof)' * f(freedof); % Norm of external force vector + 1
            err = num/denom;                    % Iterative error
        end
        % If solution diverge extremely, update nit and perr, then go to next iteration
        fprintf(1,'%10.5f %10.3e %5d %14.5e \n',curtime,timeInterval,nit,err); 
        if err/perr > divgRatio && nit > 2
            nit = maxit+1; % divgRatio is diverge ratio
        else
            perr = err; 
        end
    end
    if  nit <= maxit                    % Converged and nit is smaller than the maximum
        reit = 0;                       % Reset reduction index
        step = step + 1;                % Increase converged steps by 1
        cu = u;                         % Update converged displacements
        cnit = [cnit, nit];             % Record current iterative step
        tsteps = [tsteps;curtime];      % Record current times
        usteps = [usteps,{cu}];         % Record displacements at current step
        % If iterative steps <=5, increase time interval
        if length(cnit) >=2 && all(cnit(end-1:end) <= 5) 
            timeInterval = timeInterval*1.5;    % increase the time interval by times 1.5
        end
    else                                % If it is not converged
        if reit <= maxreit              % Decrease time interval and continue iterating
            curtime = curtime - timeInterval;   % Recover the last time step
            timeInterval = timeInterval/4;      % Decrease time interval
            % If time interval is smaller than the minimum value, terminate iteration
            if timeInterval < minInterval
                return; 
            end  
            reit = reit+1;              % Increase reduction index by 1
            u = cu;                     % Recover the last converged displacements
        else 
            return; 
        end               % Stop iteration
    end
end

end

