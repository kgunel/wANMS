function [xmin,fmin,nFeval,nRef,nExp,nIC,nOC,nShrink,iter,BestCost]=NMS(myfunction, funId, nVar, VarMin, VarMax, maxIter, maxfeval, tol)
%
% Nelder-Mead Simplex Algorithm for
% solving the unconstrained optimization problem:
%         min f(x).
%
% It uses the adaptive parameters introduced in the following paper:
%
%
% Inputs:
%  xinit--initial guess
%  tol--tolerance for termination (Recommended value: 10^-4)
%  max_feval--maximum number of function evaluations
%  myfunction--objective function. (In myfunction(x), x is a ROW vector)
%
% Outputs:
%  xmin--approximate optimal solution at termination. It is a row vector.
%  fmin--minimum function value at termination
%  ct--number of function evaluations at termination
%

x0 = VarMin + (VarMax - VarMin)*rand(1,nVar);
%x0 = ones(1,nVar);
x0=x0(:)';  % x0 is a row vector.
myfunction = fcnchk(myfunction);
dim=max(size(x0)); % dimension of the problem
% set up parameters
alpha=1; beta = 2; gamma=0.5; delta=0.5;
%
% Construct the initial simplex: Large initial simplex is used.
%
scalefactor = min(max(max(abs(x0)),1),10);

D0=eye(dim);
D0(dim+1,:)=(1-sqrt(dim+1))/dim*ones(1,dim);
for i=1:dim+1
    X(i,:)=x0 + scalefactor*D0(i,:);
    [~, FX(i), ~, ~]=feval(myfunction, X(i,:), funId);
end
ct=dim+1;
[FX,I]=sort(FX);
X=X(I,:);

%% Main iteration
BestCost = zeros(1,maxIter);
nRef = 0; nExp=0; nIC=0; nOC=0; nShrink=0;
%gcf=figure;
for iter = 1:maxIter
    %if max(max(abs(X(2:dim+1,:)-X(1:dim,:)))) <= scalefactor*tol
%     if max(max(abs(X(2:dim+1,:)-X(1:dim,:)))) <= tol
%         break;
%     end
    
%     if ct>maxfeval
%         break;
%     end
  
    % Centroid of the dim best vertices
    M = mean(X(1:dim,:));  
    
    xref=(1+alpha)*M- alpha*X(dim+1,:);
    [~, Fref, ~, ~]=feval(myfunction, xref, funId);
    nRef = nRef +1;
    ct=ct+1;
    if Fref<FX(1)
        % expansion
        xexp=(1+alpha*beta)*M-alpha*beta*X(dim+1,:);
        [~, Fexp, ~, ~]=feval(myfunction, xexp, funId);
        ct=ct+1;
        nExp=nExp+1;
        if Fexp < Fref
            X(dim+1,:)=xexp;
            FX(dim+1)=Fexp;
        else
            X(dim+1,:)=xref;
            FX(dim+1)=Fref;
        end
    else
        if Fref<FX(dim)
            % accept reflection point
            X(dim+1,:)=xref;
            FX(dim+1)=Fref;
        else
            if Fref<FX(dim+1)
                % Outside contraction
                xoc=(1+alpha*gamma)*M-alpha*gamma*X(dim+1,:);
                [~, Foc, ~, ~]=feval(myfunction, xoc, funId);
                ct=ct+1;
                nOC = nOC+1;
                if Foc<=Fref
                    X(dim+1,:)=xoc;
                    FX(dim+1)=Foc;
                else
                    % shrink
                    for i=2:dim+1
                        X(i,:)=X(1,:)+ delta*(X(i,:)-X(1,:));
                        [~, FX(i), ~, ~]=feval(myfunction, X(i,:), funId);
                        nShrink=nShrink+1;
                    end
                    ct=ct+dim;
                end
            else
                %inside contraction
                xic=(1-gamma)*M+gamma*X(dim+1,:);
                [~, Fic, ~, ~]=feval(myfunction, xic, funId);
                ct=ct+1;
                nIC=nIC+1;
                if Fic<FX(dim+1)
                    X(dim+1,:)=xic;
                    FX(dim+1)=Fic;
                else
                    % shrink
                    for i=2:dim+1
                        X(i,:)=X(1,:)+ delta*(X(i,:)-X(1,:));
                        [~, FX(i), ~, ~]=feval(myfunction, X(i,:), funId);
                        nShrink=nShrink+1;
                    end
                    ct=ct+dim;
                end
            end
        end
    end
    [FX,I]=sort(FX);
    X=X(I,:);
    BestCost(iter) = FX(1);
    %fprintf('Iteration = %4d\t :\t Best Cost = %5.4e\n', iter, BestCost(iter));
end
xmin=X(1,:);
fmin=FX(1);
nFeval=ct;