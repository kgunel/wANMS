function [xmin,fmin,nFeval,nRef,nExp,nIC,nOC,nShrink,iter,BestCost]=ANMS(myfunction, funId, nVar, VarMin, VarMax, maxIter, maxfeval, tol)
%
% Adaptive Nelder-Mead Simplex Algorithm for
% solving the unconstrained optimization problem:
%         min f(x).
%
% It uses the adaptive parameters introduced in the following paper:
%
%              Fuchang Gao and Lixing Han
% "Implementing the Nelder-Mead simplex algorithm with adaptive parameters"
%  Computational Optimization and Applications, 2012, 51: 259-277.
%
% It also uses a relatively large initial simplex (The initial simplex
% used in the numerical experiments in the above paper is small in order
% to compare ANMS and FMINSEARCH).
%
%                        REMARK:
% If you use ANMS as a local search method in combination with a
% metaheuristic method, you may want to use a smaller initial simplex.
%
%  by Fuchang Gao and Lixing Han
%  coded August, 2010
%  slightly updated August, 2015
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
load('x0.mat');
%x0 = VarMin + (VarMax - VarMin)*rand(1,nVar);
%x0 =(1/2).* ones(1,nVar);
%x0=x0(:)';  % x0 is a row vector.
myfunction = fcnchk(myfunction);
dim=max(size(x0)); % dimension of the problem
% set up adaptive parameters
alpha=1; beta=1+2/dim; gamma=0.75-0.5/dim; delta=1-1/dim;
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

% Plot surface
%fig = figure;
%set(fig, 'Position', get(0, 'Screensize'));
%{
hold on
%surf(peaks,'EdgeColor','none','FaceColor','interp','FaceLighting','phong')

syms x y
x = linspace(VarMin,VarMax,100);
y = linspace(VarMin,VarMax,100);
[xx,yy] = meshgrid(x,y);
Z = xx.^2+yy.^2
%surfaceMeshHandle = fsurf(Z,'ShowContours','on');
[M,c] = contour(x,y,Z);
c.LineWidth = 3;
xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex'); zlabel('$z$','Interpreter','latex');
%title(['$\displaystyle f(x,y) = ' latex(Z) '$' ],'Interpreter','latex','Color', 'b','FontSize',14); 
colormap(jet);
axis([VarMin VarMax VarMin VarMax]),

plot(X(:,1),X(:,2),'k*');
line(X(1:2,1),X(1:2,2));
line(X(2:3,1),X(2:3,2));
line(X([1,3],1),X([1,3],2));
%}
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
    
%     if nVar==2
%         % Plotting
%         figure(gcf);
%         hold on;
%         line([X(1,1),X(2,1)],[X(1,2),X(2,2)]);
%         line([X(1,1),X(3,1)],[X(1,2),X(3,2)]);
%         line([X(3,1),X(2,1)],[X(3,2),X(2,2)]);
%         plot(X(1,1),X(1,2),'bo','MarkerSize',8);
%         plot(X(2,1),X(2,2),'bo','MarkerSize',8);
%         plot(X(3,1),X(3,2),'bo','MarkerSize',8);
%         plot(M(1,1),M(1,2),'ro','MarkerSize',8);
%         hold off;
%         pause(0.001);
%     end
%     
    % FM=mean(FX(1:dim));
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
    % plot Adaptive Nelder-Mead simplex
    %pause(0.1)
    %clf;
    %plot(X(:,1),X(:,2),'k*');
    %line(X(1:2,1),X(1:2,2));
    %line(X(2:3,1),X(2:3,2));
    %line(X([1,3],1),X([1,3],2));
    %fprintf('Iteration = %4d\t :\t Best Cost = %5.4e\n', iter, BestCost(iter));
end
xmin=X(1,:);
fmin=FX(1);
nFeval=ct;