function [X,fX,clrCode] = wAdaptiveNelderMead(Cost,X,fX)
    global calpha cbeta cgamma1 cgamma2 cdelta ctear iter BestCost VarMin VarMax op
    dim = size(X,1)-1;
    L = 8;

    % Centroid of the best vertices
    w = (2*(dim-(1:dim))+1)/(dim^2);
    M = mean(w.*X(1:dim,:));

    alpha = 1;
    gamma1 = 0.75 + 0.1/(1+cgamma1);
    gamma2 = 0.5 - 0.1/(1+cgamma2);
    delta = 0.5 + 0.1/(1+cdelta);
    term = (dim^2 - w.*(dim-w))./(w.*(dim-w));
    if iscolumn(term)
        beta =  term';
    else
        beta = term;
    end

    % Reflection point
    R = reflection(alpha,M,X(dim+1,:));
    fR = Cost(R);
    if fR < fX(1)
        % Update parameter for expansion point
        E = expansion(beta,M,R);
        fE = Cost(E);
        if fE < fR
            X(dim+1,:) = E;
            fX(dim+1) = fE;
            op = 'Expansion'; flag = 2;
            cbeta = cbeta + 1;
        else
            X(dim+1,:) = R;
            fX(dim+1) = fR;
            op = 'Reflection'; flag = 1;
            calpha= calpha + 1;
        end
    else
        if fR < fX(dim)
            % accept reflection point
            X(dim+1,:) = R;
            fX(dim+1) = fR;
            op = 'Reflection'; flag = 1;
            calpha = calpha + 1;
        else
            if fR < fX(dim+1)
                % Update parameter for Outside Contraction point
                OC = outsideContraction(gamma1,M,R);
                fOC = Cost(OC);
                if fOC <= fR
                    X(dim+1,:) = OC;
                    fX(dim+1) = fOC;
                    op = 'Outside Contraction'; flag = 3;
                    cgamma1 = cgamma1 + 1;
                else
                    % Shrink
                    for i=2:dim+1
                        % Update parameter for shrink point
                        S = shrink(delta,X(1,:),X(i,:));
                        fS = Cost(S);
                        X(i,:) = S;
                        fX(i) = fS;
                    end
                    op = 'Shrink'; flag = 5;
                    cdelta = cdelta+1;
                end
            else
                % Update parameter for Inside Contraction point
                IC = insideContraction(gamma2,M,R);
                fIC = Cost(IC);
                if fIC < fX(dim+1)
                    X(dim+1,:) = IC;
                    fX(dim+1) = fIC;
                    op = 'Inside Contraction'; flag = 4;
                    cgamma2 = cgamma2 + 1;
                else
                    % Shrink
                    for i=2:dim+1
                        % Update parameter for shrink point
                        S = shrink(delta,X(1,:),X(i,:));
                        fS = Cost(S);
                        X(i,:) = S;
                        fX(i) = fS;
                    end
                    op = 'Shrink'; flag = 5;
                    cdelta = cdelta+1;
                end
            end
        end
    end
    [fX,I] = sort(fX);
    X = X(I,:);
    BestCost(iter) = fX(1);
    
    if (iter>L) && (BestCost(iter)== BestCost(iter-L))
        X(dim+1,:) = VarMin + (VarMax - VarMin)*rand(1,dim);
        fX(dim+1) = Cost(X(dim+1,:));
        op = 'Tearing the worst vertex';
        ctear = ctear + 1;
    end
    
    [fX, I] = sort(fX);
    X = X(I,:);
    switch flag
        case 1
            clrCode = 'magenta';
        case 2
            clrCode = 'blue';
        case 3
            clrCode = 'yellow';
        case 4
            clrCode = 'cyan';
        case 5
            clrCode = 'red';
    end
end