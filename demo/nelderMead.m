function [X,FX,clrCode] = nelderMead(Cost,X,FX)

global calpha cbeta cgamma1 cgamma2 cdelta iter BestCost VarMin VarMax op

alpha = 1; beta = 2; gamma = 0.5; delta = 0.5;
dim = size(X,1)-1;

M = mean(X(1:dim,:));

xref = (1+alpha)*M - alpha*X(dim+1,:);
Fref = Cost(xref);
if Fref<FX(1)
    % expansion
    xexp = (1+alpha*beta)*M - alpha*beta*X(dim+1,:);
    Fexp = Cost(xexp);
    if Fexp < Fref
        X(dim+1,:) = xexp;
        FX(dim+1) = Fexp;
        op = 'Expansion';
        flag=2;
        cbeta = cbeta + 1;
    else
        X(dim+1,:) = xref;
        FX(dim+1) = Fref;
        op = 'Reflection';
        flag=1;
        calpha = calpha + 1;
    end
else
    if Fref < FX(dim)
        % accept reflection point
        X(dim+1,:) = xref;
        FX(dim+1) = Fref;
        op = 'Reflection';
        flag=1;
        calpha = calpha + 1;
    else
        if Fref < FX(dim+1)
            % Outside contraction
            xoc = (1+alpha*gamma)*M - alpha*gamma*X(dim+1,:);
            Foc = Cost(xoc);
            if Foc <= Fref
                X(dim+1,:) = xoc;
                FX(dim+1) = Foc;
                op = 'Outside Contraction';
                flag = 3;
                cgamma1 = cgamma1 + 1;
            else
                % shrink
                for i=2:dim+1
                    X(i,:) = X(1,:) + delta*(X(i,:)-X(1,:));
                    FX(i) = Cost(X(i,:));

                end
                op = 'Shrink';
                flag = 5;
                cdelta = cdelta + 1;
            end
        else
            %inside contraction
            xic = (1-gamma)*M + gamma*X(dim+1,:);
            Fic = Cost(xic);
            if Fic < FX(dim+1)
                X(dim+1,:) = xic;
                FX(dim+1) = Fic;
                op = 'Inside contraction';
                flag = 4;
                cgamma2 = cgamma2 + 1;
            else
                % shrink
                for i=2:dim+1
                    X(i,:) = X(1,:) + delta*(X(i,:)-X(1,:));
                    FX(i) = Cost(X(i,:));
                end
                op ='Shrink';
                flag = 5;
                cdelta = cdelta + 1;
            end
        end
    end
end

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
[FX, I] = sort(FX);
X = X(I,:);
end