function [E,C, all_positive] = processInputs(AorV,xyCorE)

all_positive = 1;
[n,nc] = size(AorV);
[m,mc] = size(xyCorE);
C = sparse(n,n);
        if n == nc %check if square matrix
            if m == n %check row count of AorV and xyCorE, are they equal
                if m == mc % Inputs are A, COST matrix
                    A = AorV;
                    A = A - diag(diag(A)); 
                    %  diag(X) is the main diagonal of X. diag(diag(X)) is a diagonal matrix.
                    C = xyCorE;
                    all_positive = all(C(logical(A)) > 0);
                    E = a2e(A); %convert adjacency matrix 2 edge list
                else % Inputs are A, xy adjacency matrix
                    A = AorV;
                    A = A - diag(diag(A)); %this makes diagonal of A equal to 0, all diagonal members are zero
                    %  diag(X) is the main diagonal of X. diag(diag(X)) is a diagonal matrix.
                    xy = xyCorE;
                    E = a2e(A); %convert adjacency matrix 2 edge list
                    D = ve2d(xy,E); %compute euclidean distance for edges 
                    all_positive = all(D > 0);
                    for row = 1:length(D)
                        C(E(row,1),E(row,2)) = D(row);
                    end
                end
            else
                eval(['help ' mfilename]);
                error('Invalid [A,xy] or [A,cost] inputs. See help notes above.');
            end
        else
            if mc == 2 % Inputs: V,E
                V = AorV;
                E = xyCorE;
                D = ve2d(V,E);
                all_positive = all(D > 0);
                for row = 1:m
                    C(E(row,1),E(row,2)) = D(row);
                end
            elseif mc == 3 % Inputs: V,E3
                E3 = xyCorE;
                all_positive = all(E3 > 0);
                E = E3(:,1:2);
                for row = 1:m
                    C(E3(row,1),E3(row,2)) = E3(row,3);
                end
            else
                eval(['help ' mfilename]);
                error('Invalid [V,E] inputs. See help notes above.');
            end
        end
    end
