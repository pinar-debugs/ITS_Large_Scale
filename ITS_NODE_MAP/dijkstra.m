function [costs,paths] = dijkstra(AorV,xyCorE,SID,FID,iswaitbar)
%% DIJKSTRA Calculate Minimum Costs and Paths using Dijkstra's Algorithm
%
%   Inputs:
%     [AorV] Either A or V where
%         A   is a NxN adjacency matrix, where A(I,J) is nonzero (=1)
%               if and only if an edge connects point I to point J
%               NOTE: Works for both symmetric and asymmetric A
%         V   is a Nx2 (or Nx3) matrix of x,y,(z) coordinates
%     [xyCorE] Either xy or C or E (or E3) where
%         xy  is a Nx2 (or Nx3) matrix of x,y,(z) coordinates (equivalent to V)
%               NOTE: only valid with A as the first input
%         C   is a NxN cost (perhaps distance) matrix, where C(I,J) contains
%               the value of the cost to move from point I to point J
%               NOTE: only valid with A as the first input
%         E   is a Px2 matrix containing a list of edge connections
%               NOTE: only valid with V as the first input
%         E3  is a Px3 matrix containing a list of edge connections in the
%               first two columns and edge weights in the third column
%               NOTE: only valid with V as the first input
%     [SID] (optional) 1xL vector of starting points
%         if unspecified, the algorithm will calculate the minimal path from
%         all N points to the finish point(s) (automatically sets SID = 1:N)
%     [FID] (optional) 1xM vector of finish points
%         if unspecified, the algorithm will calculate the minimal path from
%         the starting point(s) to all N points (automatically sets FID = 1:N)
%     [iswaitbar] (optional) a scalar logical that initializes a waitbar if nonzero 
%
%   Outputs:
%     [costs] is an LxM matrix of minimum cost values for the minimal paths
%     [paths] is an LxM cell array containing the shortest path arrays
%

%% check Inputs
    narginchk(2,5); 
    [n,nc] = size(AorV);
    [m,mc] = size(xyCorE);
    [E,cost, all_positive] = processInputs(AorV,xyCorE);
    % all_positive = 1;
    
    if nargin < 5,   iswaitbar = 0;  end
    if nargin < 4,   FID = (1:n);    end
    if nargin < 3,   SID = (1:n);    end
    if max(SID) > n || min(SID) < 1
        eval(['help ' mfilename]);
        error('Invalid [SID] input. See help notes above.');
    end
    if max(FID) > n || min(FID) < 1
        eval(['help ' mfilename]);
        error('Invalid [FID] input. See help notes above.');
    end

    isreversed = 0;
    if length(FID) < length(SID)
        E = E(:,[2 1]);
        cost = cost';
        tmp = SID;
        SID = FID;
        FID = tmp;
        isreversed = 1;
    end

    L = length(SID);
    M = length(FID);
    costs = zeros(L,M);
    paths = num2cell(nan(L,M));

    % Find the Minimum Costs and Paths using Dijkstra's Algorithm
    if iswaitbar, wbh = waitbar(0,'Please Wait ... '); end
    for k = 1:L
        % Initializations
        if all_positive, TBL = sparse(1,n); else TBL = NaN(1,n); end
        min_cost = Inf(1,n);
        settled = zeros(1,n);
        path = num2cell(nan(1,n));
        I = SID(k);
        min_cost(I) = 0;
        TBL(I) = 0;
        settled(I) = 1;
        path(I) = {I};

        while any(~settled(FID))
            % Update the Table
            TAB = TBL;
            if all_positive, TBL(I) = 0; else TBL(I) = NaN; end
            nids = find(E(:,1) == I);
            % Calculate the Costs to the Neighbor Points and Record Paths
            for kk = 1:length(nids)
                J = E(nids(kk),2);
                if ~settled(J)
                    c = cost(I,J);
                    if all_positive, empty = ~TAB(J); else empty = isnan(TAB(J)); end
                    if empty || (TAB(J) > (TAB(I) + c))
                        TBL(J) = TAB(I) + c;
                        if isreversed
                            path{J} = [J path{I}];
                        else
                            path{J} = [path{I} J];
                        end
                    else
                        TBL(J) = TAB(J);
                    end
                end
            end

            if all_positive, K = find(TBL); else K = find(~isnan(TBL)); end
            % Find the Minimum Value in the Table
            N = find(TBL(K) == min(TBL(K)));
            if isempty(N)
                break
            else
                % Settle the Minimum Value
                I = K(N(1));
                min_cost(I) = TBL(I);
                settled(I) = 1;
            end
        end
        % Store Costs and Paths
        costs(k,:) = min_cost(FID);
        paths(k,:) = path(FID);
        if iswaitbar, waitbar(k/L,wbh); end
    end
    if iswaitbar, close(wbh); end

    if isreversed
        costs = costs';
        paths = paths';
    end

    if L == 1 && M == 1
        paths = paths{1};
    end

end
