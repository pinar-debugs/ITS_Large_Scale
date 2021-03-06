
close all
clear
clc

%% generate random n points
      % Calculate the shortest distance and path between two points
      n = 1000;      
      xy = 10*rand(n,2); %generate random n points
%% define center and clear points close to it to create single high density point
      center = ones(size(xy));
      center = center*5;
      K = sqrt( sum( (xy-center).^2, 2) ); %sum(X,DIM) sums along the dimension DIM.
      xy(K <= 0.5, :) = [];
      xy = [xy;[5,5]];
      n=size(xy,1);
      A = zeros(n); 

%% triangulate points and store point to point relations
      tri = delaunay(xy(:,1),xy(:,2)); %triangulate n points
      %tri's format: [1 3 9] means point1, p3, p9 are corner of triangle
      I = tri(:); %format as sub-bottom, makes them single column matrix
      J = tri(:,[2 3 1]);  % change the order of columns, shift first column to last
      J = J(:); %format as sub-bottom, makes them single column matrix
      % until this time, we had tri as nX3, and then we generated I and J 
      % tri stores triangular relation, we want to make point2point
      % relation using I and J which was generated by formating tri

%% delete long ways to generate more complex map       
      D = sqrt( sum( (xy(I,:)-xy(J,:)).^2, 2) ); % sum squares of columns in a row and then sqrt them
      % calculates distancce btwn points
      I(D > 0.75,:) = []; %if distance greater than 0.75 mark as []
      J(D > 0.75,:) = []; %if distance greater than 0.75 mark as []

%% adjacency and drawing       
      IJ = I + n*(J-1); % generates random adjacency,
      A(IJ) = 1;
      
      gplot(A,xy,'k.:'); hold on;
      drawCircle(5, 5, 1);hold on;% draw map center

%% dijkstra : calculate costs and paths
      [cost,path] = dijkstra(A,xy,1,n)
      xy(path,:) %coordinates of points
            
      plot(xy(path,1),xy(path,2),'ro-','LineWidth',2); hold off;
      sprintf('end');