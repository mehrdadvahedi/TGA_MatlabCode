clc
clear all

for B=1:63
[D , L , U] = SetBounds(B);
GB = 0; % global best index
% D dimention of problem
N1 = 15;
N2 = 45;
N3 = 50;
N4 = 15;
N = N3+N2+N1; % Population of swarm
teta = 0.8;
maxIter = 2000;
T = unifrnd(L,U,[N,D]);
fit = zeros(1,N);

% calculate fit
for i=1:N
   fit(i) =  f(T(i,:) , B);
end
% set N1 , N2 , N3 tree groups ans sort T
[~, sortedIndexesN] = sort(fit);
temp = T;
for i=1:N
   T(i,:) = temp(sortedIndexesN(i),:);
end

for iter=1:maxIter
    % N1 trees
    for i=1:N1
       Tj = T;
       T(i,:) = (T(i,:)/teta)+rand .* T(i,:);
       if f(T(i,:) , B) > fit(i)
           T(i,:) = Tj(i,:);
       end
    end
    % N2 trees
    d = zeros((N2-N1),(N1+N2));
    for i=(N1+1):(N2)
       for k=1 : (N2)
            if i~=k
                d(i,k) = sqrt(sum((T(i,:)-T(k,:)).^2));
            else
                d(i,k) = inf;
            end
       end
       [tempSorted,tempIndex] = sort(d(i,:));
       x1 = T(tempIndex(1),:);
       x2 = T(tempIndex(2),:);
       landa = unifrnd(0,1);
       y = landa*x1 + (1-landa)*x2;
       alpha = unifrnd(0,1);
       T(i,:) = T(i,:) + alpha*y;
    end
    % N3 trees
    for i=(N1+N1+1):N
       T(i,:) =unifrnd(-10,10,[1,D]);
    end
    % N4 trees
    S = zeros(N4,D); % S(N,D)
    for k=1:N4
        S(k,:) = unifrnd(L,U,[1,D]);
        randIndex = randi(N4);
    % mask oprator
        for dim=1:D
           binRand = unifrnd(0,1);
           if binRand >0.5
               S(k,dim) = T(randIndex,dim);
           end
        end
    end
    % choos N best trees from T .+ S
%   Set fit and sort T
    tempT = zeros((N4+N),D); % S(N,D)
    tempFit = zeros(1,(N4+N));
    for i=1:N
        tempT(i,:) = T(i,:);
        tempFit(i) = f(tempT(i,:) , B);
    end
    for i=1:N4
       tempT(i+N,:) = S(i,:);
       tempFit(i+N) = f(tempT(i+N,:) , B);
    end
    [sortedTempFit, sortedIndexN4] = sort(tempFit);
    for i=1:N
       T(i,:) = tempT(sortedIndexN4(i),:);
       fit(i) = sortedTempFit(i);
    end
%     disp(['t = ' num2str(iter) ' , b = ' num2str(fit(1))]);
end
disp(['test = ' num2str(B) ' ans = ' num2str(fit(1)) ' D = ' num2str(D)]);
end