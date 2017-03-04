function yout = rbfregress(xin,yin,xout,knots,lambda)
    if isempty(lambda)
        myfun = @(lambda)DS(xin,yin,knots,lambda);
        lambda = fminbnd(myfun,0,1)
    end
    [n,p] = size(xin);
    m = size(xout,1);
    tps = @(r)r.^2.*log(r);
    D = pdist2(xin,knots);
    D(D==0) = realmin;
    X = [ones(n,1),tps(D)];
    b = (X'*X+lambda*eye(size(knots,1)+1))\(X'*yin);
    D = pdist2(xout,knots);
    D(D==0) = realmin;
    X = [ones(m,1),tps(D)];
    yout = X*b;
end

function err = DS(xin,yin,knots,lambda)
    n = length(yin);
    [Train, Test] = crossvalind('HoldOut', n, 0.1);
    yout = rbfregress(xin(Train,:),yin(Train),xin(Test,:),knots,lambda);
    err = mean((yin(Test)-yout).^2);
end