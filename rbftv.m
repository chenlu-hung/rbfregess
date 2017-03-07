function [yout] = rbftv(xin,yin,xout,knots,lambda)
    [n,p] = size(xin);
    e = ones(n,1);
    % Diff = spdiags([e -2*e e], 0:2, n-2, n);
    Diff = spdiags([e -e], 0:1, n-1, n);
    m = size(xout,1);
    D = pdist2(xin,knots);
    X = [ones(n,1),tps(D)];
    cvx_begin
        variable b(size(knots,1)+1)
        minimize norm(yin-X*b,2) + lambda*norm(Diff*X*b,1)
    cvx_end
    D = pdist2(xout,knots);
    X = [ones(m,1),tps(D)];
    yout = X*b;
end

function z = tps(r)
    z = r.^2.*log(r);
    z(r==0) = 0;
end