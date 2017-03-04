clear
n = 1000;
t = (linspace(0,1,n)');
[t1,t2] = ndgrid(t);
ytrue = exp(sqrt(0.5*(t1-0.6).^2+0.8*t2));
ytrue = ytrue(:);
yin = ytrue + normrnd(0,0.1,n^2,1);
xin = [t1(:),t2(:)];
[to1,to2] = ndgrid(linspace(0,1,n));

tt = (linspace(0,1,10)');
[tt1,tt2] = ndgrid(tt);
knots = [tt1(:),tt2(:)];
tic; yout = rbfregress(xin,yin,xin,knots,0); toc
mean((yout-ytrue).^2)
F = griddedInterpolant(to1,to2,reshape(yout,n,n));
yout = F(t1,t2);

tic; yhat = llinear_openmp(xin',yin',xin',[0.1,0.1]); toc
mean((yhat-ytrue).^2)
F = griddedInterpolant(to1,to2,reshape(yhat,n,n));
yhat = F(t1,t2);

mesh(t1,t2,reshape(ytrue,n,n))
figure, mesh(t1,t2,yout)
figure, mesh(t1,t2,yhat)