n = 1000000;
x = sort(unifrnd(0,4,n,1));
f = zeros(n,1); f(mod(floor(x),2)==1) = 1;
y = f + normrnd(0,.5,n,1);

xout = linspace(0,4,10000)';
knots = quantile(x,100)';
yy = interp1(x,y,xout,'linear','extrap');
mu = rbftv(xout,yy,xout,knots,1);

scatter(xout,yy)
hold on;
plot(xout,interp1(x,f,xout,'linear','extrap'),'k')
plot(xout,mu)