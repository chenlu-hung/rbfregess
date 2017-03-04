A = dlmread('motorcycle.txt');
x = A(:,1); y = A(:,2);
xout = linspace(min(x),max(x),101)';
f = llinear_openmp(x',y',xout',4);
knots = quantile(x,11)';
g = rbfregress(x,y,xout,knots,0);
scatter(x,y)
hold, plot(xout,f,'b')
plot(xout,g,'r')