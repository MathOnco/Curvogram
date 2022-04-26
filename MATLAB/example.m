clc;clear;close all;

% bins of histogram (also used for curvogram)
X = 0.5:0.2:2.5;

% normal distribution:
mu = 1.5; sigma = 0.5;
Y = 1./(sigma*sqrt(2*pi)).*exp((-1/2).*((X-mu)/sigma).^2);


%% example 1: f(x) = 0 (no transformation)
f = @(x) 0;
fprime = @(x) 0;

figure(1);
curvogram(X,Y,f,fprime);

%% example 2: sin(x)
f = @(x) sin(x);
fprime = @(x) cos(x);

figure(2);
curvogram(X,Y,f,fprime,'ScaleHeight',1,'BinWidth',1,'Color',blue());


%% example 3: cos(x)
f = @(x) cos(x);
fprime = @(x) -sin(x);

figure(3);
curvogram(X,Y,f,fprime,'ScaleHeight',1,'BinWidth',1);


%% example 4: x.^2
f = @(x) x.^2;
fprime = @(x) 2*x;

figure(4);
curvogram(X,Y,f,fprime,'ScaleHeight',15,'BinWidth',1);



