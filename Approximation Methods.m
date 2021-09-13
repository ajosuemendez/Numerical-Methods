%% Authors
%Alejandro Josue Mendez Lopez
%Hashim Alsayed

%% TASK 1
x = linspace(-2,2,100);
y = x.^2.*sin(x) + cos(4*x);

figure(1);
plot(x,y);
title("Original Function");
xlabel("x");
ylabel("y(x)");
legend("x^2sin(x) + cos(4x)");
 
%% TASK 2
N = 100; %Number of Samples
x_in = linspace(-2,2,N);
y_in = x_in.^2.*sin(x_in) + cos(4*x_in); %Given function
%% TASK 3
n_6 = 6; %Order
a_6 = approximation(x_in,y_in,n_6); %Coefficients for the 6th order p(x)

n_16 = 16; %Order
a_16 = approximation(x_in,y_in,n_16); %Coefficients for the 6th order p(x)

xd = linspace(-2,2,100);

yd_6  = DensePoints_yd(xd, a_6); %y values for the dense function of the 6th order
yd_16 = DensePoints_yd(xd, a_16);%y values for the dense function of the 16th order

figure(2);
plot(x, y, '-', x_in, y_in, 'o', xd, yd_6, xd, yd_16);
title("Approximation Method");
xlabel("x");
ylabel("y(x)");
legend("Given Function f(x)", "Input Data Points", "p(x) 6th order", "p(x) 16th order");

%{
Analisis:
The plot shows that the higher the order of the approximation function, the
more accurate is related to the given function. Nevertheless,using lower
orders give pretty much good results with a relative small error. As is
expected the p(x) functions try to go through the closest path in between
the input data points. The higher the p(x) order the higher the number of
coefficients to be calculated. In order to get high accurate results is
necessary to increase the order. This also increase the computational power needed.
 
It is worth mentioning that the higher the order of p(x) the more over feed
it will be. This produces sometimes inconsistencies when extrapolating.
Therefore it is convenient sometimes to stay with lower orders of p(x).
%}

%% TASK 4
e_6 = abs(y - yd_6);    %Error for the 6th p(x) related to the given function
e_16 = abs(y - yd_16);  %Error for the 16th p(x) related to the given function

fprintf('The Average Error for the P(x) 6th order is: %d\nThe Maximum Error for the P(x) 6th order is: %d\n',averError(e_6), maxError(e_6));
fprintf('The Average Error for the P(x) 16th order is: %d\nThe Maximum Error for the P(x) 16th is: %d\n',averError(e_16), maxError(e_16));

%{
Analisis:
As it was written above the higher the order of p(x) the more accurate the
results. Therefore we see that the average error and the maximal error is
considerably smaller for the 16th p(x) rather than the 6th p(x).
%}
%% TASK 5
%The following table represents the average errors and max errors for
%diferrent number of input data points
%{
   order      AverageError  MaxError
     5        5.685631e-01  1.090190e+00
     8        7.429732e-02  2.524183e-01
     11       1.144458e-02  4.026405e-02
     13       1.193051e-03  4.106611e-03
     15       9.148874e-05  2.914814e-04
     18       2.427540e-07  6.254774e-07
     21       1.227839e-07  5.079608e-07

Analisis:
The table shows the expected results. The higher the p(x) order the more
accurate it is in comparison with the given function. In both cases
(Average and Max) the error decreases significantly. Most likely from p(x)
orders equal or higher than 18 the errors are so small that there is no
need to keep increasing the order. So we can save some computational power.
%}
%% TASK 6
p = {@p0,@p1,@p2,@p3,@p4,@p5,@p6,@p7}; %Pointer of functions

A_legendre = zeros(length(p)); %Pre-initializing the Matrix A for storing values
B_legendre = zeros(length(p),1);%Pre-initializing the Matrix B for storing values

scale  = 2/(max(x_in) - min(x_in)); %Scale factor
offset = 1-scale*min(x_in); 
xps = scale*x_in + offset; %Adapting the input data points between -1 and 1

for r = 1:length(p)
    for c = 1:length(p)
        A_legendre(r,c) = sumuj(xps, p{r}, p{c}); %Allocating the values for each specific row and column in A MAtrix
    end
    B_legendre(r) = sumujb(xps, y_in, p{r}); %Allocating the values for each specific row and column in B Matrix
end
a_legendre = A_legendre\B_legendre;         %Calculating the Coefficients a

y_legendre = legendre(a_legendre, p, xd*scale + offset); %Calculating the y values for the corresponding least sqaure method function

figure(3);
plot(x,y,'-',xd, y_legendre,'--');
title("Least Square Method Approximation");
xlabel("x");
ylabel("y(x)");
legend("Original Function", "Approximation Function");
%{
It is seen that the plot is similar to the 6th order of p(x). The good
thing about using this method is that the number of samples can be
significantly reduced and still get pretty much good results. On the other
hand, this method is a little bit more complex to perform since here is
necessary to define the basis functions.
%}
%% Functions
function a = approximation(x_in,y_in,n)
A = zeros(n+1);
B = zeros(n+1,1);
N = length(x_in);

for r = 1:n+1
    for c = 1:n+1
        s = 0;
        for j=1:N
            s = s + x_in(j)^(r -1 + c-1);
        end
        A(r,c) = s;
    end
    k = 0;
    for i= 1:N
        k = k + y_in(i)*x_in(i)^(r-1);
    end
    B(r,1) = k;
end
a = A\B;
end

function yd = DensePoints_yd(xd, a)
yd = zeros(size(xd));
    for j = 1:length(xd)
            s = 0;
            for i = 1:length(a)
                s = s + a(i)*xd(j)^(i-1);
            end
            yd(j)= s;
    end
end

function result = averError(e) %Average Error
    N = length(e);
    result=0;
    for i=1:N-1
        result = result+ ((e(i) + e(i+1))/2);
    end
    result = result/(N-1);
end

function max = maxError(e) %Max Error
    N = length(e);
    max = 0;
    for j=1:N
        if e(j)>max
        max = e(j);
        end
    end
end

function y = p0(x)
    y = sin(x);
end

function y = p1(x)
    y = cos(x);
end

function y = p2(x)
    y = sin(4*x);
end

function y = p3(x)
    y = cos(4*x);
end

function y = p4(x)
    y = x.*sin(x);
end

function y = p5(x)
    y = x.*cos(x);
end

function y = p6(x)
    y = (x.^2).*sin(x);
end

function y = p7(x)
    y = (x.^2).*cos(x);
end

function s = sumuj(x_in, pa, pb)
    s = 0;
    for i = 1:length(x_in)
        s = s + pa(x_in(i))*pb(x_in(i));
    end
end

function s = sumujb(x_in, y_in, pa)
    s = 0;
    for j=1:length(x_in)
        s = s + pa(x_in(j))*y_in(j);
    end
end

function y = legendre(a, p, xd)
    y = zeros(size(xd));
    for k = 1:length(a)
        y = y + a(k)*p{k}(xd);
    end
end
