%%Report
%% Taks 1
len =100; %Number of samples for the original function
x = linspace(-5,5,len); %Sampling x interval
y = x.^3 -2*x.^2 + x -5; %Original Function

figure(1);
plot(x,y);
title('Polynomial Function');
xlabel('x');
ylabel('y');
legend('x^3 -2x^2 + x -5');

%% Task 2
n = 3;  %Number of input points (samples)
xin = linspace(-5,5,n); %Input points x
yin = xin.^3 -2*xin.^2 + xin -5; %Input points y
%% Task 3
xx = linspace(-5,5,len); %We use same amount of samples as the orginal function
yy = zeros(length(xx),1);

for t=1:length(xx)
    yy(t) = lagr(xin,yin, xx(t)); %Interpolation p(x)
end

figure(2);
plot(x,y,'--',xin,yin,'o', xx,yy,'-');
title("Lagrange Interpolation");
xlabel("x");
ylabel("y");
legend("Original Function", "Input Points","Interpolation p(x)");

%{
using only 3 sample points is very inaccurate.(we see a 2nd order function)
The closer we get to the input points the more accurate it becomes. 
However, it is necessary to increase the number of input samples to 6 to obtain satisfactory results.
%}

%% Task 4
e = abs(y' - yy); %Calculating the error, y is the original function and yy the interpolated p(x)

%Average and Max error can be calculated using the built-in Matlab func
%averError = sum(e)/length(e);
%maxError = max(e);
%However, for this lab separated functions were developped

fprintf('The Average Error is: %d\nThe Maximum Error is: %d',averError(e), maxError(e));

%% Task 5
%The following table represents the average errors and max errors for
%diferrent number of input data points
%{
   input      AverageError  MaxError
     5        5.192106e-15  5.684342e-14
     6        4.744390e-15  2.842171e-14
     7        7.150059e-15  8.526513e-14
     8        7.966548e-15  8.526513e-14
     9        6.537693e-15  8.526513e-14
     10       7.709495e-15  8.526513e-14
     11       7.756986e-15  1.136868e-13

In general the error is very small (both the maximum and the average).
Theoretically, if we increase the number of sample points we should see that the error is reduced.
However, the error is so minimal that matlab is not able to recognize it. 
It is very likely that there is an error when calculating it if we do not use variables with adequate capacities to store the result.
%}

%% Functions

function p = lagr(xin,yin,xx)
    p=0;
    for i=1:length(xin)
        L=1;
        for k=1:length(xin)
            if k ~=i
                L = L * (xx-xin(k))/(xin(i)-xin(k));
            end
        end
        p = p + yin(i)*L;
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

