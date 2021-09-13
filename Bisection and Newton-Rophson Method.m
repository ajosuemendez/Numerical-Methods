%% Task 1 Nonlinear Equation
x_original = linspace(-10,10,2000);
y_original = 0.05*x_original.^3 -50 +0.01*exp(x_original);

figure(1);
plot(x_original,y_original);
title('Polynomial Function To Be Analysed');
xlabel('x');
ylabel('y');
xlim([-10 10]);
%{
The main goal of this report is to implement and analyse both methods
(Bisection and Newton-Rophson) which allow us to determine roots values
from a non-linear equation.
%}
%% Task 2 Bisection Method
%First we display the original function
figure(2);
plot(x_original,y_original);
title('Bisection Method');
xlabel('x');
ylabel('y');
xlim([-10 10]);

%We need to define to boundaries where the bisection algorithm will take place
lower = -10;
upper = 10;
[a,b,i] = bisection(lower,upper); %We call the function for the bisection method

%{
Now we display the calculated values from the Bisection Method.
It's good practice to determine the minimum value from both of them (a and b) in
order to get the most accurate result. However, after testing the program
both of them seem to be quite the same. therefore we only display the
obtained value of b (bisection root) and the number of iterations needed.
%}
fprintf("\nBisection Root found for x = %d",b); 
fprintf("\nBisection Root value = %d",f(b));
fprintf("\nBisection Number of iterations = %d",i);
%{
In general this method is quite easy to implement since the idea behind is
trivial. It gives good result with a relatively small computational power.
Nevertheless if we want to optimize the efficiency it is necessesary to
implement the Newton-Rophson Method.
%}

%% Task 3 Newton-Rophson Method
x = -10; %Starting Point
allowError = 1e-7; %chosen Allowed Error according the task
dx = 1e-1; %chosen differential size for the derivative

[x,k] = Rophson(x,allowError,dx); %We call the Rophson-Method Function

figure(3);
plot(x_original,y_original,'-',x,f(x),'o');
title('Newton-Rophson Method');
xlabel('x');
ylabel('y');
legend('Given Function','Found Root');
xlim([-10 10]);

fprintf("\n\nBisection Root found for x = %d",x);
fprintf("\nBisection Root value = %d",f(x));
fprintf("\nBisection Number of iterations = %d\n",k);
%{
In general this method tends to be more precise and is has a better performance in terms of computational efforts.
Using this algorithm we can decrease significantly the number of iterations
needed. However, choosing the right approach often depends on the function or equation under
investigation.
%}


%% Task 4 Influence of Allowable Error and Numerical dx
x = -10; %Starting Point
vec_AllowError = [1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8,1e-9,1e-10,1e-11,1e-12,1e-13,1e-14];
vec_dx = [1e-1,1e-2,1e-3,1e-4,1e-5];

Matrix = zeros(length(vec_dx),length(vec_AllowError));
for r=1:length(vec_dx)
    for c =1:length(vec_AllowError)
        [x,k] = Rophson(x,vec_AllowError(c),vec_dx(r));
        Matrix(r,c) = k;
        x=-10;
    end
end

for n=4:8
    figure(n);
    bar(Matrix(n-3,:));
    xlabel('Allowable error in 10 to the power of minus');
    ylabel('Number Of Iterations');
    if n==4 
        legend('dx=1e-1');
    elseif n==5
        legend('dx=1e-2');
    elseif n==6
        legend('dx=1e-3');
    elseif n==7
        legend('dx=1e-4');
    else
        legend('dx=1e-5');
    end
end
%{
The results show a clear trend that the smaller the allowable error the more iterations
are needed in order to find the root value. This makes sense since the
margin of error is reduced and therefore more trials are needed. We can
also see that the smaller the differential dx the lower the number of
iterations. By making the dx smaller we are derivating more precisely
which helps to reduce the number of iterations needed.
%}

%% Functions

function [x,k] = Rophson(x,allowError,dx) 
    for k=1:100
        x = x - f(x)/df(x,dx);
        if (abs(f(x))<allowError)
            break;
        end
    end
end

function dy = df(x,dx)
    dy = (f(x+dx) - f(x))/dx;
end


function [a,b,i] = bisection(a,b)
    hold on
    for i=1:100
        c = (a+b)/2;
        if f(a)*f(c)<0
            b=c;
        else
            a=c;
        end
        figure(2);
        plot(c,f(c),'o');
        if ((abs(f(a))<1e-7) || (abs(f(b))<1e-7))
            break;
        end
    end
    hold off
end

function y = f(x)
    y = 0.05*x.^3 -50 +0.01*exp(x);
end
