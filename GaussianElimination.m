A =[1 3 4;-2 2 3;1 1 2];
B =[2 -1 3];
B = B';

AUG = gaussian(A,B);
x = back(AUG(:,1:end-1),AUG(:,end));

%A \ B - x % the '\' implement the built-in elimination function in matlab
%so if we rest our calculated x values it should return 0 since 

function Ag = gaussian(A,B)
Ag =[A B];
n = size(Ag,1);
    for c=1:n-1
        for r = c+1:n
            l = Ag(r,c) / Ag(c,c);
            for j = c:n+1
                Ag(r,j) = Ag(r,j)- l*Ag(c,j);
            end
        end
    end
end

function a = back(A,B)
new_n = size(A,1);
a = zeros(size(A,1),1);
    for k = new_n:-1:1
        s=0;
        for j=k+1:new_n
            s = s +A(k,j)*a(j);
        end
        a(k) = (B(k)-s)/A(k,k);
    end
end        
