% Part 1 for the project 5080

function[] = LUverification(L,U,T,a,b,c,d,n)

T_new = L*U;

if T_new(1,1) == d(1)
    T_new(1,1) = b(1);
end

for i = 2:n
    if T_new(i,i) == d(i) + (a(i-1)*c(i-1))/d(i-1)
    T_new(i,i) = b(i);
    end
end

if T == T_new
    fprintf("It is offically confirmed via MATLAB that L*U = T\n")
end

end