% Finding solutions using LU Factorizations

function[x_LUFactorization] = LUFactorization_Solution(L,U,b)

y = inv(L)*b;

x_LUFactorization = inv(U)*y;


end