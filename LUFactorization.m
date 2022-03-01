% L and U factorization for Project 5080

function[L,U] = LUFactorization(a,d,c_over_d,n)

L = diag(ones(1,n)) + diag(c_over_d,-1);
U = diag(d) + diag(a,1);

end