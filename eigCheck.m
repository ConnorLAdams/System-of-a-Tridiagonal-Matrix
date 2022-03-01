% eigenvalue check for project 5080

function[] = eigCheck(eig_H,eig_cosine,n)
true = 0;
for i = 1:n
   for j = 1:n
      if eig_H(i) == eig_cosine(j)
         true = 1;
         break
      end
   end
   if true ~= 1
       tru(i) = 0;
   else
       tru(i) = 1;
   end
   tru = tru';
end

if sum(tru) == n
   fprintf("It is confirmed that our eigenvalues are equal to the equation.\n") 
end

end