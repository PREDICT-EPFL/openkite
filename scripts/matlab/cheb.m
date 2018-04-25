function[x, Dn] = cheb(N)
%compute Chebyshev interpolation points and corresponding

%proper implementation
x = cos(pi*(0:N)/N)';
x = (x + 1) / 2;
c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
X = repmat(x,1,N+1);
dX = X-X';
Dn  = (c*(1./c)')./(dX+(eye(N+1)));        % off-diagonal entries
Dn  = Dn - diag(sum(Dn'));                 % diagonal entries

diag(sum(Dn'))

end

%helper function
function[ret] = c(i,N)
if(isequal(i,1) || isequal(i,N+1))
    ret = 2;
else
    ret = 1;
end
end