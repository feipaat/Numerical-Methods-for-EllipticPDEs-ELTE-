function [x,it] = jacobi(A,b,x0,TOL,maxit)

%% Jacobi-iteration with stopping criteria (for other criteria see: Farago-Horvath Section 3.6.5. in Hungarian)
% Input list
% A       matrix
% b       column vector 
% xo      initial (starting) vector
% TOL     given tolerance; relative error of numerical solution vector in maximum norm; stopping criterion
% maxit   maximum iteration steps; stopping criterion


relerr = inf; %At the beginning we give a huge relative error
it = 1;
X(:,1) = x0;

%% Splitting
D=diag(diag(A));	
L=tril(sparse(A),-1);
U=triu(sparse(A),1);

%% Jacobi iteration with stopping criteria
while (relerr>TOL) && (it<maxit) 	% Iterating until maxiter or the required error
X(:,it+1) = D\(-(L+U)*X(:,it)+b);
relerr = norm(X(:,it+1)-X(:,it),inf)/norm(X(:,it),inf);
it = it+1;
end
x= X(:, it);
