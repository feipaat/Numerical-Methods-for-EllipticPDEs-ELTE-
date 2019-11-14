function [x,itt] = dampedjacobi(N,b,x0,it)


%% Damped Jacobi iteration for multigrid; Poisson (The parameter it is the intersting one because of smoothing)
% Input list
% A       matrix
% b       column vector 
% xo      initial (starting) vector
% maxit   maximum iteration steps; stopping criterion
% omega   optimal value is 0.8 for this problem

itt=1;
X(:,1) = x0;
omega=2/3;

%% Splitting
e=ones(N+1,1);
A=spdiags([-e 2*e -e],[-1 0 1],N+1,N+1);
I=eye(length(A));
D=diag(diag(A));	
Romega=I-omega/2*A;

%% Damped Jacobi iteration with stopping criteria
while itt<it	% maxiterig 
X(:,itt+1) = Romega*X(:,itt)+omega*D\b;
itt = itt+1;
end
x= X(:, itt);
