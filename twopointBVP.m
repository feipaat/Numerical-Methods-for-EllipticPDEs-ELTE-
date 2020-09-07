function [y_h]=twopointBVP(inta,intb,alpha,beta,N)
%% Numerical solution for two-point BVPs
%
%     a(t)u''(t)+c(t)u(t)+d(t)u'(t)=f(t) c\in\mathbb{R}, f(t)\inC[a,b]
%     u(a)=\alpha u(b)=\beta
%

%% Input:

%     inta Start of the interval
%     intb End of the interval
%     N    Number of intervals


%% Preliminaries

% Step size

h=(intb-inta)/N;

% Forming the discretization matrix

%Use index j instead of i because of the imaginary unit meaning in MATLAB/Octave
for j=1:N-1
    a(j)=a1(inta+j*h); %a_i
end

for j=1:N-1 
    c(j)=c1(inta+j*h); %c_i
end

for j=1:N-1 
    d(j)=d1(inta+j*h); %d_i
end
e=ones(N-1,1);
A_h=(1/h^2)*spdiags([a'-0.5*h*d' -2*a'+h^2*c' a'+0.5*h*d'],-1:1,N-1,N-1);

%% Calculating the numerical solution and plotting it

b_h=zeros(N-1,1);
b_h(1)=f(inta+h)-alpha*(a(1)/h^2-d(1)/(2*h));
b_h(N-1)=f(inta+(N-1)*h)-beta*(a(N-1)/h^2+d(N-1)/(2*h));
for i=2:N-2 
    b_h(i)=f(inta+i*h);
end
y=A_h\b_h;
y_i=linspace(alpha,beta,N+1);
y_i(2:N)=y;
y_h=y_i';

x_i=inta:h:intb;
plot(x_i,y_i,'r+')
hold on;

%% The problem's right hand side function  f
function rhs_f=f(t)
rhs_f=0;
%% The problem's left hand side function a(t) 
function lhs_a=a1(t)
lhs_a=1;
%% The problem's left hand side function c(t) 
function lhs_c=c1(t)
lhs_c=-1;
%% The problem's left hand side function d(t) 
function lhs_d=d1(t)
lhs_d=0;
