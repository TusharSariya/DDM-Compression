% Copywrite 2008, Maria Paola Clarizia
% license: GPL, see gpl.txt

function y=lambda(x)
%		y=lambda(x) generates a triangular signal
y=(x+1).*(x>-1 & x<0)+(x==0)+(-x+1).*(x>0 & x<1);