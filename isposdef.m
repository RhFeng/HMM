function b = isposdef(a)

%  From Tom Minka's lightspeed toolbox

[R,p] = chol(a);
b = (p == 0);
