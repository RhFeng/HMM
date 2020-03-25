function y = logdet(A)

%  From Tom Minka's lightspeed toolbox

U = chol(A);
y = 2*sum(log(diag(U)));
