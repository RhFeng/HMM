function z = calcov2(mu, Sigma, varargin)

% toolbox from Mark A. Paskin

if size(Sigma) ~= [2 2], error('Sigma must be a 2 by 2 matrix'); end
if length(mu) ~= 2, error('mu must be a 2 by 1 vector'); end


p=0.9;
n=100;
h = [];
holding = ishold;
if (Sigma == zeros(2, 2))
  z = mu;
else
  % Compute the Mahalanobis radius of the ellipsoid that encloses
  % the desired probability mass.
  k = chi2inv(p, 2);
  % The major and minor axes of the covariance ellipse are given by
  % the eigenvectors of the covariance matrix.  Their lengths (for
  % the ellipse with unit Mahalanobis radius) are given by the
  % square roots of the corresponding eigenvalues.
  if (issparse(Sigma))
    [V, D] = eigs(Sigma);
  else
    [V, D] = eig(Sigma);
  end
  % Compute the points on the surface of the ellipse.
  t = linspace(0, 2*pi, n);
  u = [cos(t); sin(t)];
  w = (k * V * sqrt(D)) * u;
  z = repmat(mu, [1 n]) + w;

end

% cmap=hsv(lithos_num);
% 
% % h = plot(z(1, :), z(2, :), plot_opts{:},'color',cmap(mod(lithos_index,lithos_num)+1,:),'LineWidth',4);
% 
% h=plot(z(1, :), z(2, :),'color',cmap(mod(lithos_index,lithos_num)+1,:),'LineWidth',4);
% xlim([X(1) X(2)]);
% ylim([Y(1) Y(2)]);
% grid on;
