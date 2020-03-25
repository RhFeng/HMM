function [converged, decrease] = em_converged(loglik, previous_loglik, threshold, check_increased)
% EM_CONVERGED

% toolbox from Kevin Murphy

if nargin < 3, threshold = 1e-4; end
if nargin < 4, check_increased = 1; end

converged = 0;
decrease = 0;

if check_increased
  if loglik - previous_loglik < -1e-3 % allow for a little imprecision
    fprintf(1, '******likelihood decreased from %6.4f to %6.4f!\n', previous_loglik, loglik);
    decrease = 1;
converged = 0;
return;
  end
end

delta_loglik = abs(loglik - previous_loglik);
avg_loglik = (abs(loglik) + abs(previous_loglik) + eps)/2;
if (delta_loglik / avg_loglik) < threshold, converged = 1; end
