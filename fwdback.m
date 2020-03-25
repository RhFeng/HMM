function [alpha, beta, gamma, loglik, xi_summed, gamma2] = fwdback(init_state_distrib, ...
   transmat, obslik, varargin)
% FWDBACK Compute the posterior probs. in an HMM using the forwards backwards algo.
% toolbox from Kevin Murphy

if 0 % nargout >= 5
  warning('this now returns sum_t xi(i,j,t) not xi(i,j,t)')
end

if nargout >= 5, compute_xi = 1; else compute_xi = 0; end
if nargout >= 6, compute_gamma2 = 1; else compute_gamma2 = 0; end

[obslik2, mixmat, fwd_only, scaled, act, maximize, compute_xi, compute_gamma2] = ...
   process_options(varargin, ...
       'obslik2', [], 'mixmat', [], ...
       'fwd_only', 0, 'scaled', 1, 'act', [], 'maximize', 0, ...
                   'compute_xi', compute_xi, 'compute_gamma2', compute_gamma2);

[Q T] = size(obslik);

if isempty(obslik2)
 compute_gamma2 = 0;
end

if isempty(act)
 act = ones(1,T);
 transmat = { transmat } ;
end

scale = ones(1,T);


loglik = 0;

alpha = zeros(Q,T);
gamma = zeros(Q,T);
if compute_xi
 xi_summed = zeros(Q,Q);
else
 xi_summed = [];
end

%%%%%%%%% Forwards %%%%%%%%%%

t = 1;
alpha(:,1) = init_state_distrib(:) .* obslik(:,t);
if scaled
 %[alpha(:,t), scale(t)] = normaliseC(alpha(:,t));
 [alpha(:,t), scale(t)] = normalise(alpha(:,t));
end
%assert(approxeq(sum(alpha(:,t)),1))
for t=2:T
 %trans = transmat(:,:,act(t-1))';
 trans = transmat{act(t-1)};
 if maximize
   m = max_mult(trans', alpha(:,t-1));
   %A = repmat(alpha(:,t-1), [1 Q]);
   %m = max(trans .* A, [], 1);
 else
   m = trans' * alpha(:,t-1);
 end
 alpha(:,t) = m(:) .* obslik(:,t);
 if scaled
   %[alpha(:,t), scale(t)] = normaliseC(alpha(:,t));
   [alpha(:,t), scale(t)] = normalise(alpha(:,t));
 end
 if compute_xi & fwd_only  % useful for online EM
   %xi(:,:,t-1) = normaliseC((alpha(:,t-1) * obslik(:,t)') .* trans);
   xi_summed = xi_summed + normalise((alpha(:,t-1) * obslik(:,t)') .* trans);
 end
 %assert(approxeq(sum(alpha(:,t)),1))
end
if scaled
 if any(scale==0)
   loglik = -inf;
 else
   loglik = sum(log(scale));
 end
else
 loglik = log(sum(alpha(:,T)));
end

if fwd_only
 gamma = alpha;
 beta = [];
 gamma2 = [];
 return;
end

%%%%%%%%% Backwards %%%%%%%%%%

beta = zeros(Q,T);
if compute_gamma2
  if iscell(mixmat)
    M = size(mixmat{1},2);
  else
    M = size(mixmat, 2);
  end
 gamma2 = zeros(Q,M,T);
else
 gamma2 = [];
end

beta(:,T) = ones(Q,1);
%gamma(:,T) = normaliseC(alpha(:,T) .* beta(:,T));
gamma(:,T) = normalise(alpha(:,T) .* beta(:,T));
t=T;
if compute_gamma2
 denom = obslik(:,t) + (obslik(:,t)==0); % replace 0s with 1s before dividing
 if iscell(mixmat)
   gamma2(:,:,t) = obslik2(:,:,t) .* mixmat{t} .* repmat(gamma(:,t), [1 M]) ./ repmat(denom, [1 M]);
 else
   gamma2(:,:,t) = obslik2(:,:,t) .* mixmat .* repmat(gamma(:,t), [1 M]) ./ repmat(denom, [1 M]);
 end
 %gamma2(:,:,t) = normaliseC(obslik2(:,:,t) .* mixmat .* repmat(gamma(:,t), [1 M])); % wrong!
end
for t=T-1:-1:1
 b = beta(:,t+1) .* obslik(:,t+1);
 %trans = transmat(:,:,act(t));
 trans = transmat{act(t)};
 if maximize
   B = repmat(b(:)', Q, 1);
   beta(:,t) = max(trans .* B, [], 2);
 else
   beta(:,t) = trans * b;
 end
 if scaled
   %beta(:,t) = normaliseC(beta(:,t));
   beta(:,t) = normalise(beta(:,t));
 end
 %gamma(:,t) = normaliseC(alpha(:,t) .* beta(:,t));
 gamma(:,t) = normalise(alpha(:,t) .* beta(:,t));
 if compute_xi
   %xi(:,:,t) = normaliseC((trans .* (alpha(:,t) * b')));
   xi_summed = xi_summed + normalise((trans .* (alpha(:,t) * b')));
 end
 if compute_gamma2
   denom = obslik(:,t) + (obslik(:,t)==0); % replace 0s with 1s before dividing
   if iscell(mixmat)
     gamma2(:,:,t) = obslik2(:,:,t) .* mixmat{t} .* repmat(gamma(:,t), [1 M]) ./ repmat(denom,  [1 M]);
   else
     gamma2(:,:,t) = obslik2(:,:,t) .* mixmat .* repmat(gamma(:,t), [1 M]) ./ repmat(denom,  [1 M]);
   end
   %gamma2(:,:,t) = normaliseC(obslik2(:,:,t) .* mixmat .* repmat(gamma(:,t), [1 M]));
 end
end

