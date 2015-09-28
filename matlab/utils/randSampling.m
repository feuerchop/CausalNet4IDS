function [ds_tr, ds_t] = randSampling( org_d, tr_s, t_s)
%RANDSAMPLING Summary of this function goes here
%   Detailed explanation goes here
   dsize = size(org_d, 1);
   if (nargin > 2)
      % we want also sample test dataset
      if (tr_s + t_s > dsize)
         fprintf('SAMPLE SIZE SHOULD NOT EXCEED %d!\n', dsize);
         return;
      else
         tr_idx = randsample(dsize, tr_s);
         t_idx = randsample(mysetdiff(1:dsize, tr_idx), t_s);
         ds_tr = org_d(tr_idx, :);
         ds_t = org_d(t_idx, :);
      end
   else
      if (tr_s > dsize)
         fprintf('TRAINING SIZE SHOULD NOT EXCEED %d!\n', dsize);
         return;
      end
      if (nargout > 1)
         fprintf('YOU SHOULD ALSO INPUT TESTSET SIZE!\n');
         return;
      end
      tr_idx = randsample(dsize, tr_s);
      ds_tr = org_d(tr_idx, :);
   end
end

