function d = DistOfBNs( bna, bnb )
%DISTOFBNS Summary of this function goes here
%   compute the hemming distance of two BNs among the same set of
%   variables
     if size(bna, 1) ~= size(bnb, 1)
         % different variables, we can not compare them
         d = -1;
         return;
     end
     if size(bna, 1) == 1
         d = 0;
     else
         missing = 0;
         extra = 0;
         inversion = 0;
         for i = 2:size(bna, 1)
             if bna(1, i) == 0 && bna(i, 1) == 0
                 % no such an edge in original
                 % but there is one in second bn
                 if bnb(1, i) == 1 || bnb(i, 1) == 1
                    extra = extra + 1;
                 end
             else
                 % there is an edge in orginal
                 % but there is none in second bn
                 if bnb(1, i) == 0 && bnb(i, 1) == 0
                     missing = missing + 1;
                 else
                     % two entries in two graphs should be exactly the same
                     % unless, inversion is invoked
                     if bna(1, i) ~= bnb(1, i)
                         inversion = inversion + 1;
                     end
                 end
             end
         end
         d = missing + extra + inversion + DistOfBNs(bna(2:size(bna,1), 2:size(bna,1)), bnb(2:size(bnb,1), 2:size(bnb,1)));
                        
end

