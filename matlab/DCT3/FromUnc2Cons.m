function f = FromUnc2Cons(x, a, flag)

% Let x be an unconstrained parameter. This function takes x and
% transforms it to a parameter f, that belongs to a set, based on a proper
% function. The supported reparametrizations are:
% 
% if a > f > 0  then f = a/(1+exp(-x)             x  belongs to R       (1)
% if f > a      then f = a+exp(x)                 x  belongs to R       (2)
% if a > f > -a then f = -a+2aexp(x)/(1+exp(x))   x  belongs to R       (3)
% if a+b > f >a then f = a + bexp(x)/(1+exp(x))   x  belongs to R       (4)

% df and d2f are the first and second derivative of f, wrt x

% USAGE: For reparametrization to estimate a constrained model with
%        fminunc. The derivatives are needed for the chain rule.
% INPUTS: 
% x:        A vector of constrained values
% a:        Parameters (see above)
% flag:     integer, with values 1,2 or 3 (see above)

switch flag
    
    case 1 % x belongs to (0,a)
        f = a./(1+exp(-x));
        
    case 2 % x belongs to (a, +inf)
        f = a + exp(x);
        
    case 3 % x belongs to (-a,a)
        f = -a + 2*a*exp(x)./(1+exp(x));
        
    case 4 % x belongs to (a,a+b) for b = 200;
        if isvector(a) == 1
            b = a(2); a = a(1);
        else
            b = 98;
        end
        f = a + b*exp(x)./(1+exp(x));
end