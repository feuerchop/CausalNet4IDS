function [f, df, d2f] = FromCons2Unc(x,a,flag)

% Let x be a constrained parameter. This function takes x and
% transforms it to a parameter f, that belongs to R, based on the interval
% that x belongs. For example
% 
% if a > x > 0 then f = log(x/(a-x)) belongs to R        (1)
% if x > a then f = log(x - a) belongs to R              (2)
% if a > x > -a then f = log((a+x)/(a-x)) belongs to R   (3)
% if a+b > x > a then f = log((a-x)/(x-a-b)) belongs 2 R (4)

% df and d2f are the first and second derivative of f, wrt x

% USAGE: For reparametrization to estimate a constrained model with
%        fminunc. The derivatives are needed for the chain rule.
% INPUTS: 
% x:        A vector of constrained values
% a:        Parameters (see above)
% flag:     integer, with values 1,2 or 3 (see above)

switch flag
    
    case 1 % x belongs to (0,a)
        % first check if x belongs in (0,a)
        if max(max(x))>a || min(min(x))<0
            error('the constrained parameter should belong to (0,a) but it is not')
        end
        f = log(x./(a-x));
        df = 1./x + 1./(a - x);
        d2f = 1./(a - x).^2 - 1./x.^2;
        
    case 2 % x belongs to (a, +inf)
        if min(min(x))<a
            error('x should be larger than a, but it is not')
        end
        f = log(x - a);
        df = 1./(x - a);
        d2f = -1./((x - a).^2);
        
    case 3 % x belongs to (-a,a)
        if max(max(x))> a || min(min(x))<-a
            error('x should lie in (-a,a) but it is not')
        end
        f = log((a+x)./(a-x));
        df = 2*a./(a^2 - x.^2);
        d2f = 4*a*x./(a^2 - x.^2).^2;
        
    case 4
        if isvector(a)==1
            b = a(2); a = a(1);
        else
            b = 98;
        end
        f = log((a-x)./(x-b-a));
        df = -1./(a + b - x) - 1./(a - x);
        d2f = -1./(a + b - x).^2 - 1./(a - x).^2;
end