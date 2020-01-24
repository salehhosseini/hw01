function dfdx = fjacobian(f, x0, dx)
% FGRADIENT   Computes Jacobian (partial derivative or gradient) of function
% dfdx = fgradient(@f, x0 [, dx])
%   Uses finite differences to compute partial derivative of function f
%   evaluated at vector x0. Optional argument dx specifies finite difference
%   step. Note that argument f should typically be entered as a
%   function handle, e.g. @f.
%   When dx is not given or empty, a default of 1e-6 is used.
  if nargin < 3 || isempty(dx)
    dx = 1e-6;   % with no input provided
  end
  f0 = f(x0);
  nX = length (x0);
  nY = length (f0);
  df = zeros(nY,nX); % store row vector for gradient (m x n matrix)
  
   for i = 1: nX
    xperturbed = x0;
    xperturbed(i) = xperturbed(i) + dx;
    df(:,i) = f(xperturbed) - f0;
  end
  dfdx = df / dx;
  
end