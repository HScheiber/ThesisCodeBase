function w = Wigner3j( j123, m123 )
% Compute the Wigner 3j symbol using the Racah formula. 
%
% W = Wigner3j( J123, M123 ) 
%
% J123 = [J1, J2, J3].
% M123 = [M1, M2, M3].
% All Ji's and Mi's have to be integeres or half integers (correspondingly).
%
% According to seletion rules, W = 0 unless:
%   |Ji - Jj| <= Jk <= (Ji + Jj)    (i,j,k are permutations of 1,2,3)
%   |Mi| <= Ji    (i = 1,2,3)
%    M1 + M2 + M3 = 0
% 
% Reference: 
% Wigner 3j-Symbol entry of Eric Weinstein's Mathworld:
% http://mathworld.wolfram.com/Wigner3j-Symbol.html
%
% Inspired by Wigner3j.m by David Terr, Raytheon, 6-17-04
%  (available at www.mathworks.com/matlabcentral/fileexchange).
%
% By Kobi Kraus, Technion, 25-6-08.
% Updated 1-8-13.
j1 = j123(1); j2 = j123(2); j3 = j123(3);
m1 = m123(1); m2 = m123(2); m3 = m123(3);
% Input error checking
if any( j123 < 0 ),
    error( 'The j must be non-negative' )
elseif any( rem( [j123, m123], 0.5 ) ),
    error( 'All arguments must be integers or half-integers' )
elseif any( rem( (j123 - m123), 1 ) )
    error( 'j123 and m123 do not match' );
end
% Selection rules
if ( j3 > (j1 + j2) ) || ( j3 < abs(j1 - j2) ) ... % j3 out of interval
   || ( m1 + m2 + m3 ~= 0 ) ... % non-conserving angular momentum
   || any( abs( m123 ) > j123 ), % m is larger than j
    w = 0;
    return
end
    
% Simple common case
if ~any( m123 ) && rem( sum( j123 ), 2 ), % m1 = m2 = m3 = 0 & j1 + j2 + j3 is odd
    w = 0;
    return
end
% Evaluation
t1 = j2 - m1 - j3;
t2 = j1 + m2 - j3;
t3 = j1 + j2 - j3;
t4 = j1 - m1;
t5 = j2 + m2;
tmin = max( 0,  max( t1, t2 ) );
tmax = min( t3, min( t4, t5 ) );
t = tmin : tmax;
w = sum( (-1).^t .* exp( -ones(1,6) * gammaln( [t; t-t1; t-t2; t3-t; t4-t; t5-t] +1 ) + ...
                         gammaln( [j1+j2+j3+1, j1+j2-j3, j1-j2+j3, -j1+j2+j3, j1+m1, j1-m1, j2+m2, j2-m2, j3+m3, j3-m3] +1 ) ...
                         * [-1; ones(9,1)] * 0.5 ) ) * (-1)^( j1-j2-m3 );
         
% Warnings
if isnan( w )
    warning( 'MATLAB:Wigner3j:NaN', 'Wigner3J is NaN!' )
elseif isinf( w )
    warning( 'MATLAB:Wigner3j:Inf', 'Wigner3J is Inf!' )
end