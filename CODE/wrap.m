function wrapped = wrap(imatrix)
% WRAP(MATRIX)  wrap (phase values) to principal interval.
%   MAT = WRAP(MATRIX); wrap to principal interval
%   return MAT wrapped MATRIX to principal interval [-pi,pi)
%   wrapped = mod(imatrix+pi,2*pi) - pi; (note: not rem).
%
%   See also RESIDUES, MOD, REM.
%

% $Revision: 1.6 $  $Date: 2001/09/28 14:24:34 $
% Bert Kampes, 01-Mar-2000

if (~isreal(imatrix)) helphelp; return; end;


% old: slow?
% old: complex: 
%wrapped = atan2(sin(imatrix),cos(imatrix));
%wrapped = mod(imatrix,2*pi);
%wrapped = rem(imatrix+pi,2*pi) - pi;
% not rem! (signed)
wrapped = mod(imatrix+pi,2*pi) - pi;

% might be more efficient for only few 2b wrapped, but not in general.
% xxx          = find(abs(imatrix)<pi);
% imatrix(xxx) = mod(imatrix(xxx)+pi,2*pi) - pi;
% wrapped      = imatrix(xxx);



%%% EOF
