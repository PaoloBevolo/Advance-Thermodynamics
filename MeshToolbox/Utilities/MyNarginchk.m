function MyNarginchk(args, minargs, maxargs)
%Validate number of input arguments, similarly to narginchk, for use in older MATLAB versions
%From Mathworks we site: 
%narginchk(minargs, maxargs) throws an error if the number of inputs 
%specified in the call to the currently executing function is less than 
%minargs or greater than maxargs. If the number of inputs is between 
%minargs and maxargs (inclusive), narginchk does nothing.
%
%Based on the narginchk octave function, Copyright (C) 2011 Carnë Draug

%Version 2017.1
%Copyright 2014-2017 Paolo Bardella

  if (nargin ~= 3)
    error('MyNarginchk:help','%s: Usage: MyNarginchk(args, minargs, maxargs)',upper(mfilename));
  elseif (~isnumeric (minargs) || ~isscalar (minargs))
    error ('MyNarginchk:WrongInputArguments','minargs must be a numeric scalar');
  elseif (~isnumeric (maxargs) || ~isscalar (maxargs))
    error ('MyNarginchk:WrongInputArguments','maxargs must be a numeric scalar');
  elseif (minargs > maxargs)
    error ('MyNarginchk:WrongInputArguments','minargs cannot be larger than maxargs')
  end
  %args = evalin ('caller', 'nargin;');
  if (args < minargs)
    error ('MyNarginchk:notEnoughInputs','Not enough input arguments.');
  elseif (args > maxargs)
    error ('MyNarginchk:tooManyInputs','Too many input arguments.');
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Original version from 
%https://gitorious.org/octave-carandraug/inputparser/source/839bbd4ee5e5d0c7d02c5c89580e45b73e549aed:narginchk.m#L1-2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ## Copyright (C) 2011 Carne Draug <carandraug+dev@gmail.com>
% ##
% ## This program is free software; you can redistribute it and/or modify
% ## it under the terms of the GNU General Public License as published by
% ## the Free Software Foundation; either version 3 of the License, or
% ## (at your option) any later version.
% ##
% ## This program is distributed in the hope that it will be useful,
% ## but WITHOUT ANY WARRANTY; without even the implied warranty of
% ## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% ## GNU General Public License for more details.
% ##
% ## You should have received a copy of the GNU General Public License
% ## along with this program; if not, see <http://www.gnu.org/licenses/>.
% 
% ## -*- texinfo -*-
% ## @deftypefn {Function File} {} narginchk (@var{minargs}, @var{maxargs})
% ## Checks for correct number of arguments.
% ##
% ## This function returns an error unless the number of arguments in its caller
% ## is between the values of @var{minargs} and @var{maxargs}. It does nothing
% ## otherwise.
% ##
% ## Both @var{minargs} and @var{maxargs} need to be a numeric scalar. Zero, Inf
% ## and negative are all valid, and they can hev the same value.
% ##
% ## Note that this function evaluates the value of @code{nargin} on the caller so
% ## its value must have not been tampered with.
% ##
% ## @seealso{nargin, nargchk, nargoutchk, inputParser}
% ## @end deftypefn
% 
% function narginchk (minargs, maxargs)
% 
%   ## it requires always two arguments (can't specify only min)
%   ## zero, negative and inf are all valid arguments and they can be equal
%   ## thanks to Oldak in ##matlab for the help in checking these corner cases
%   ## tested compatibility in version 2011b
% 
%   if ( nargin != 2 )
%     print_usage;
%   elseif ( !isnumeric (minargs) || !isscalar (minargs) )
%     error ("minargs must be a numeric scalar");
%   elseif ( !isnumeric (maxargs) || !isscalar (maxargs) )
%     error ("maxargs must be a numeric scalar");
%   elseif ( minargs > maxargs )
%     error ("minargs cannot be larger than maxargs")
%   endif
% 
%   args = evalin ('caller', ...
%                  'nargin;' ...
%                  );
% 
%   if ( args < minargs )
%     error ("Not enough input arguments.");
%   elseif ( args > maxargs )
%     error ("Too many input arguments.");
%   endif
% 
% endfunction
% 
