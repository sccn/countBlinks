% eegplugin_countBlinks() - A solution for manually counting blink-related IC activations in a single channel/IC. 

% Copyright (C) 2018 Makoto Miyakoshi, SCCN,INC,UCSD.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1.07  USA

function vers = eegplugin_countBlinks(fig, trystrs, catchstrs)
    
vers = 'beta';
if nargin < 3
    error('eegplugin_arfitStudio requires 3 arguments');
end;
    
% create menu
toolsmenu = findobj(fig, 'tag', 'tools');
uimenu( toolsmenu, 'label', 'Count blinks', 'separator','on','callback', 'countBlinks');