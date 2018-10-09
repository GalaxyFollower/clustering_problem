%fixme: colorbar sometimes randomly gets associatedd with figure 2
% idont know how to tie colorbar to figure 1. this method might be more
% accurate and 65 would not look like a magic number as its the # of discrete
% colour levels
%c = colorbar;
%uSteps = linspace(c.Limits(1), c.Limits(2), 65)';

uSteps = linspace(min(u), max(u), 65)';

%recommendation:
% choose centre band and bandwidth such that broken rings can be avoided 
centreBand = 45; %ideal 45-50
nSideBands = 2; %ideal 2-3
uMax = uSteps( centreBand + nSideBands + 1 );
uMin = uSteps( centreBand - nSideBands );

i1 = find( u <= uMax );
i2 = find( u >= uMin );
tempIndex = intersect( i1, i2 );

filteredData = [cpx(tempIndex) cpy(tempIndex) cpz(tempIndex)];% u value is irrelevant now
showData( filteredData, 4 )
