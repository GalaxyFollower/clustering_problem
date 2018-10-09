[temp, I] = sort(u);
I = flipud(I);

% gradual scatter with highest points first. u visual aid
speed = 200; % can be set to 0.01 for ultra slow scatter
cap = length(u);

for i = 1:speed:cap
    consideredDataPoints = I(1:ceil(i));
    
    CPX = cpx(consideredDataPoints);
    CPY = cpy(consideredDataPoints);
    CPZ = cpz(consideredDataPoints);
    CPU = u(consideredDataPoints);

    figure(3); clf;
    scatter3(CPX, CPY, CPZ, '.');%, [], U, 200
    axis equal;
end

