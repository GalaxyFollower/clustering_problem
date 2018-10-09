% random density checks
nTests = 1000;

%tabulated storage of results
centres = zeros( nTests, 3 );
radii = zeros( nTests, 1 );
densities = zeros( nTests, 1 );
points = zeros( nTests, 1 );
areas = zeros( nTests, 1 );

for i = 1:nTests
	    centre = randPointOnSphere( R );
	    radius = getRandDistance( R );
	    scannedArea = pi * radius^2;
	    
        % scan area for number of data points
	    nPoints = 0;
        for j = 1:length(u)
            somePoint = cpdata( j, 1:3 );
            if distance( somePoint, centre ) < radius
                nPoints = nPoints + 1;
            end
        end  
        
	    pointDensity = nPoints/scannedArea;      
	    %tabulate
        centres(i,:) = centre;
        radii(i) = radius;
        densities(i) = ceil(pointDensity);
        points(i) = nPoints;
        areas(i) = scannedArea;
        i
end

[~, I] = sort(densities);
densitySort = [centres(I,:) radii(I) areas(I) points(I) densities(I)];
[~, I] = sort(radii);
radiiSort = [centres(I,:) radii(I) areas(I) points(I) densities(I)];

% emphasis on u at centre of check area
	
%view(-10, 60)

function answer = getRandDistance( sphereRadius )
	answer = 0 + 2 * sphereRadius * rand;
end

function pointOnSphere = randPointOnSphere( radius )
    randPoint = -1 + 2 .* rand(1,3);
    pointOnSphere = projectToSphere( randPoint, radius );
end

function pointOnSphere = projectToSphere(point, radius)
    if(length(point) ~= 3)
        disp('error projecting to sphere');
    else
        scaleFactor = radius/norm(point);
        pointOnSphere = point * scaleFactor;
    end
end

