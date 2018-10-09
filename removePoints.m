% this calculation assumes uniform square grid of data points
avgDensity = length(u) / (4 * pi * R^2);
closeDistance = ( sqrt(2) / sqrt(avgDensity) );

% some definitions:
% knnRadius: any knn serarch returns points that lie within a certain distance of the root point


Multiplier = 20; 
closenessThreshold = Multiplier * closeDistance;
% closeness tolerance i.e multiplier must be such that "jumpy/broken rings"
% (which may occur because of narrow bandwidth ans non uniform sample point
% density on sphere surface) dont cause problems.
% To see broken rings, reduce nSideBands parameter in filteru.m

k = 17;
% 1. AS HIGHER K IS CHOSEN, HIGHER closenessThreshold MUST ALSO BE CHOSEN!
% 	BECAUSE THE FARTHEST POINT AMONG THE K nearest neighbours WILL BE FATHER AWAY
%	for safety, choose k such that closenessThreshold >= 2 * knnRadius
% 2. k must be such that "random walks" dont occur within closed bands of
% 	points and the bands are traversed in a linear fashion

rng('default');
% this is important: if our parameters from filteru.m or k or multiplier dont give desired results, 
% we would like to watch exactly what happened again so we can tweak parameters. so pseudo randomnumbers are ideal

nRings = 0;
% we count the rings in filteredData instead of number of spots on sphere.

while ~isempty(filteredData)

	group = [0 0 0]; %stores a group of points what belong to the same ring.
	groupingStatus = 1;

	%  choose a random point from ungrouped data aka raw data
	nPoints = size(filteredData, 1);
	rootPoint = filteredData(randi(nPoints, 1), :);

	while groupingStatus == 1 && ~isempty(filteredData)
	    
	    % scan the ring to which rootPoint belongs for k nearese neighbours
	    % from second iteration onward this works like a ray sweep of the thin ring.
	    [idx, distance] = knnsearch(filteredData, rootPoint, 'k', k);
	    if ~isempty(distance)
	       	closestNeiborsDist = distance(1);
		farthestNeigborsDist = distance(length(distance));
	    end
	    
	    % if farthest knn point lies "far" away, we have swept the whole ring.
	    % its time to stop putting them into the same group
	    if farthestNeigborsDist > closenessThreshold
		groupingStatus = 0;
	    end
	    
	    %for termination, some points knn search returns will lie in other rings.
	    closePointsIdx = find(distance <= closenessThreshold);
	    if ~isempty(closePointsIdx)
		% append points that are "close" to a group
		group = [group; filteredData(idx(closePointsIdx), :)];

		% shift root to the farthest point among knn data for the next ray sweep
		rootPoint = filteredData(idx(length(closePointsIdx)), :);

		% remove data we have just grouped from raw data
		filteredData(idx(closePointsIdx), :) = [];
	    end
	    
	    % visual aid
	    showData(group, 5)
	    pause(delay);
	end
	
	nRings = nRings + 1
end
