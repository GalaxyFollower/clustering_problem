function showData( Data, figNo )
    figure(figNo);
    scatter3( Data(:,1), Data(:,2), Data(:,3), '.' );
    hold on;
    axis equal;
end
