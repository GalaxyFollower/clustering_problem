%use data = shuffle(data); for overwrite
function data = shuffle( data )
    x = data(:,1);
    y = data(:,2);
    z = data(:,3);
    group = data(:,4);

    s = size(data);
    order = randperm(s(1));

    x = x(order);
    y = y(order);
    z = z(order);
    group = group(order);

    data = [x y z group];
end