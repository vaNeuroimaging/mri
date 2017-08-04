done = 0;

base = [repmat([1],1,40) repmat([0],1,160)];

while done == 0
    
    order = randperm(200);
    test = base(order);
    
    locations = find(test);
    itis = locations(2:end) - locations(1:end-1);
    
    if mean(itis) > 4.99 && mean(itis) < 5.01 && std(itis) > 3.99 && std(itis) < 4.01
        done = 1;
    end
end
