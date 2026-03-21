function bool = isbetween(value,boolrange)
    bool = value >= boolrange(1) & value <= boolrange(2);
end