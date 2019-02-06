function Y = normalized(X)
maxX = max(X(:));
minX = min(X(:));
Y    = (X-minX)/(maxX-minX);
end