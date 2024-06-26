function [lower,upper,weight] = getWeights(x,y,meth)

% The function returns three vectors:
% - lower: position in y of the nearest element below each x.
% - upper: position in y of the nearest element above each x.
% - weight: weight of lower in the linear combination that gives x.

% x: A vector representing values to be interpolated.
% y: A vector representing the grid or sample points.
% metodo: A string indicating the method, used for capping values beyond the boundaries.

% Find nearest elements in y above and below of x
    % Find the lowest element
    sizeX= length(x);
    sizeY= length(y);

    % Determine the indices of the bins in y to which each element in x belongs. 
    % The resulting indices are stored in the lower vector.
    [~,~,lower] = histcounts(x,y);

    %  Adjust the lower indices for elements in x that are beyond the boundaries of y
    if strcmp(meth, 'cap')
        lower(x>y(end)) = sizeY-1;
        lower(x<y(1)) = 1;
    end

    % Upper bound
    upper = lower+1;

    % Computing the weight of the element below
        weight = (y(upper) - x) ./ (y(upper) - y(lower));
        % upper element weight: 1 - weight







