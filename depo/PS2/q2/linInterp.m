function [linInterp,auxInterp] = linInterp(x,y,z,auxInterp)
% This function uses linear interpolation for a function at x, using
% coordinates (vector y) and their images (vector z). To save time, we can use 
% auxInterp to get the nearest elements and weights. 

% - lower: position in y of the nearest element below each x.
% - upper: position in y of the nearest element above each x.

if nargin ==3
    % Find nearest elements in y above and below of x
        size =length(x);
        last = length(y);
        % Find upper and lower elements for x and y
                for i=1:size
                    try
                        upper(i) = find(y>=x(i),1,"first");
                    catch
                        upper(i) = size;
                    end
                end
                lower = upper-1;
                lower(lower==0) = 2;
            % Computing the weight of the element below
                weight = (y(upper) - x) ./ (y(upper) - y(lower));
    else % the information has been externally provided
        lower = auxInterp{1};
        upper = auxInterp{2};
        weight = auxInterp{3};
end

% Interpolating
    linInterp = (weight .* z(lower)) + (1 - weight).*z(upper);

% Other returns
    auxInterp = {lower, upper, weight};






