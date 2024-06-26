function [linInterp,auxInterp] = linInterp(x,y,z,auxInterp)
% This function provides a linear interpolation for a function at x, using
% known pairs of coordinates (vector y) and their images (vector z).
% x can be a vector. Optionally, to save computing time it is possible to
% directly provide the nearest elements and weight in auxInterp.

% Besides the interpolation value, the function returns three vectors:
% - lower: position in y of the nearest element below each x.
% - upper: position in y of the nearest element above each x.
% - weight: weight of lower in the linear combination that gives x.


% % Check size of vectors
%     if length(y) ~= length(z)
%         error(message("Each element in y must be paired with a single element in z"))
%     end


% Obtaining nearest values and weights
    if nargin == 3 % no help has been provided
        % Finding elements in y immediately above and below x
            % Number of elements in x:
                siz = length(x);
            % Number of elements in y:
                last = length(y);
            % Find upper and lower elements for each of them
                for i=1:siz
                    lower = find(y<x(i),1,"last");
                    if isempty(lower)
                        lower(i)=0;
                    else
                        lower(i)=lower;
                    end
                end
                upper = lower+1;
                % When x(i) is the first element in y, the lower neighbour
                % will not matter, all the weight will be on upper. But we
                % need lower to remain on the grid, so it cannot be 0. It
                % cannot be 1 either, or we will have lower=upper and
                % weight will not be well defined.
                lower(lower==0) = 2;
            % Computing the weight of the element below
                weight = (y(upper) - x) ./ (y(upper) - y(lower));
                % the weight for the upper element is (1 - weight)
    else % the information has been externally provided
        lower = auxInterp{1};
        upper = auxInterp{2};
        weight = auxInterp{3};
    end

% Interpolating
    linInterp = (weight .* z(lower)) + (1 - weight).*z(upper);

% Other returns
    auxInterp = {lower, upper, weight};