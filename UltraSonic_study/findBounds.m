function [bounds] = findBounds(anchorLoc)
for i = 1:size(anchorLoc, 2)
  bounds(1, i) = min(anchorLoc(:, i)); % minimum boundary of ith axis
  bounds(2, i) = max(anchorLoc(:, i)); % maximum boundary of ith axis
end
% hard coded minimum height (0 m) of search boundary
bounds(1, end) = 0;