function [out] = shuffleVector (in)

% Vector shuffle
out = in(randperm(length(in)), :);

end