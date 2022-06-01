function [loc, A, B] = findinSorted(x,range,err)
%%%%%%%%%%%%%%% inputs %%%%%%%%%%%%%%%
% x = value to find                  %
% range = (sorted) vector to search  %
% err = tolerence on search          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Outputs
A   = 1;            %left boundary
B   = numel(range); %right boundary
loc = 0;            %initially set to 0 implying not found
%% Binary Search Algorithm
while B-A > 1
    mid = floor((A+B)/2);
    if x > range(mid)
        B = mid;
    else
        A = mid;
    end
end
%% Returning returned value if found (within tolerance)
if abs(range(A) - x) < err
    loc = A;
    return
elseif abs(range(B) - x) < err
    loc = B;
    return
end