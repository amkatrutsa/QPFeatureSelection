function [idx] = FilterMasks(masks, maxItems)
% Function filters masks and leaves only one of which type with items 
% no more than maxItems.
%
% Input:
% masks    - matrix with masks as columns
% maxItems - masks with more items than max will be droped
% 
% Output:
% idx - indices of rest columns

idx = 1:size(masks, 2); % all masks at first
iMask1 = 1; % first iterator
while iMask1 <= length(idx) % for all masks in result
    if (length(find(masks(:, idx(iMask1)))) < 1) || ... 
       (length(find(masks(:, idx(iMask1))))) > maxItems % too many items
        idx(iMask1) = []; % drop it
    else
        iMask2 = iMask1 + 1; % second iterator
        while iMask2 <= length(idx) % for all rest masks
            if masks(:, idx(iMask1)) == masks(:, idx(iMask2)) % masks are equivalent
                idx(iMask2) = []; % drop second mask
            else
                iMask2 = iMask2 + 1; % next second mask
            end
        end
        iMask1 = iMask1 + 1; % next first mask
    end
end
end