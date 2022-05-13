function RepairedString = parseArgValue(PassedValue, AllowedLowercaseStrings)
% Remove hyphens from the Passed Value, then lowercase it, then compare it
% to AllowedLowercaseStrings using parseArgs. Return the RepairedString if
% found, otherwise error (parseArgs throws the error).
% AllowedLowercaseStrings should be the full set of hyphen-free, lowercase,
% full-length, allowed strings.

%   Copyright 2016 The MathWorks, Inc.

Defs = cell(1,numel(AllowedLowercaseStrings));
Values = cell(1,numel(AllowedLowercaseStrings));
[Values{:},set] = internal.stats.parseArgs(AllowedLowercaseStrings, Defs, ...
    lower(removeHyphens(PassedValue)), 'dummyval');
SetIdx = find(cell2mat(struct2cell(set)),1,'first');
if ~isempty(SetIdx)
    RepairedString = AllowedLowercaseStrings{SetIdx};
else
    RepairedString = [];
end
end

function S = removeHyphens(S)
S(S=='-') = '';
end
