function mets=findExtMet(model,suffix)
%mets=findExtMet(model,suffix)
%find External metabolites based on a suffix
%INPUTS
%  model    A CBModel
%  suffix     The suffix of the external metabolites
%OUTPUT
%mets      List of external metabolites
if nargin<2
    suffix='[e]';
end
mets=model.mets(~cellfun('isempty',regexp(model.mets,suffix,'end')));
