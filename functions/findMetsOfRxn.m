function mets = findMetsOfRxn(model, rxns)
%Find metabolites involved in reaction(s)
%
% mets = findMetsOfRxn(model, rxns)
%
%INPUT
%  model         CBModel
%  rxns          reaction name (possibly a cell of several reaction
%   names)
%OUTPUT
% mets           list of metabolites



mets = model.mets(sum(abs(model.S(:,findRxnIDs(model, rxns))),2) ~= 0);

