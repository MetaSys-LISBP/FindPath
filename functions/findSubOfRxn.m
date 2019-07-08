function mets = findSubOfRxn(model, rxns)
%  findSubOfRxn    Find metabolites substrate involved in reaction(s)
%
%mets = findSubOfRxn(model, rxns)
%
%INPUTS
% model         CBModel
% rxns            List of reactions
%
%
%OUTPUTS
% mets  List of metabolites

mets = model.mets(model.S(:,findRxnIDs(model, rxns))< 0);


