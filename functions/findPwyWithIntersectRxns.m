function pwyIds=findPwyWithIntersectRxns(mnet,rxns)
%Find the intersection of pathways having specific reactions
%
% pwyIds=findPwyWithIntersectRxns(mnet,rxns)
%
%INPUTS
% mnet structure results from EFMTools
% rxns  list of reactions
%OUTPUTS
% pwyIds  list of pathway ID 

if ~iscell(rxns)
    rxns={rxns};
end
pwyIds=1:size(mnet.efms,2);
for i=1:length(rxns)
    rxn=rxns{i};
    rxnId= strcmp(mnet.reactionNames,rxn)==1;
    tmpPwyIds=find(mnet.efms(rxnId,:)~=0);
    pwyIds=intersect(pwyIds,tmpPwyIds);
end
pwyIds=unique(pwyIds);