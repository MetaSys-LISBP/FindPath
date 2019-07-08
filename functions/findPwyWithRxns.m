function pwyIds=findPwyWithRxns(mnet,rxns)
%Find the intersection of pathways having specific reactions
%
% pwyIds=findPwyWithRxns(mnet,rxns)
%
%INPUTS
% mnet structure results from EFMTools
% rxns  list of reactions
%OUTPUTS
% pwyIds  list of pathway ID 

if ~iscell(rxns)
    rxns={rxns};
end
pwyIds=[];
for i=1:length(rxns)
    rxn=rxns{i};
    rxnId= strcmp(mnet.reactionNames,rxn)==1;
    tmpPwyIds=find(mnet.efms(rxnId,:)~=0);
    pwyIds=union(pwyIds,tmpPwyIds);
end
pwyIds=unique(pwyIds);