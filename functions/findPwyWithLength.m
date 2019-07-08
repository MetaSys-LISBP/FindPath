function pwyIds=findPwyWithLength(mnet,maxLength)
%Find pathways having a length under or equal to the given value
%
%pwyIds=findPwyWithLength(mnet,maxLength)
%
%INPUTS
% mnet structure results from EFMTools
% maxLength  maximal size of the pathways
%OUTPUTS
% pwyIds  list of pathways ID 
index=1:size(mnet.efms,2);
vectLength=cellfun(@(x) length(find(mnet.efms(:,x)~=0)),num2cell(1:size(mnet.efms,2)), 'UniformOutput',false);
vectLength=cell2mat(vectLength);
pwyIds=index(vectLength<=maxLength);
%[vals,index]=sort(vectLength);
%pwyIds=index(vals<=maxLength);

