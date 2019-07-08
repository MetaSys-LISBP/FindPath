function pwys=findPwyFromPwyIds(mnet,pwyIds)
%Extract pathways from the pathway IDs
%
% pwys=findPwyFromPwyIds(mnet,pwyIds)
%
%INPUTS
% mnet structure results from EFMTools
% pwyIds  list of the pathway IDs
%OUTPUTS
% pwys  list of pathways 
pwys=cellfun(@(x) mnet.reactionNames(mnet.efms(:,pwyIds(x))~=0),num2cell(1:length(pwyIds)), 'UniformOutput',false);