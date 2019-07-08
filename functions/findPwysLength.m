function [pwys,lengthPwys]=findPwysLength(pwys)
% [pwys,lengthPwys]=findPwysLength(pwys)
%  Find the length of the pathways
%
%INPUTS
% pwys list of pathways
%OUTPUTS
% pwys pathways sorted by length
% lengthPwys  list of pathways length

lengthTmp=cellfun(@(x) length(pwys{x}),num2cell(1:length(pwys)), 'UniformOutput',false);
lengthPwys=cell2mat(lengthTmp);
