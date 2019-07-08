function pwyModel=addPwy(model,pwy,dbModel)
% Add the pathway in a model
%
% pwyModel=addpwys(model,dbModel,pwy)
%
%INPUTS
% model         CBM host for metabolic modifications (sbml file or model object)
% pwy             List of reactions or sbml file 
% dbModel      CBM model if pwy is a list of reactions
%
%
%OUTPUTS
% pwyModel  model structure with the new pathway
if ischar(pwy)
    tmpModel=readCbModel(pwy);
    pwyModel=mergeTwoModels(model,tmpModel);
    
else
    
    rxns=setdiff(pwy,model.rxns);
    pwyModel=model;
    for i=1:length(rxns)
        rxn=rxns{i};
        mets=findMetsOfRxn(dbModel,rxn);
        metsId=findMetIDs(dbModel,mets);
        rxnId=findRxnIDs(dbModel,rxn);
        coeff=dbModel.S(metsId,rxnId);
        revF=dbModel.rev(rxnId);
        [pwyModel,~] = addReaction(pwyModel,rxn,mets,coeff,revF);
    end
    
end
end