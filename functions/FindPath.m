function [pwys,rankFinal]=FindPath(model,dbModel,EFMToolPath,options)
%Design pathways from a set of reactions
%
% [pwys,score]=FindPath(model,dbModel,options)
%
%INPUTS
% model         CBM host for metabolic modifications (sbml file or model object)
% dbModel     Set of reactions to design the new pathway (sbml file or model object)
% EFMToolPath      Path to the EFMTool package
% options          structure with fields
%   rxnIn               reactions starting the metabolic conversion process
%   rxnMust          reactions included in the metabolic conversion process
%   rxnOut            reactions ending the metabolic conversion process
%   rxnForb           reactions to avoid in the metabolic conversion process
%   metIn              metabolites starting the metabolic conversion process
%   metMust         metabolites included in the metabolic conversion process
%   metOut           metabolites ending the metabolic conversion process
%   metForb          metabolites to avoid in the metabolic conversion process
%   maxLength      maximal number of reactions in the pathways
%   namePwy        name for pathways model
%   maxPwy         maximal number of solution pathways
%   ignoreRxns    List of reactions to not take into account
%   weightEfficiency   Weight of the pathway efficeincy in the final ranking
%   weightLength  Weight of the pathway length in the final ranking
%
%OUTPUTS
% pwys list of pathways
% score score of the pathways


%%
%options
%t=clock;

if ischar(model)
    model=readCbModel(model);
    
end
if ischar(dbModel)
    
    dbModel=readCbModel(dbModel);
    
end
if nargin<4
    disp('No options')
    return;
    
end

if ~isfield(options,'weightLength')
    options.weightLength=1;
end

if ~isfield(options,'weightEfficiency')
    options.weightEfficiency=1;
end

%exchange flux
if ~isfield(options,'exchangeFlux')
    
    if  ~isfield(options,'rxnIn') && ~isfield(options,'metIn')
        
        disp('no pathway input')
        return
    else
        metTmp={};
        if isfield(options,'rxnIn')
            metTmp= findSubOfRxn(dbModel,options.rxnIn);
            
        end
        
        if isfield(options,'metIn')
            metTmp=unique(union(metTmp,options.metIn));
        end
        exRxns=dbModel.rxns(findExcRxns(dbModel));
        
        tmpExFlux={};
        
        for i=1:length(metTmp)
            
            tmpMet=regexprep(metTmp(i),'\[e\]$','(e)');
            if isempty(intersect(exRxns,findRxnsFromMets(dbModel,tmpMet)))
                dbModel=addExchangeRxn(dbModel,tmpMet);
                tmpExFlux{end+1}=dbModel.rxns{end};
            else
                
                tmpExFlux{end+1}=intersect(exRxns,findRxnsFromMets(dbModel,tmpMet));
                
            end
            
        end
        
    end
    
    options.exchangeFlux=tmpExFlux;
end
if ~isfield(options,'maxPwy')
    options.maxPwy=1000;
end



[pwys,totalPwy,topoPwy,lengthPwy,rankLength]=computePwy(dbModel,EFMToolPath,options);
[scoreOptimized,rankEfficiency]=evaluatePwy(model,dbModel,pwys,options);
if isempty(rankLength)  || isempty(scoreOptimized)
    disp('No pathway found');
    pwys={};
    score=[];
    return
end

options.weightEfficiency

rankFinal=options.weightEfficiency.*rankEfficiency'+options.weightLength.*rankLength./min(rankLength);
[rankFinal,index]=sort(rankFinal);
pwys=cellfun(@(x) pwys{index(x)},num2cell(1:length(pwys)), 'UniformOutput',false);

exportPwys(dbModel,pwys,options)
disp(['Total pathways: ' num2str(totalPwy)]);
disp(['Metabolic composition filter: ' num2str(topoPwy)]);
disp(['Length filter: ' num2str(lengthPwy)]);
disp(['Export filter: ' num2str(options.maxPwy)]);

fileID = fopen('findPath.log','w');
fprintf(fileID,'Total pathways: %4.0f\n',totalPwy);
fprintf(fileID,'Metabolic composition filter: %4.0f\n',topoPwy);
fprintf(fileID,'Length filter: %4.0f\n',lengthPwy);
fprintf(fileID,'Export filter: %4.0f\n',options.maxPwy);
fprintf(fileID,'\npathways information:\n' );
for i=1:min(length(index),options.maxPwy)

    fprintf(fileID,'pathways name: %s, length size: %3.0f, efficiency value: %2.6f\n',strcat(options.namePwy,'_',num2str(i)),rankLength(index(i)),scoreOptimized(index(i)));
        for j=1:length(pwys{i})
        tmp=printRxnFormula(dbModel,pwys{i}{j},false);
        fprintf(fileID,'%s: %s \n',strcat(pwys{i}{j}),strcat(tmp{:}));
        end
     fprintf(fileID,'\n');
end

fclose(fileID);
end


function exportPwys(dbModel,pwys,options)
for i=1:min(length(pwys),options.maxPwy)
    pwy=pwys{i};
    pwy={pwy{:} options.exchangeFlux{:}};
    pwyModel=extractSubNetwork(dbModel,pwy);
    strcat(options.resPath,'/',options.namePwy,'_',num2str( i))
    
    writeCbModel(pwyModel,'sbml',strcat(options.resPath,'/',options.namePwy,'_',num2str( i)));
    
end
end

function [pwys,totalPwy,topoPwy,lengthPwy,rankLength]=computePwy(dbModel,EFMToolPath,options)
localFolder=pwd;
cd(EFMToolPath);
if ~isfield(options,'rxnIn')
    rxnIn={};
else
    if ~iscell(options.rxnIn)
        rxnIn={options.rxnIn};
    else
        rxnIn=options.rxnIn;
    end
end
if ~isfield(options,'rxnOut')
    rxnOut={};
else
    if ~iscell(options.rxnOut)
        rxnOut={options.rxnOut};
    else
        rxnOut=options.rxnOut;
    end
end
if ~isfield(options,'rxnMust')
    rxnMust={};
else
    if ~iscell(options.rxnMust)
        rxnMust={options.rxnMust};
    else
        rxnMust=options.rxnMust;
    end
end
if ~isfield(options,'rxnForb')
    rxnForb={};
else
    if ~iscell(options.rxnForb)
        rxnForb={options.rxnForb};
    else
        rxnForb=options.rxnForb;
    end
end


if ~isfield(options,'metIn')
    metIn={};
else
    if ~iscell(options.metIn)
        metIn={options.metIn};
    else
        metIn=options.metIn;
    end
end
if ~isfield(options,'metOut')
    metOut={};
else
    if ~iscell(options.metOut)
        metOut={options.metOut};
    else
        metOut=options.metOut;
    end
end

if ~isfield(options,'metMust')
    metMust={};
else
    if ~iscell(options.metMust)
        metMust={options.metMust};
    else
        metMust=options.metMust;
    end
end
if ~isfield(options,'metForb')
    metForb={};
else
    if ~iscell(options.metForb)
        metForb={options.metForb};
    else
        metForb=options.metForb;
    end
end

if ~isfield(options,'maxLength')
    maxLength=0;
else
    
    maxLength=options.maxLength;
    
end
if ~isfield(options,'ignoreRxns')
    ignoreRxns={};
else
    ignoreRxns=options.ignoreRxns;
    
end




stru=struct();
stru.stoich=full(dbModel.S);
stru.reversibilities=dbModel.rev;
stru.metaboliteNames=dbModel.mets;
stru.reactionNames =dbModel.rxns;
clear('mnet')
mnet = CalculateFluxModes(stru);
clear('mnet')
tmp=open('./tmp/efms_0.mat');
mnet=tmp.mnet;
mnet.efms(find(findExcRxns(dbModel)==1),:)=0;

cd(localFolder)
pause(1);

defautID=1:size(mnet.efms,2);
totalPwy=length(defautID);

if ~isempty(metIn)
    pwyMetIdIn=findPwyWithRxns(mnet,findRxnsFromMets(dbModel,metIn));
else
    pwyMetIdIn=[];
end

if ~isempty(metOut)
    pwyMetIdOut=findPwyWithRxns(mnet,findRxnsFromMets(dbModel,metOut));
else
    pwyMetIdOut=[];
end
if~isempty(metForb)
    pwyMetIdForb=findPwyWithIntersectRxns(mnet,findRxnsFromMets(dbModel,metForb));
else
    pwyMetIdForb=[];
    
end
if ~isempty(metMust)
    pwyMetIdMust=findPwyWithIntersectRxns(mnet,findRxnsFromMets(dbModel,metMust));
else
    pwyMetIdMust=defautID;
end




if ~isempty(rxnIn)
    pwyRxnIdIn=findPwyWithRxns(mnet,rxnIn);
else
    pwyRxnIdIn=[];
end

if ~isempty(rxnOut)
    pwyRxnIdOut=findPwyWithRxns(mnet,rxnOut);
else
    pwyRxnIdOut=[];
end
if~isempty(rxnForb)
    pwyRxnIdForb=findPwyWithIntersectRxns(mnet,rxnForb);
else
    pwyRxnIdForb=[];
    
end
if ~isempty(rxnMust)
    pwyRxnIdMust=findPwyWithIntersectRxns(mnet,rxnMust);
else
    pwyRxnIdMust=defautID;
end



pwyInId=union(pwyMetIdIn,pwyRxnIdIn);

pwyOutId=union(pwyMetIdOut,pwyRxnIdOut);

pwyMustId=intersect(pwyMetIdMust,pwyRxnIdMust);

pwyForbId=union(pwyMetIdForb,pwyRxnIdForb);


pwyIds= setdiff(intersect(intersect(pwyInId,pwyOutId),pwyMustId),pwyForbId);
topoPwy=length(pwyIds);
if maxLength>0
    pwyLengthId=findPwyWithLength(mnet,maxLength);
    
else
    pwyLengthId=defautID;
end

pwyIds=intersect(pwyLengthId, pwyIds);

pwys=findPwyFromPwyIds(mnet,pwyIds);

pwys=cellfun(@(x) setdiff(pwys{x},ignoreRxns),num2cell(1:length(pwys)), 'UniformOutput',false);
pwys=cellfun(@(x) sort(pwys{x}),num2cell(1:length(pwys)), 'UniformOutput',false);
pwy2str=cellfun(@(x) strcat(pwys{x}{:}),num2cell(1:length(pwys)), 'UniformOutput',false);
[~,IA]=unique(pwy2str);
pwys=pwys(sort(IA));
lengthPwy=length(pwyIds);
[pwys,rankLength]=findPwysLength(pwys);

end

function [scoreOptimized,rankEfficiency]=evaluatePwy(model,dbModel,pwys,options)

scoreOptimized=zeros(length(pwys),1);
if ~isfield(options,'exchangeFlux')
    return;
else
    for i=1:length(pwys)
        pwy=pwys{i};
        pwy=union(pwy,options.exchangeFlux);
        modelTmp=addPwy(model,pwy,dbModel);
        modelTmp.lb(findRxnIDs(modelTmp,options.exchangeFlux))=-10;
        opti=optimizeCbModel(modelTmp);
        scoreOptimized(i)=opti.f;
    end
end

rankEfficiency=max(scoreOptimized)-scoreOptimized;
end

