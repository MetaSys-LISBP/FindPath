function model=makeCBModelFromCSV(fileName,delimiter,modelName,metExt,lbExt,ubExt)
%model=makeCBModelFromCSV(filename,delimiter)
%Create a model from a CSV file
%INPUTS
%   fileName   name of the CSV file
%   delimiter   delimiter of the CSV file
%   modelName Name of the ouput CBModel file
%   metExt     list of external metabolites
%   lbExt        lower bound for the exchange flux of external metabolites
%   ubExt       upper bound for the exchange flux of external metabolites
%OUPUT
%   model    A CBModel
data=csvimport(fileName,'delimiter',delimiter);
data=data(2:end,:);
model=createCBModel(data(:,1),data(:,2),data(:,3),cell2mat(data(:,4)));
if nargin<4
 model=   addExchangeRxn(model,findExtMet(model))
else
model=addExchangeRxn(model,metExt,lbExt,ubExt);
end
writeCbModel(model,'sbml',modelName);