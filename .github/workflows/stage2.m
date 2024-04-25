adapterLocation = fullfile('D:\xyModelGEM\xyModelGEMAdapter.m');
ModelAdapterManager.setDefault(adapterLocation);
ModelAdapter = ModelAdapterManager.getDefault();
params = ModelAdapter.getParameters();
ecModel = loadEcModel('ecModel.yml');

% 17-19
ecModel = getECfromGEM(ecModel);
noEC = cellfun(@isempty, ecModel.ec.eccodes);
ecModel = getECfromDatabase(ecModel, noEC);
kcatList_fuzzy = fuzzyKcatMatching(ecModel);

% 20-25
%[ecModel, noSMILES] = findMetSmiles(ecModel);
%writeDLKcatInput(ecModel);
%runDLKcat();
kcatList_DLKcat = readDLKcatOutput(ecModel);

%26-27
kcatList_merged = mergeDLKcatAndFuzzyKcats(kcatList_DLKcat, kcatList_fuzzy);
ecModel = selectKcatValue(ecModel, kcatList_merged);

%28 (optional?)
%[ecModel, rxnUpdated, notMatch] = applyCustomKcats(ecModel);

%29 (optional)
%ecModel = getKcatAcrossIsozymes(ecModel);

%30 (optional)
%[ecModel, rxnsMissingGPR, standardMW, standardKcat] = getStandardKcat(ecModel);

%31 
ecModel = applyKcatConstraints(ecModel);

%32
%f = calculateFfactor(ecModel);
ecModel = setProtPoolSize(ecModel, [], 0.5);

saveEcModel(ecModel,'stage2model.yml');