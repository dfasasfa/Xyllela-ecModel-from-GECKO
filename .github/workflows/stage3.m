adapterLocation = fullfile('D:\xyModelGEM\xyModelGEMAdapter.m');
ModelAdapterManager.setDefault(adapterLocation);
ModelAdapter = ModelAdapterManager.getDefault();
params = ModelAdapter.getParameters();
ecModel = loadEcModel('stage2model.yml');

%33-37
ecModel = setParam(ecModel, 'lb', 'EX_C00031__dra', -1000);

constructEquations(ecModel, 'EX_C00031__dra');

%biomass_reaction_index = find(contains(ecModel.rxns, 'Biomass__cytop'));
ecModel = setParam(ecModel, 'obj', params.bioRxn, 1);
sol = solveLP(ecModel);
bioRxnIdx = getIndexes(ecModel, params.bioRxn, 'rxns');
fprintf('Growth rate: %f /hour\n', sol.x(bioRxnIdx))

%38
model = loadConventionalGEM();
model = setParam(model, 'lb', params.c_source, -1000);
model = setParam(model, 'obj', params.bioRxn, 1);
sol = solveLP(model);
bioRxnIdx = getIndexes(model, params.bioRxn, 'rxns');
fprintf('Growth rate: %f /hour\n', sol.x(bioRxnIdx))

%39-42 
ecModel = setParam(ecModel, 'lb', 'prot_pool_exchange', -1000);

ecModel = setParam(ecModel, 'lb', 'EX_e_Biomass__cytop', -0.11);
ecModel = setParam(ecModel, 'obj', 'prot_pool_exchange', 1);
sol = solveLP(ecModel);
protPoolIdx = strcmp(ecModel.rxns, 'prot_pool_exchange');
fprintf('Protein pool usage is: %.0f mg/gDCW\n', abs(sol.x(protPoolIdx)))

ecModel = setParam(ecModel, 'lb', protPoolIdx, sol.x(protPoolIdx));
ecModel = setParam(ecModel, 'EX_e_Biomass__cytop', 0);
ecModel = setParam(ecModel, 'EX_e_Biomass__cytop', 1);

%43-44
ecModel = setProtPoolSize(ecModel);
[ecModel, tunedKcats] = sensitivityTuning(ecModel);
struct2table(tunedKcats)

%45-52
rxnIdx = find(strcmp(kcatList_merged.rxns, 'r_0079'));
kcatList_merged.wildcardLvl(rxnIdx)
kcatList_merged.eccodes(rxnIdx)
kcatList_merged.origin(rxnIdx)

convKcat = 2.15;
convKcat = convKcat / 1000;
convKcat = convKcat / 60;

enzMW = ecModel.ec.mw(strcmp(ecModel.ec.enzymes, 'P38972'));
convKcat = convKcat * enzMW

ecModel = setKcatForReactions(ecModel, 'r_0079', 5.34);
ecModel = applyKcatConstraints(ecModel);

saveEcModel(ecModel);