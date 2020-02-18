outputs = dir(['output','*']);
num_files = length(outputs);


discovery = [1 5 20 200 1000]*10^(-3);
preparation = [10 20]*10^(-3);
densityBL = [0.01 0.1];
densityBS = [200 300 400 500]*10^(-6);
connectivity = [1 2 3 4];

results_array = zeros(length(discovery),length(preparation),length(densityBS),length(connectivity),length(densityBL),num_files);
 

for ii=1:num_files
    ii
    aa = load(outputs(ii).name);
    results_array(:,:,:,:,:,ii) = aa.finaldata;
end

final_results = mean(results_array,6);

save(strcat('finalresults.mat'),'final_results','discovery','preparation','densityBL','densityBS','connectivity')
