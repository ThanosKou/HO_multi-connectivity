outputs_1 = dir(['blockages1/blockages','*']);
outputs_2 = dir(['blockages2/blockages','*']);
outputs = [outputs_1; outputs_2];
num_files = length(outputs);


discovery = [1 5 20 200 1000]*10^(-3);
preparation = [10 20]*10^(-3);
densityBL = [0.01 0.1];
densityBS = [200 300 400 500]*10^(-6);
connectivity = [1 2 3 4];

% blockages = cell(length(discovery),length(preparation),length(densityBS),length(connectivity),length(densityBL),num_files);
mean_blockages = zeros(length(discovery),length(preparation),length(densityBS),length(connectivity),length(densityBL));
% all the u0s initalized as zeros
k = zeros(length(discovery),length(preparation),length(densityBS),length(connectivity),length(densityBL));
% all the ks initalized as 0; when got  first value will become one.

for ii=1:num_files
    ii
    aa = load([outputs(ii).folder,'/',outputs(ii).name]);
    for idxD = 1:length(discovery)
        for idxP = 1:length(preparation)
            for idxBL = 1:length(densityBL)
                for idxBS = 1:length(densityBS)
                    for idxK = 1:length(connectivity)
                        if isempty(aa.blockageDurations{idxBS,idxK,idxBL}{idxD,idxP})
                            % do nothing
                        else
                            k(idxD,idxP,idxBS,idxK,idxBL) = k(idxD,idxP,idxBS,idxK,idxBL) + 1;
                            xk = mean(aa.blockageDurations{idxBS,idxK,idxBL}{idxD,idxP});
                            uk_1 = mean_blockages(idxD,idxP,idxBS,idxK,idxBL);
                            mean_blockages(idxD,idxP,idxBS,idxK,idxBL) = uk_1 + (1/k(idxD,idxP,idxBS,idxK,idxBL)) * (xk - uk_1);
                        end
                    end
                end
            end
        end
    end
end



Description = 'The blockages is a cell array of 6 dimension where 1st dimension is index for discovery, 2nd for preperation, 3rd for Base Station Desity, 4th is Connectivity, 5th is  blocker density, 6th is the individual files. Each element of this cell is an array where each blockage event duration is recorded. The mean blockages is 5 dimensional and the average blockage duration over all the files recorded according to the first 5 parameters of aforomentioned cell.'
save(strcat('BlockageData_RunAvg.mat'),'mean_blockages','discovery','preparation','densityBL','densityBS','connectivity','Description')
% save(strcat('AllBlockageData.mat'),'blockages','mean_blockages','discovery','preparation','densityBL','densityBS','connectivity','Description' , '-v7.3')







