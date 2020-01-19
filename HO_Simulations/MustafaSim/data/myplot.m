
discovery = [1 5 20 200 1000]*10^(-3);
preparation = [10 20]*10^(-3);
densityBL = [0.01 0.1];
densityBS = [200 300 400 500]*10^(-6);
connectivity = [1 2 3 4];

load('finalresults.mat')
figure()
semilogy(densityBS, reshape(final_results(1,1,:,1,1),[],1),'-+r')
hold on;
semilogy(densityBS, reshape(final_results(2,1,:,1,1),[],1),'-xr')
semilogy(densityBS, reshape(final_results(3,1,:,1,1),[],1),'-*r')
semilogy(densityBS, reshape(final_results(4,1,:,1,1),[],1),'-or')
semilogy(densityBS, reshape(final_results(5,1,:,1,1),[],1),'-sr')


semilogy(densityBS, reshape(final_results(1,2,:,1,1),[],1),'-+b')
semilogy(densityBS, reshape(final_results(2,2,:,1,1),[],1),'-xb')
semilogy(densityBS, reshape(final_results(3,2,:,1,1),[],1),'-*b')
semilogy(densityBS, reshape(final_results(4,2,:,1,1),[],1),'-ob')
semilogy(densityBS, reshape(final_results(5,2,:,1,1),[],1),'-sb')


figure()
semilogy(densityBS, reshape(final_results(1,1,:,2,1),[],1),'-+r')
hold on;
semilogy(densityBS, reshape(final_results(2,1,:,2,1),[],1),'-xr')
semilogy(densityBS, reshape(final_results(3,1,:,2,1),[],1),'-*r')
semilogy(densityBS, reshape(final_results(4,1,:,2,1),[],1),'-or')
semilogy(densityBS, reshape(final_results(5,1,:,2,1),[],1),'-sr')


semilogy(densityBS, reshape(final_results(1,2,:,2,1),[],1),'-+b')
semilogy(densityBS, reshape(final_results(2,2,:,2,1),[],1),'-xb')
semilogy(densityBS, reshape(final_results(3,2,:,2,1),[],1),'-*b')
semilogy(densityBS, reshape(final_results(4,2,:,2,1),[],1),'-ob')
semilogy(densityBS, reshape(final_results(5,2,:,2,1),[],1),'-sb')



figure()
semilogy(connectivity, reshape(final_results(1,1,3,:,1),[],1),'-+r')
hold on;
semilogy(connectivity, reshape(final_results(2,1,3,:,1),[],1),'-xr')
semilogy(connectivity, reshape(final_results(3,1,3,:,1),[],1),'-*r')
semilogy(connectivity, reshape(final_results(4,1,3,:,1),[],1),'-or')
semilogy(connectivity, reshape(final_results(5,1,3,:,1),[],1),'-sr')

semilogy(connectivity, reshape(final_results(1,2,3,:,1),[],1),'-+b')
semilogy(connectivity, reshape(final_results(2,2,3,:,1),[],1),'-xb')
semilogy(connectivity, reshape(final_results(3,2,3,:,1),[],1),'-*b')
semilogy(connectivity, reshape(final_results(4,2,3,:,1),[],1),'-ob')
semilogy(connectivity, reshape(final_results(5,2,3,:,1),[],1),'-sb')