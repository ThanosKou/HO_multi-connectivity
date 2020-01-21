figure()
semilogy(lambda_BS,P_OS(:,1,2,3,1),'--bx')
grid on;
hold on;
semilogy(lambda_BS,P_OS(:,2,2,3,1),'--b^')
semilogy(lambda_BS,P_OS(:,3,2,3,1),'--bo')
semilogy(lambda_BS,P_OS(:,4,2,3,1),'--bp')
legend('1 Connectivity','2 Connectivity','3 Connectivity','4 Connectivity')
title('20 ms Association and 20 ms Discovery Time')

figure()
semilogy(lambda_BS,P_OS(:,1,2,1,2),'-rx')
grid on;
hold on;
semilogy(lambda_BS,P_OS(:,1,2,2,2),'-r+')
semilogy(lambda_BS,P_OS(:,1,2,3,2),'-r*')
semilogy(lambda_BS,P_OS(:,1,2,4,2),'-rd')
semilogy(lambda_BS,P_OS(:,1,2,5,2),'-r^')

% semilogy(lambda_BS,P_OS(:,1,2,1,2),'-bx')
% semilogy(lambda_BS,P_OS(:,1,2,2,2),'-b+')
% semilogy(lambda_BS,P_OS(:,1,2,3,2),'-b*')
% semilogy(lambda_BS,P_OS(:,1,2,4,2),'-bd')
% semilogy(lambda_BS,P_OS(:,1,2,5,2),'-b^')

semilogy(lambda_BS,P_OS(:,2,2,1,2),'-gx')
semilogy(lambda_BS,P_OS(:,2,2,2,2),'-g+')
semilogy(lambda_BS,P_OS(:,2,2,3,2),'-g*')
semilogy(lambda_BS,P_OS(:,2,2,4,2),'-gd')
semilogy(lambda_BS,P_OS(:,2,2,5,2),'-g^')

semilogy(lambda_BS,P_OS(:,3,2,1,2),'-mx')
semilogy(lambda_BS,P_OS(:,3,2,2,2),'-m+')
semilogy(lambda_BS,P_OS(:,3,2,3,2),'-m*')
semilogy(lambda_BS,P_OS(:,3,2,4,2),'-md')
semilogy(lambda_BS,P_OS(:,3,2,5,2),'-m^')

semilogy(lambda_BS,P_OS(:,4,2,1,2),'-bx')
semilogy(lambda_BS,P_OS(:,4,2,2,2),'-b+')
semilogy(lambda_BS,P_OS(:,4,2,3,2),'-b*')
semilogy(lambda_BS,P_OS(:,4,2,4,2),'-bd')
semilogy(lambda_BS,P_OS(:,4,2,5,2),'-b^')




figure()
semilogy(K_list,P_OS(3,:,2,1,2),'-rx')
hold on;
semilogy(K_list,P_OS(3,:,2,2,2),'-r+')
semilogy(K_list,P_OS(3,:,2,3,2),'-r*')
semilogy(K_list,P_OS(3,:,2,4,2),'-rd')
semilogy(K_list,P_OS(3,:,2,5,2),'-r^')
grid on;

semilogy(K_list,P_OS(3,:,2,1,1),'-bx')
semilogy(K_list,P_OS(3,:,2,2,1),'-b+')
semilogy(K_list,P_OS(3,:,2,3,1),'-b*')
semilogy(K_list,P_OS(3,:,2,4,1),'-bd')
semilogy(K_list,P_OS(3,:,2,5,1),'-b^')
