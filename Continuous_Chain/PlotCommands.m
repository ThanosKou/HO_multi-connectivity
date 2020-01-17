figure()
semilogy(lambda_BS,P_OS(:,1,1,2,1),'-rx')
grid on;
hold on;
semilogy(lambda_BS,P_OS(:,1,1,2,2),'-r+')

semilogy(lambda_BS,P_OS(:,1,1,5,1),'-bx')
hold on;
semilogy(lambda_BS,P_OS(:,1,1,5,2),'-b+')
