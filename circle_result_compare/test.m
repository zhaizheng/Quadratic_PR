    
marker1 = {'-.s','-.*','-.o','-->','--d'};
sigma = 0.1:0.1:0.9;
for k = 1:5
    semilogy(sigma,result(k,:),marker1{k},'linewidth',1.2); 
    hold on; 
end