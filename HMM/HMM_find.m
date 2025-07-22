% Arrangment by sorting
posimu = zeros(1,Q);

mu_result = zeros(1,Q);
mixmat_result = zeros(1,Q);
prior_result = zeros(1,Q);
sigsize = size(sigma);
sigma_result = 0;
transmat_result = zeros(Q,Q);


arrangemu = sort(mu,'descend');
   
for i=1:Q
    posimu(i)= find(arrangemu(i) == mu);
end

for i = 1:Q
    mu_result(i) = mu(posimu(i));
    mixmat_result(i) = mixmat(posimu(i));
    prior_result(i) = prior(posimu(i));
    sigma_result(:,:,i) = sigma(:,:,posimu(i));
    for j =1:Q
        transmat_result(i,j)=transmat(posimu(i),posimu(j))*fps;
    end
end

HMMtrace = {};
Gaussprob = {};
for  i = 1:size(data,2)
    datab = data{i};
    B = mixgauss_prob(datab,mu_result,sigma_result,mixmat_result);
    path = viterbi_path(prior_result,  transmat_result, B);
    
    Gaussprob{i} = B;
    
    for  j = 1:length(path)
        HMMtrace{i}(j) =  mu_result(path(j));
    end
end

