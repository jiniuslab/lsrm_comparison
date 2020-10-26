using DelimitedFiles
#using Pkg
#Pkg.add("Distributions")
#Pkg.update("Distributions")
using Distributions
using Printf

fname = "data/item.txt"
dataset = readdlm(fname, ' ', Int64)
dataset_origin = readdlm(fname, ' ', Int64)

#index for NA values
tmp_idx=findall(x->x==999,dataset)
dataset_origin[tmp_idx]

ndim = 2
niter = 30000
nburn = 5000
nthin = 5
nprint = 100

pr_mean_beta = 0.0
pr_mean_theta = 0.0
pr_mean_z = 0.0
pr_mean_w = 0.0

prior_a = 0.001
prior_b = 0.001

#pr_sd_beta = 10.0
#pr_sd_theta = 10.0

jump_beta = 0.25
jump_theta = 1.0
jump_z = 0.50
jump_w = 0.25

old_sigma_w = 1.0
old_sigma_z = 1.0

nsamp = size(dataset,1)
nitem = size(dataset,2)

count_item = sum(dataset,dims=1)
count_samp = sum(dataset,dims=2)

old_beta = randn(nitem)
old_theta = randn(nsamp)
old_z = randn(nsamp,ndim)
old_w = randn(nitem,ndim)

new_beta = zeros(nitem)
new_theta = zeros(nsamp)
new_z = zeros(nsamp,ndim)
new_w = zeros(nitem,ndim)

samp_like = zeros(nitem)
old_samp_dist = zeros(nitem)
new_samp_dist = zeros(nitem)

item_like = zeros(nsamp)
old_item_dist = zeros(nsamp)
new_item_dist = zeros(nsamp)

nmcmc = Int64((niter-nburn)/nthin)
sample_beta = zeros(nmcmc,nitem)
sample_theta = zeros(nmcmc,nsamp)
sample_z = zeros(nmcmc,nsamp,ndim)
sample_w = zeros(nmcmc,nitem,ndim)
sample_var_beta  = zeros(nmcmc)
sample_var_theta = zeros(nmcmc)
sample_mle = zeros(nmcmc)

accept_beta = zeros(nitem)
accept_theta = zeros(nsamp)
accept_z = zeros(nsamp)
accept_w = zeros(nitem)

distance = zeros(nsamp,nitem)


sample_impute = zeros(nmcmc,size(tmp_idx,1))


for iter = 1:niter
    ##Imputation step============================================================================================================================================================================
    impute_idx=0
    for i = 1:nitem
        for k = 1:nsamp
            if dataset_origin[k,i] == 999 # if missing value
                impute_dist=0
                for j = 1:ndim
                    impute_dist += (old_z[k,j]-old_w[i,j])^2
                end                
                p_ki=1/(1+exp(-(old_beta[i]+old_theta[k]-impute_dist)))
                impute_value=rand(Binomial(1,p_ki))
                dataset[k,i]=impute_value
                
                if iter > nburn && iter % nthin == 0 #to check samples of imputed value
                    impute_idx=impute_idx+1
                    count = (Int64)((iter-nburn)/nthin)                
                    sample_impute[count,impute_idx]=impute_value
                end
            end
        end
    end
#     @printf "\n"
#     @printf "imputation done" 
#     @printf "\n"
    ##Estimation step============================================================================================================================================================================
    for i = 1:nitem
        for j = 1:ndim
            new_w[i,j] = rand(Normal(old_w[i,j],jump_w))
        end
        for k = 1:nsamp
            item_like[k] = old_item_dist[k] = new_item_dist[k] = 0.0
        end

        for k = 1:nsamp
            for j = 1:ndim
                old_item_dist[k] += (old_z[k,j]-old_w[i,j])^2
                new_item_dist[k] += (old_z[k,j]-new_w[i,j])^2
            end
            old_item_dist[k] = sqrt(old_item_dist[k])
            new_item_dist[k] = sqrt(new_item_dist[k])
            if  dataset[k,i] == 1
                item_like[k] -= -log(1.0 + exp(-(old_beta[i]+old_theta[k]-old_item_dist[k])))
                item_like[k] += -log(1.0 + exp(-(old_beta[i]+old_theta[k]-new_item_dist[k])))
            else
                item_like[k] -= -log(1.0 + exp(old_beta[i]+old_theta[k]-old_item_dist[k]))
                item_like[k] += -log(1.0 + exp(old_beta[i]+old_theta[k]-new_item_dist[k]))
            end
        end

        update_like_item = 0.0
        for k = 1:nsamp
            update_like_item += item_like[k]
        end
        num = den = 0.0
        for j = 1:ndim
            num += logpdf(Normal(pr_mean_w,old_sigma_w),new_w[i,j])
            den += logpdf(Normal(pr_mean_w,old_sigma_w),old_w[i,j])
        end
        ratio = update_like_item + (num - den)

        if ratio > 0.0
            accept = 1
        else
            un = rand()
            if log(un) < ratio
                accept = 1
            else
                accept = 0
            end
        end

        if accept == 1
            for j = 1:ndim
                old_w[i,j] = new_w[i,j]
            end
            accept_w[i] += 1.0 / niter
        else
            for j = 1:ndim
                new_w[i,j] = old_w[i,j]
            end
        end
    end

    for k = 1:nsamp
        for j = 1:ndim
            new_z[k,j] = rand(Normal(old_z[k,j],jump_z))
        end
        for i = 1:nitem
            samp_like[i] = old_samp_dist[i] = new_samp_dist[i] = 0.0
        end

        for i = 1:nitem
            for j = 1:ndim
                old_samp_dist[i] += (old_z[k,j]-old_w[i,j])^2
                new_samp_dist[i] += (new_z[k,j]-old_w[i,j])^2
            end
            old_samp_dist[i] = sqrt(old_samp_dist[i])
            new_samp_dist[i] = sqrt(new_samp_dist[i])
            if  dataset[k,i] == 1
                samp_like[i] -= -log(1.0 + exp(-(old_beta[i]+old_theta[k]-old_samp_dist[i])))
                samp_like[i] += -log(1.0 + exp(-(old_beta[i]+old_theta[k]-new_samp_dist[i])))
            else
                samp_like[i] -= -log(1.0 + exp(old_beta[i]+old_theta[k]-old_samp_dist[i]))
                samp_like[i] += -log(1.0 + exp(old_beta[i]+old_theta[k]-new_samp_dist[i]))
            end
        end

        update_like_samp = 0.0
        for i = 1:nitem
            update_like_samp += samp_like[i]
        end
        num = den = 0.0
        for j = 1:ndim
            num += logpdf(Normal(pr_mean_z,old_sigma_z),new_z[k,j])
            den += logpdf(Normal(pr_mean_z,old_sigma_z),old_z[k,j])
        end
        ratio = update_like_samp + (num - den)

        if ratio > 0.0
            accept = 1
        else
            un = rand()
            if log(un) < ratio
                accept = 1
            else
                accept = 0
            end
        end

        if accept == 1
            for j = 1:ndim
                old_z[k,j] = new_z[k,j]
            end
            accept_z[k] += 1.0 / niter
        else
            for j = 1:ndim
                new_z[k,j] = old_z[k,j]
            end
        end
    end

    mean_beta = mean(old_beta)
    post_a = prior_a + 0.5 * nitem
    post_b = prior_b
    for i = 1:nitem
        post_b += 0.5 * (old_beta[i]-mean_beta)^2
    end
    old_var_beta = rand(InverseGamma(post_a,post_b))
    pr_sd_beta = sqrt(old_var_beta)

    for i = 1:nitem
        new_beta[i] = rand(Normal(old_beta[i],jump_beta))
        old_like_beta = new_like_beta = 0.0
        for k = 1:nsamp
            beta_dist = 0.0
            for j = 1:ndim
                beta_dist += (old_z[k,j]-old_w[i,j])^2
            end
            beta_dist = sqrt(beta_dist)

            if dataset[k,i] == 1
                old_like_beta += -log(1.0 + exp(-(old_beta[i]+old_theta[k]-beta_dist)))
                new_like_beta += -log(1.0 + exp(-(new_beta[i]+old_theta[k]-beta_dist)))
            else
                old_like_beta += -log(1.0 + exp(old_beta[i]+old_theta[k]-beta_dist))
                new_like_beta += -log(1.0 + exp(new_beta[i]+old_theta[k]-beta_dist))
            end
        end

        num = logpdf(Normal(pr_mean_beta,pr_sd_beta),new_beta[i])
        den = logpdf(Normal(pr_mean_beta,pr_sd_beta),old_beta[i])
        ratio = (new_like_beta-old_like_beta) + (num-den)

        if ratio > 0.0
            accept = 1
        else
            un = rand()
            if log(un) < ratio
                accept = 1
            else
                accept = 0
            end
        end

        if accept == 1
            old_beta[i] = new_beta[i]
            accept_beta[i] += 1.0 / niter
        else
            new_beta[i] = old_beta[i]
        end
    end

    mean_theta = mean(old_theta)
    post_a = prior_a + 0.5 * nsamp
    post_b = prior_b
    for k = 1:nsamp
        post_b += 0.5 * (old_theta[k]-mean_theta)^2
    end
    old_var_theta = rand(InverseGamma(post_a,post_b))
    pr_sd_theta = sqrt(old_var_theta)

    for k = 1:nsamp
        new_theta[k] = rand(Normal(old_theta[k],jump_theta))
        old_like_theta = new_like_theta = 0.0
        for i = 1:nitem
            theta_dist = 0.0
            for j = 1:ndim
                theta_dist += (old_z[k,j]-old_w[i,j])^2
            end
            theta_dist = sqrt(theta_dist)

            if dataset[k,i] == 1
                old_like_theta += -log(1.0 + exp(-(old_beta[i]+old_theta[k]-theta_dist)))
                new_like_theta += -log(1.0 + exp(-(old_beta[i]+new_theta[k]-theta_dist)))
            else
                old_like_theta += -log(1.0 + exp(old_beta[i]+old_theta[k]-theta_dist))
                new_like_theta += -log(1.0 + exp(old_beta[i]+new_theta[k]-theta_dist))
            end
        end

        num = logpdf(Normal(pr_mean_theta,pr_sd_theta),new_theta[k])
        den = logpdf(Normal(pr_mean_theta,pr_sd_theta),old_theta[k])
        ratio = (new_like_theta-old_like_theta) + (num-den)

        if ratio > 0.0
            accept = 1
        else
            un = rand()
            if log(un) < ratio
                accept = 1
            else
                accept = 0
            end
        end

        if accept == 1
            old_theta[k] = new_theta[k]
            accept_theta[k] += 1.0 / niter
        else
            new_theta[k] = old_theta[k]
        end
    end

    if iter > nburn && iter % nthin == 0
        count = (Int64)((iter-nburn)/nthin)

        mle = 0.0
        for k = 1:nsamp
            for i = 1:nitem
                distance[k,i] = 0.0
                for j = 1:ndim
                    distance[k,i] += (old_z[k,j]-old_w[i,j])^2
                end
                distance[k,i] = sqrt(distance[k,i])
                if dataset[k,i] == 1
                    mle += -log(1.0 + exp(-(old_beta[i]+old_theta[k]-distance[k,i])))
                else
                    mle += -log(1.0 + exp(old_beta[i]+old_theta[k]-distance[k,i]))
                end
            end
        end
        for k = 1:nsamp
            mle += logpdf(Normal(pr_mean_theta,pr_sd_theta),old_theta[k])
        end
        for i = 1:nitem
            mle += logpdf(Normal(pr_mean_beta,pr_sd_beta),old_beta[i])
        end

        for k = 1:nsamp
            for j = 1:ndim
                sample_z[count,k,j] = old_z[k,j]
            end
        end

        for i = 1:nitem
            for j = 1:ndim
                sample_w[count,i,j] = old_w[i,j]
            end
        end

        for i = 1:nitem
            sample_beta[count,i] = old_beta[i]
        end

        for k = 1:nsamp
            sample_theta[count,k] = old_theta[k]
        end

        sample_var_beta[count] = old_var_beta
        sample_var_theta[count] = old_var_theta

        sample_mle[count] = mle
    end

    if iter % nprint == 0
        @printf "%5d\n" iter
        for i = 1:nitem
            @printf "% .3f " old_beta[i]
        end
        @printf "\n"
        for k = 1:nitem
            @printf "% .3f " old_theta[k]
        end
        @printf "\n"
        for j = 1:ndim
            @printf "% .3f " old_z[1,j]
        end
        for j = 1:ndim
            @printf "% .3f " old_w[1,j]
        end
        @printf "\n"
        @printf "% .3f " old_var_beta
        @printf "% .3f\n" old_var_theta
    end
end


writedlm("result/impute.txt", sample_impute)
writedlm("result/beta.txt", sample_beta)
writedlm("result/theta.txt", sample_theta)
for j = 1:ndim
    oname = string("result/z",j,".txt")
    writedlm(oname, sample_z[:,:,j])
end
for j = 1:ndim
    oname = string("result/w",j,".txt")
    writedlm(oname, sample_w[:,:,j])
end
writedlm("result/var_beta.txt", sample_var_beta)
writedlm("result/var_theta.txt", sample_var_theta)
writedlm("result/mle.txt", sample_mle)
writedlm("result/impute.txt", sample_impute)
writedlm("result/accept_beta.txt", accept_beta)
writedlm("result/accept_theta.txt", accept_theta)
writedlm("result/accept_z.txt", accept_z)
writedlm("result/accept_w.txt", accept_w)