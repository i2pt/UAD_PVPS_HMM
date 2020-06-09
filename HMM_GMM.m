function [y, m0, v0] = HMM_GMM(X,imax)
%% Fixed Parameters
T = length(X); % Samples
m = 2; % Classes
mylim = [0 0];
%% Initial Values of Parameters to Estimate
del = [0.5 0.5]; % initial probability of state
A = [0.5 0.5; 0.5 0.5]; % state transition matrix
m0 = mylim; % observation mean
v0 = [std(X) 0.5*max(abs(X))];
%% Simulation Parameters
LL = -Inf; % log-likelihood
for itr = 1:imax
    %%
    alpha = zeros(T,m); % forward probability
    beta = zeros(T,m); % backward probability
    c = zeros(T,1); % normalizing factor
    %% Scaled Forward Probabilities
    alpha_1 = del.*normpdf(X(1),m0,v0); % Unscaled
    c(1) = 1/sum(alpha_1); % c_1
    alpha(1,:) = alpha_1.*c(1); % Scaled
    
    for t = 2:T
        abar = zeros(1,m); % Intermediate
        for i = 1:m
            for j = 1:m
                abar(i) = abar(i) + alpha(t-1,j)*A(j,i)*normpdf(X(t),m0(i),v0(i));
            end
        end
        c(t) = 1/sum(abar);
        alpha(t,:) = abar.*c(t); % Scaled
    end
    %% Scaled Backward Probabilities
    beta(T,:) = c(T)*ones(1,m);
    
    for t = (T-1):-1:1
        bbar = zeros(1,m);
        for i = 1:m
            for j = 1:m
                bbar(i) = bbar(i) + A(i,j)*normpdf(X(t+1),m0(j),v0(j))*beta(t+1,j);
            end
        end
        beta(t,:) = bbar.*c(t);
    end
    %%
    d = zeros(T,m);
    dd = zeros(m,m,T);    
    d(T,:) = alpha(T,:);
%%    
    for t = 1:(T-1)
        for i = 1:m
            d(t,i) = 0;
            for j = 1:m
                dd(i,j,t) = alpha(t,i)*A(i,j)*normpdf(X(t+1),m0(j),v0(j))*beta(t+1,j);
                d(t,i) = d(t,i) + dd(i,j,t);
            end
        end
    end
    %% Estimate
    del = d(1,:)/sum(d(1,:));
    
    for i = 1:m
        den = sum(d(1:(T-1),i));
        for j = 1:m
            A(i,j) = sum(dd(i,j,1:(T-1)))/den;
        end
    end
    
    for i = 1:m
        den = sum(d(:,i));
        m0(i) = sum(X.*d(:,i)')/den;
        v0(i) = sqrt(sum(((X-m0(i)).^2).*d(:,i)')/den);
    end
    %%  Log-likelihood    
    LL(itr+1) = -sum(log(c));
    
    if LL(itr+1) < LL(itr)
        disp('mahayat : LL cannot decrease < SOMETHING WRONG >');
        break;
    elseif LL(itr+1)-LL(itr) < 1e-5
        disp('mahayat : LL not increasing');
        break;
    elseif isnan(LL(itr+1))
        disp('mahayat : NaN value < TAKE ACTION > ');
        break;
    else
        disp(['mahayat : Iteration #',num2str(itr),' :: LL = ',num2str(LL(itr+1))]);
    end
end
%% Distribution Parameters
[~, anomaly_class_index] = max(v0);
%% Viterbi Algorithm
[~, temp2] = max(d,[],2);
y = temp2-1;
end