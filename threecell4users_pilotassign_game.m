% 3 cell 4 users/cell pilot sequence choice game
clear
clc
%% simualtion parameters radndomly generated
N = 8;
Ns = 500;
K = 3;%3 cells
I = 4;%4 users/cell
tau_p = 4;% number of available pilot sequence in each cell
Pp = 19+randi(21,I*K,1);% pilot power of all users ranging from 20 to 40
sigma = sqrt(db2lin(-174)*180000);%standard deviation of noise
gamma = 0.5.*eye(N);
cov_mat = eye(N);
AOA = 20;% the angle of arrival (assumed to be same for all users)
seq = ones(tau_p,1);
seqT = transpose(seq);

% pathloss randomization
for i = 1:K
    for j = 1:K
        if j == i
            PL(I*(i-1)+1:I*i,j) = 39+randi(21,I,1);%pathloss within the own cell ranging from 40 to 60
        else
            PL(I*(i-1)+1:I*i,j) = 59+randi(21,I,1);%pathloss between the cells ranging from 60 to 80
        end
    end
end

%% random pilot sequence allocation
S_pool = eye(tau_p);
for i = 1:K
    choice = randperm(tau_p,I);
    for j = 1:I
        S_user(j,:) = S_pool(choice(j),:);
        SS((i-1)*I+j,:) = S_user(j,:);
    end
    S(i,:) = reshape(S_user',1,tau_p*I);
end
%% channel estimation simulation
for i = 1:K
    for j = 1:I
        MS_index = (i-1)*I+j;
        pilot_usage = SS(:,SS(MS_index,:)==1);
        num_interfer = sum(pilot_usage)-1;%number of MS that interfering
        PL_cell = PL(:,i);
        Pp_MS = Pp(MS_index);% pilot power of the tagged user
        PL_MS = PL_cell(MS_index);% pathloss of the tagged user
        Pp_nMS = Pp;
        Pp_nMS((i-1)*I+1:i*I,1) = 0;
        PL_nMS = PL_cell;
        PL_nMS((i-1)*I+1:i*I,1) = 0;
        alpha_MS = db2lin(-PL_MS);
        h_MS = mat2cell(Rchanneladv(Ns,N,AOA,gamma),N,ones(1,Ns));% channel of tagged user
        yp_MS = cellfun(@(x) x*seqT*alpha_MS*sqrt(Pp_MS),h_MS,'un',0);
        for z = 1:Ns
            n_p{z} = Gmultinoise(N,tau_p,sigma);% pilot noise
        end
        if num_interfer ~= 0
            Pp_interfer = Pp_nMS.*pilot_usage;
            PL_interfer = PL_nMS.*pilot_usage;
            Pp_interfer = Pp_interfer(Pp_interfer ~= 0);
            PL_interfer = PL_interfer(PL_interfer ~= 0);
            alpha_interfer = db2lin(-PL_interfer);
            yp_interfer = mat2cell(zeros(N,tau_p),N,tau_p);
            yp_interfer = repmat(yp_interfer,1,Ns);
            for v = 1:num_interfer
                h_interfer = transpose(Gchannel(Ns,N,cov_mat));
                h_interfer = mat2cell(h_interfer,N,ones(1,Ns));
                yp_interfer_temp = cellfun(@(x) x*seqT*alpha_interfer(v)*sqrt(Pp_interfer(v)),h_interfer,...
                    'un',0);
                yp_interfer = cellfun(@plus,yp_interfer,yp_interfer_temp,...
                    'un',0);
            end
            yp_nnp = cellfun(@plus,yp_MS,yp_interfer,'un',0);
            yp = cellfun(@plus,yp_nnp,n_p,'un',0);
        else
            yp = cellfun(@plus,yp_MS,n_p,'un',0);
        end
        s_MS_key = 1/(alpha_MS*sqrt(Pp_MS)).*conj(seq)*inv((transpose(seq)*conj(seq)));
        h_MS_e = cellfun(@(x) x*s_MS_key,yp,'un',0);
        h_MS_e = cell2mat(h_MS_e);
        est_error = h_MS_e-cell2mat(h_MS);
        normsqu_error = sum(abs(est_error).^2);
        normsqu_h1 = sum(abs(cell2mat(h_MS)).^2);
        normSE_random(MS_index,:) = normsqu_error./normsqu_h1;
    end
end

%%%plotting
% for i = 1:I*K
%     fig = cdfplot(lin2db(normSE_random(i,:)));
%     set(fig,'Color',color_pool(i,:),'Linewidth',1.5);
%     hold on
% end

% cdfplot(lin2db(reshape(normSE_random,1,[])));hold on

%% pilot assignment game
Nit = 50;
S_game{1} = S;% random assignment

%interference matrix costruction
alpha = db2lin(-PL);
for i = 1:K
    PL_self(i,:) = PL((i-1)*I+1:i*I,i)';
    Pp_game(i,:) = Pp((i-1)*I+1:i*I);
end
LS_par = Pp_game.*(db2lin(-PL_self).^2);

alpha_game = reshape(alpha',K,I,[]);
for i = 1:K
    alpha_i = alpha_game(:,:,i);
    Pp_i = Pp_game(i,:);
    Pp_ni = Pp_game;
    Pp_ni(i,:) = [];
    alpha_game_i = alpha_game(:,:,i);
    alpha_game_i(i,:) = [];
    creat_interfer(:,:,i) = repmat(Pp_i,K-1,1).*(alpha_game_i).^2;
    alpha_game_ni = alpha_game(i,:,:);
    alpha_game_ni = reshape(permute(alpha_game_ni,[1 3 2]),K,[]);
    alpha_game_ni(i,:) = [];
    suff_interfer(:,:,i) = Pp_ni.*(alpha_game_ni.^2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% game asynchronous
for i = 2:20
    i
    S_game{i} = S_game{i-1};
    for j = 1:K
        S_game_other = S_game{i};
        S_game_other(j,:) = [];
        LS_omega_par = LS_par;
        LS_omega_par(j,:) = [];
        LS_delta_par = LS_par(j,:);
        seq_choice = BRPilot(S_game_other,suff_interfer(:,:,j),creat_interfer(:,:,j),tau_p,LS_omega_par,LS_delta_par);
        S_game{i}(j,:) = seq_choice;
    end
    if isequal(S_game{i},S_game{i-1})
        break% another way to terminate the game
    end
    clc
end

S_game = S_game{end};
SS_game = reshape(S_game',tau_p,[])';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% channel estimation simulation
for i = 1:K
    for j = 1:I
        MS_index = (i-1)*I+j;
        pilot_usage = SS_game(:,SS_game(MS_index,:)==1);
        num_interfer = sum(pilot_usage)-1;%number of MS that interfering
        PL_cell = PL(:,i);
        Pp_MS = Pp(MS_index);% pilot power of the tagged user
        PL_MS = PL_cell(MS_index);% pathloss of the tagged user
        Pp_nMS = Pp;
        Pp_nMS((i-1)*I+1:i*I,1) = 0;
        PL_nMS = PL_cell;
        PL_nMS((i-1)*I+1:i*I,1) = 0;
        alpha_MS = db2lin(-PL_MS);
        h_MS = mat2cell(Rchanneladv(Ns,N,AOA,gamma),N,ones(1,Ns));% channel of tagged user
        yp_MS = cellfun(@(x) x*seqT*alpha_MS*sqrt(Pp_MS),h_MS,'un',0);
        for z = 1:Ns
            n_p{z} = Gmultinoise(N,tau_p,sigma);% pilot noise
        end
        if num_interfer ~= 0
            Pp_interfer = Pp_nMS.*pilot_usage;
            PL_interfer = PL_nMS.*pilot_usage;
            Pp_interfer = Pp_interfer(Pp_interfer ~= 0);
            PL_interfer = PL_interfer(PL_interfer ~= 0);
            alpha_interfer = db2lin(-PL_interfer);
            yp_interfer = mat2cell(zeros(N,tau_p),N,tau_p);
            yp_interfer = repmat(yp_interfer,1,Ns);
            for v = 1:num_interfer
                h_interfer = transpose(Gchannel(Ns,N,cov_mat));
                h_interfer = mat2cell(h_interfer,N,ones(1,Ns));
                yp_interfer_temp = cellfun(@(x) x*seqT*alpha_interfer(v)*sqrt(Pp_interfer(v)),h_interfer,...
                    'un',0);
                yp_interfer = cellfun(@plus,yp_interfer,yp_interfer_temp,...
                    'un',0);
            end
            yp_nnp = cellfun(@plus,yp_MS,yp_interfer,'un',0);
            yp = cellfun(@plus,yp_nnp,n_p,'un',0);
        else
            yp = cellfun(@plus,yp_MS,n_p,'un',0);
        end
        s_MS_key = 1/(alpha_MS*sqrt(Pp_MS)).*conj(seq)*inv((transpose(seq)*conj(seq)));
        h_MS_e = cellfun(@(x) x*s_MS_key,yp,'un',0);
        h_MS_e = cell2mat(h_MS_e);
        est_error = h_MS_e-cell2mat(h_MS);
        normsqu_error = sum(abs(est_error).^2);
        normsqu_h1 = sum(abs(cell2mat(h_MS)).^2);
        normSE_game(MS_index,:) = normsqu_error./normsqu_h1;
    end
end



%%% plotting
% for i = 1:I*K
%     fig = cdfplot(lin2db(normSE_random(i,:)));
%     set(fig,'LineStyle', '--','Color',color_pool(i,:),'Linewidth',1.5);
%     hold on
% end
% cdfplot(lin2db(reshape(normSE_game,1,[])));hold on

% legend('MS1-1','MS2-1','MS3-1','MS4-1',...
% 'MS1-2','MS2-2','MS3-2','MS4-2',...
% 'MS1-3','MS2-3','MS3-3','MS4-3');
% xlabel('MSE of channel estimation error (dB)')
% ylabel('Individual CDF')
% title('Individual MSE of channel estimation error when using random and game pilot assignment')
%% post processing

normSE_random_db = lin2db(normSE_random);
normSE_game_db = lin2db(normSE_game);

load('SErandom.mat')
load('SEgame.mat')
normSE_random_mat = [normSE_random_mat,normSE_random_db];
normSE_game_mat = [normSE_game_mat,normSE_game_db];

save('SErandom.mat','normSE_random_mat')
save('SEgame.mat','normSE_game_mat')



% MSE_random_db = sort(mean(transpose(lin2db(normSE_random))));
% MSE_game_db = sort(mean(transpose(lin2db(normSE_game))));
% 
% load('MSErandom.mat')
% load('MSEgame.mat')
% MSE_random_mat = [MSE_random_mat;MSE_random_db];
% MSE_game_mat = [MSE_game_mat;MSE_game_db];
% 
% save('MSErandom.mat','MSE_random_mat')
% save('MSEgame.mat','MSE_game_mat')
