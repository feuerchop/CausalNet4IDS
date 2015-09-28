%%
% classification on nsl_kdd99: class label -> 16 (buffer_overflow attack type)
addpath(genpath('./'));
load datasets/nsl_kdd99.mat;
% preprocessing
tr_all = nsl_kdd99_tr(:, 1:42);
ts_all = nsl_kdd99_ts(:, 1:42);
tr_bin_all = tr_all;
ts_bin_all = ts_all;
tr_bin_all(find(tr_bin_all(:, 42) == 16), 42)

%%
% test logit distribution
close all; clear all;
q = [0, 1, 2, 3, 4, 5];
w = [-1, 1, -2, 2, 2, 1];
b = [1, 2, 1, 0, -2, 0];
x = 0:0.5:10;
p_nom = exp(w' * x + repmat(b', 1, 21));
p_dnom = sum(exp(w' * x + repmat(b', 1, 21)), 1);
p = bsxfun(@ldivide, p_dnom, p_nom);
xlim([0,10]); ylim([0,1]);
hold on;

plot(x, p(1,:), '.-k');
plot(x, w(1)' * x + repmat(b(1)', 1, 21), '--k');
plot(x, p(2,:), '.-b');
plot(x, w(2)' * x + repmat(b(2)', 1, 21), '--b');

plot(x, p(3,:), '.-r');
plot(x, w(3)' * x + repmat(b(3)', 1, 21), '--r');

plot(x, p(4,:), '.-g');
plot(x, w(4)' * x + repmat(b(4)', 1, 21), '--g');

plot(x, p(5,:), '.-c');
plot(x, w(5)' * x + repmat(b(5)', 1, 21), '--c');

plot(x, p(6,:), '.-y');
plot(x, w(6)' * x + repmat(b(6)', 1, 21), '--y');

%%
% A test for MLE
gauss = @(x, m, sigma) (1/((2*pi)^(1/2)))*(1/sigma)*exp(-(x-m).^2/(2*sigma^2));
data = mvnrnd(1.55, 2.0, 500);
plot(data, normpdf(data, 1.55, 2.0), 'ro', 'DisplayName', 'original');
hold on;
phat = mle(data, 'pdf', gauss, 'start', [0, 2]);
plot(data, normpdf(data, phat(1), phat(2)), 'bs', 'DisplayName', 'MLE');
phat
%%
% Learn a causal network with PC or PICM-CBN on kdd99 data and draw the
% graphical structure
load('datasets/nsl_kdd99.mat');
feature_names = {
   'duration'
   'protocol_type'
   'service'
   'flag'
   'src_bytes'
   'dst_bytes'
   'land'
   'wrong_fragment'
   'urgent'
   'hot'
   'num_failed_logins'
   'logged_in'
   'num_compromised'
   'root_shell'
   'su_attempted'
   'num_root'
   'num_file_creations'
   'num_shells'
   'num_access_files'
   'num_outbound_cmds'
   'is_host_login'
   'is_guest_login'
   'count'
   'srv_count'
   'serror_rate'
   'srv_serror_rate'
   'rerror_rate'
   'srv_rerror_rate'
   'same_srv_rate'
   'diff_srv_rate'
   'srv_diff_host_rate'
   'dst_host_count'
   'dst_host_srv_count'
   'dst_host_same_srv_rate'
   'dst_host_diff_srv_rate'
   'dst_host_same_src_port_rate'
   'dst_host_srv_diff_host_rate'
   'dst_host_serror_rate'
   'dst_host_srv_serror_rate'
   'dst_host_rerror_rate'
   'dst_host_srv_rerror_rate'
   'attack_type'
   'diff_level'
   };
feature_idx = [1:42];
train = nsl_kdd99_tr(1:1000, feature_idx);
test = nsl_kdd99_ts(1:1000, feature_idx);
node_names = feature_names(feature_idx);
% C = corrcoef(train);
% nsamples = size(train, 1);
% alpha = 0.05;
% nsize = size(train, 2);
% max_fan_in = 2;
% pdag = learn_struct_pdag_pc('cond_indep_fisher_z', nsize, max_fan_in, C, nsamples, alpha);
[eq, rho, pcdag] = BNLearningByCopula(train, 0.13);
% lp = CopulaLogProbabilityPerInstance(train, test, 'gauss', rho);
graph_to_dot(pcdag, 'filename', './virtualize/cpl_ids.dot', 'node_label', node_names);
%%
% An example of BNT: http://bnt.googlecode.com/svn/trunk/docs/usage.html#cg_model
% define the indexes for nodes
F = 1; W = 2; E = 3; B = 4; C = 5; D = 6; Min = 7; Mout = 8; L = 9;
% node size of the net
n = 9;
% construct the net in matrix-like format
dag = zeros(n);
dag(F,E)=1;
dag(W,[E Min D]) = 1;
dag(E,D)=1;
dag(B,[C D])=1;
dag(D,[L Mout])=1;
dag(Min,Mout)=1;
% node sizes - all cts nodes are scalar, all discrete nodes are binary
ns = ones(1, n);
dnodes = [F W B];
cnodes = mysetdiff(1:n, dnodes);
ns(dnodes) = 2;
% create bnet structure
bnet = mk_bnet(dag, ns, 'discrete', dnodes);
% The parameters of the discrete nodes can be specified as follows.
bnet.CPD{B} = tabular_CPD(bnet, B, 'CPT', [0.85 0.15]); % 1=stable, 2=unstable
bnet.CPD{F} = tabular_CPD(bnet, F, 'CPT', [0.95 0.05]); % 1=intact, 2=defect
bnet.CPD{W} = tabular_CPD(bnet, W, 'CPT', [2/7 5/7]); % 1=industrial, 2=household
% The parameters of the continuous nodes can be specified as follows.
bnet.CPD{E} = gaussian_CPD(bnet, E, 'mean', [-3.9 -0.4 -3.2 -0.5], ...
   'cov', [0.00002 0.0001 0.00002 0.0001]);
bnet.CPD{D} = gaussian_CPD(bnet, D, 'mean', [6.5 6.0 7.5 7.0], ...
   'cov', [0.03 0.04 0.1 0.1], 'weights', [1 1 1 1]);
bnet.CPD{C} = gaussian_CPD(bnet, C, 'mean', [-2 -1], 'cov', [0.1 0.3]);
bnet.CPD{L} = gaussian_CPD(bnet, L, 'mean', 3, 'cov', 0.25, 'weights', -0.5);
bnet.CPD{Min} = gaussian_CPD(bnet, Min, 'mean', [0.5 -0.5], 'cov', [0.01 0.005]);
bnet.CPD{Mout} = gaussian_CPD(bnet, Mout, 'mean', 0, 'cov', 0.002, 'weights', [1 1]);
save('examples/cond_gaussian.bnet', 'bnet');
%%
p = 0.35;
v = mnrnd(1000, [0.15, 0.85]);
data = [1*ones(v(1), 1); 2*ones(v(2), 1)];
options = statset('Display','final');
obj = gmdistribution.fit(data,2,'Options',options);
