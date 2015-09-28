clc;close all;clear all;
% load normal data first
normal_xxx = load('datasets/kdd99_normal_10_percent.mat');
normal_data = normal_xxx.kdd_normal_data;
tbl_names = normal_xxx.tbl_names;
abnormal_back_xxx = load('datasets/kdd99_back_10_percent.mat');
back_data = abnormal_back_xxx.kdd_back_data;
% random sampling 2000 of each
sample_size = 1000;
xtr = [ normal_data(randsample(length(normal_data), sample_size), :); ...
        back_data(randsample(length(back_data), sample_size), :) ];
% labels = xtr(:, end);
% xtr = xtr(:, 1:end-1);
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
node_names = feature_names(feature_idx);
[eq, rho, pcdag] = BNLearningByCopula(xtr, 0.15);
graph_to_dot(pcdag, 'filename', './virtualize/on_back_atk.dot', 'node_label', node_names);
