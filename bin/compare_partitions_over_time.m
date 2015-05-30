function [ vi_mat ] = compare_partitions_over_time( C1, C2, T)
%COMPARE_PARTITIONS_OVER_TIME Compute a VI(t,t') matrix for two given
%partition sequences and a Time vector

nr_steps = length(T);
vi_mat = zeros(nr_steps);

for i = 1:nr_steps
    for j= i:nr_steps
        vi_mat(i,j) = varinfo([C1(:,i),C2(:,j)]');       
    end
end

vi_mat = vi_mat+vi_mat' - diag(diag(vi_mat));

surf(T,T,vi_mat);
xlabel('Time t_1');
ylabel('Time t_2');
zlabel('Variation of information');

end

