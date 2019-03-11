function [Asorted,I1,I2] = sort_Adj_rel_labels(A,label1,label2)

if label1
    [~,I1] = sort(label1);
else
    I1 = 1:size(A,1);
end

if label2
    [~,I2] = sort(label2);
else
    I2 = 1:size(A,2);
end

Asorted = A(I1,I2);