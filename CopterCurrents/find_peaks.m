function peak_flag = find_peaks(A)


peak_flag = zeros(size(A));
% peak_x = zeros(size(A));
% peak_y = peak_x;

% peak_x(:,2:end-1) = (A(:,3:end) - 2*A(:,2:end-1)+A(:,1:end-2))<0 &...
%     A(:,2:end-1)>A(:,1:end-2) & A(:,2:end-1)>A(:,3:end);
% peak_y(2:end-1,:) = (A(3:end,:) - 2*A(2:end-1,:)+A(1:end-2,:))<0 &...
%     A(2:end-1,:)>A(1:end-2,:) & A(2:end-1,:)>A(3:end,:);
% peak_flag = peak_x.*peak_y;

%Corners (3 nearest neighbors)
peak_flag(1,1) = A(1,1)>A(1,2) & A(1,1)>A(2,1) & A(1,1)>A(2,2);
peak_flag(1,end) = A(1,end)>A(1,end-1) & A(1,end)>A(2,end) & A(1,end)>A(2,end-1);
peak_flag(end,end) = A(end,end)>A(end,end-1) & A(end,end)>A(end-1,end) & A(end,end)>A(end-1,end-1);
peak_flag(end,1) = A(end,1)>A(end-1,1) & A(end,1)>A(end-1,2) & A(end,1)>A(end,2);
peak_flag = logical(peak_flag);


%Edges (5 nearest neighors)
peak_flag(1,2:end-1) = A(1,2:end-1)>A(1,1:end-2) & A(1,2:end-1)>A(1,3:end) &...
    A(1,2:end-1)>A(2,1:end-2) & A(1,2:end-1)>A(2,2:end-1) & A(1,2:end-1)>A(2,3:end);%Top edge

peak_flag(end,2:end-1) = A(end,2:end-1)>A(end,1:end-2) & A(end,2:end-1)>A(end,3:end) &...
    A(end,2:end-1)>A(end-1,1:end-2) & A(end,2:end-1)>A(end-1,2:end-1) & A(end,2:end-1)>A(end-1,3:end);%Bottom edge

peak_flag(2:end-1,1) = A(2:end-1,1)>A(1:end-2,1) & A(2:end-1,1)>A(3:end,1) &...
    A(2:end-1,1)>A(1:end-2,2) & A(2:end-1,1)>A(2:end-1,2) &  A(2:end-1,1)>A(3:end,2);%Left edge

peak_flag(2:end-1,end) = A(2:end-1,end)>A(1:end-2,end) & A(2:end-1,end)>A(3:end,end) & ...
    A(2:end-1,end)>A(1:end-2,end-1) & A(2:end-1,end)>A(2:end-1,end-1) & A(2:end-1,end)>A(3:end,end-1);%Right edge

%Interior region (8 nearest neighbors)
peak_flag(2:end-1,2:end-1) = ...
    A(2:end-1,2:end-1)>A(1:end-2,1:end-2) &...
    A(2:end-1,2:end-1)>A(2:end-1,1:end-2) &...
    A(2:end-1,2:end-1)>A(3:end,1:end-2) &...
    A(2:end-1,2:end-1)>A(3:end,2:end-1) &...
    A(2:end-1,2:end-1)>A(3:end-0,3:end) &...
    A(2:end-1,2:end-1)>A(2:end-1,3:end) &...
    A(2:end-1,2:end-1)>A(1:end-2,3:end) &...
    A(2:end-1,2:end-1)>A(1:end-2,2:end-1);



end