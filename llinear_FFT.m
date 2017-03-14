function yout = llinear_FFT(xin,yin,bw)
% -------------------------------------------------------------------------
% xin: 1 x d cell array; each cell contains M1 x 1 vector
% yin: M1 x M2 x ... x Md array
% bw: 1 x d array
% -------------------------------------------------------------------------
    n = numel(yin);
    d = length(xin);
    M = zeros(1,d); a = zeros(1,d); b = zeros(1,d); L = zeros(1,d);
    for i = 1:d
        M(i) = length(xin{i});
        a(i) = min(xin{i}); b(i) = max(xin{i});
        L(i) = min(floor(bw(i)*(M(i)-1)),(M(i)-1));
        ksize = L(1); m = M(1); % assume diagnal
        K = @(u)0.75*(1-u.^2).*(abs(u)<1);
        l = (-ksize:ksize)';
        if isa(yin,'single')
            l = single(l);
        end
        x{i} = l*(b(i)-a(i))/(M(i)-1);
        k{i} = K(x{i}/bw(i))/bw(i)/M(i);
    end
    if d>1
        kk = ndgrid_cell(k);
        kernel = ones(size(kk{1}));
        for i = 1:d
            kernel = kernel.*kk{i};
        end
        xx = ndgrid_cell(x);
    else
        kernel = k{1};
        xx = x;
    end
    xx = [ones(size(xx{1}),class(yin)),xx];
    S = zeros(d+1,d+1,n,class(yin));
    T = zeros(d+1,n,class(yin));
    for i = 1:(d+1)
        for j = i:(d+1)
            S(i,j,:) = reshape(fastconvn(ones(size(yin),class(yin)),kernel.*xx{i}.*xx{j}),1,1,n);
            S(j,i,:) = S(i,j,:);
        end
    end
    for i = 1:(d+1)
        T(i,:) = reshape(fastconvn(yin,kernel.*xx{i}),1,n);
    end
    % yout = T(1,:)./reshape(S(1,1,:),1,n);
    yout = zeros(n,1,class(yin));
    for i = 1:n
        b = S(:,:,i)\T(:,i);
        yout(i) = b(1);
    end
    