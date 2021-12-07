function times =  ms4_detect_on_channel(data, detect_threshold,detect_interval,detect_sign,margin)
    % hanlin implmementation of matlab version of mountainsort detection
    % step, follow code from ms4. 
    %# Adjust the data to accommodate the detect_sign
    %# After this adjustment, we only need to look for positive peaks
    if detect_sign<0
        data=data*(-1);
    elseif detect_sign==0
        data=abs(data);
    elseif detect_sign>0
    end
    N=numel(data);
    S2=floor(N/detect_interval);
    N2=S2*detect_interval;
    data2=reshape(data(1:N2),detect_interval,S2);
    [~,max_inds2]=max(data2,[],1);
    max_inds=max_inds2+detect_interval*[0:(S2-1)];
    max_vals=data(max_inds);
    max_vals(((max_inds(1:end-1))>=max_inds(2:end)-detect_interval) & (max_vals(1:end-1)<max_vals(2:end)))=-1;
    max_vals(1+find((max_inds(2:end)<=max_inds(1:end-1)+detect_interval) & (max_vals(2:end)<=max_vals(1:end-1))))=-1;
    times=max_inds(max_vals>=detect_threshold); 
    if margin>0
        times=times((times>=margin)&(times<N-margin));
    end
    times = int64(times);
end