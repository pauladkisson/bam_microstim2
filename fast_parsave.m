function fast_parsave(savepath, Vm)
    keys = 1:size(Vm, 2);
    t_idx = 1:size(Vm, 1);
    recspikes = containers.Map;
    for key = keys
        recspikes(int2str(key)) = t_idx(Vm(:, key)==0);
    end
    save(savepath, 'recspikes')
end