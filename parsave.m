function parsave(savepath, pop_frs, Vm)
    n_down = 10;
    Vm_down = Vm(1:n_down:end);
    keys = 1:size(Vm, 2);
    t_idx = 1:size(Vm, 1);
    recspikes = containers.Map;
    for key = keys
        recspikes(int2str(key)) = t_idx(Vm(:, key)==0);
    end
    save(savepath, 'pop_frs', 'Vm_down', 'recspikes')
end