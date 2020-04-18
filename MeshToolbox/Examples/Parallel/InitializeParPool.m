
if exist('parpool','file')
    pools=parpool('size'); 
    cpus=feature('numCores');
    parpool('open',cpus-1);
    if pools ~=cpu-1
        if pools>0
            parpool('close');
        end
        parpool('open', cpus-1);
    end
else
    warning([mfilename ':ParToolboxNotInstalled'],'The Parallel Computing Toolbox is NOT installed');
end