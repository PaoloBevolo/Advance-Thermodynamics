%script to add the mesh folders to the MATLAB path 
dirname=fileparts(mfilename('fullpath'));
if ~isempty(dirname)
    addpath(genpath(dirname));
    try
        Sh=regions.rect();
        clear ('Sh');
        disp (['Folder ' dirname ' and subfolders have been added to the MATLAB path']);
    catch exc
        error ('setpath:error','A problem occurred adding %s and subfolders to the MATLAB path',dirname );  
    end
end
clear ('dirname');