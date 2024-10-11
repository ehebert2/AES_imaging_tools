[fname,path] = uigetfile('.tif');
if (isequal(path,0))
    return;
end

ind = 66000;
tin = Tiff(fullfile(path,fname),"r");
%tin.setDirectory(ind);
for ii=1:ind
    if (mod(ii,10000)==0)
        disp(ii);
    end
    tin.nextDirectory();
end
close(tin);