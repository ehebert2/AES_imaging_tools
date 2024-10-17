%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Object for parsing output from motion correction app 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef AESOutputParser < handle
    properties (GetAccess = public)
        % flagged false if parser fails to complete
        status

        % traces
        aesTraces
        smplTraces
        bgTraces
        overlapTraces
        mtnTraces
        aesNumPx
        smplNumPx
        traceZeroed

        % roi names
        aesNames
        smplNames

        % masks
        aesMasks
        smplMasks
        maskHandler

        % channel info
        splitChannels
        channels

        % video and file info
        dimension
        basenames
        vidFiles
        projImages
        frames
    end

    methods
        function obj = AESOutputParser(path,interruptDialog)
            obj.status = false;
            useDlg = ~isempty(interruptDialog);

            temp = dir(path);
            temp = temp([temp.isdir]);
            obj.basenames = cell((length(temp)-3),1);
            numFiles = length(obj.basenames);
            ind = 0;
            for ii=3:(length(temp)-1)
                if (~strcmp(temp(ii).name,'masks'))
                    ind = ind+1;
                    obj.basenames{ind}=temp(ii).name;
                end
            end

            if (ind==length(obj.basenames))
                if (~strcmp(temp(end).name,'masks'))
                    disp('Mask folder not found.');
                    return;
                end
            else
                obj.basenames{numFiles}=temp(length(temp)).name;
            end
            
            maskPath = fullfile(path,'masks');

            if (~exist(strcat(path,filesep,obj.basenames{1},filesep,'bin_files'),'dir'))
                disp('No bin folder found.');
                return;
            end

            %% determine if split channel, number of channels, get names
            hasAES = false;
            hasSmpl = false;
            obj.channels = 1;
            obj.splitChannels = false;
            if (useDlg)
                interruptDialog.Title = "Parsing file structure...";
            end
            if (exist(fullfile(maskPath,'channel_1'),'dir'))
                obj.splitChannels = true;
                while (exist(fullfile(maskPath,strcat('channel_',num2str(obj.channels+1))),'dir'))
                    obj.channels = obj.channels + 1;
                end
                if (exist(strcat(maskPath,filesep,'channel_1',filesep,'AES'),'dir'))
                    hasAES = true;
                    obj.aesNames = cell(obj.channels,1);
                    for ch=1:obj.channels
                        obj.aesNames{ch} = obj.getFileList(strcat(maskPath,filesep,'channel_',num2str(ch),filesep,'AES'),'*.tif');
                    end
                end

                tempFile = strcat(path,filesep,obj.basenames{1},filesep,'bin_files',filesep,obj.basenames{1},'_roi_index.txt');
                if (exist(tempFile,'file'))
                    hasSmpl = true;
                    obj.smplNames = obj.readRoiIndex(tempFile,true,obj.channels);
                elseif (exist(strcat(maskPath,filesep,'channel_1',filesep,'sample'),'dir'))
                    hasSmpl = true;
                    obj.smplNames = cell(obj.channels,1);
                    for ch=1:obj.channels
                        obj.smplNames{ch} = obj.getFileList(strcat(maskPath,filesep,'channel_',num2str(ch),filesep,'sample'),'*.tif');
                    end
                end

                if (hasAES)
                    ch=1;
                    while(isempty(obj.aesNames{ch}))
                        ch=ch+1;
                    end
                    tempTiff = double(imread(strcat(maskPath,filesep,'channel_',num2str(ch),filesep,'AES',filesep,obj.aesNames{ch}{1},'.tif'),1));
                    obj.dimension = [size(tempTiff,1),size(tempTiff,2)];
                elseif (hasSmpl)
                    ch=1;
                    while(isempty(obj.smplNames{ch}))
                        ch=ch+1;
                    end
                    tempTiff = double(imread(strcat(maskPath,filesep,'channel_',num2str(ch),filesep,'sample',filesep,obj.smplNames{ch}{1},'.tif'),1));
                    obj.dimension = [size(tempTiff,1),size(tempTiff,2)];
                end
            else
                aesMaskPath = strcat(maskPath,filesep,'AES');
                if (exist(aesMaskPath,'dir'))
                    hasAES = true;
                    obj.aesNames = obj.getFileList(aesMaskPath,'*.tif');
                end

                tempFile = strcat(path,filesep,obj.basenames{1},filesep,'bin_files',filesep,obj.basenames{1},'_roi_index.txt');
                if (exist(tempFile,'file'))
                    hasSmpl = true;
                    obj.smplNames = obj.readRoiIndex(tempFile,false,0);
                elseif (exist(fullfile(maskPath,'sample'),'dir'))
                    hasSmpl = true;
                    obj.smplNames = obj.getFileList(fullfile(maskPath,'sample'),'*.tif');
                end

                tempFile = strcat(path,filesep,obj.basenames{1},filesep,'bin_files');
                obj.channels = obj.checkBinChannels(strcat(tempFile,filesep,obj.basenames{1},'_bg_mean'));
                if (obj.channels == 0)
                    if (exist(strcat(tempFile,filesep,'aes_roi_traces'),'dir'))
                        obj.channels = obj.checkBinChannels(strcat(tempFile,filesep,'aes_roi_traces',filesep,obj.basenames{1},'_',obj.aesNames{1}));
                        if (obj.channels==0)
                            obj.channels = obj.checkBinChannels(strcat(tempFile,filesep,'aes_roi_traces',filesep,obj.basenames{1},'_',obj.aesNames{1},'_zeroed'));
                        end
                    elseif (exist(strcat(tempFile,filesep,'smpl_roi_traces'),'dir'))
                        obj.channels = obj.checkBinChannels(strcat(tempFile,filesep,'smpl_roi_traces',filesep,obj.basenames{1},'_',obj.smplNames{1}));
                        if (obj.channels==0)
                            obj.channels = obj.checkBinChannels(strcat(tempFile,filesep,'smpl_roi_traces',filesep,obj.basenames{1},'_',obj.smplNames{1},'_zeroed'));
                        end
                    end
                    if (obj.channels==0)
                        disp('no meaningful traces found.');
                        return;
                    end
                end
                
                if (hasAES)
                    tempTiff = double(imread(strcat(maskPath,filesep,'AES',filesep,obj.aesNames{1},'.tif'),1));
                    obj.dimension = [size(tempTiff,1),size(tempTiff,2)];

                    temp = obj.aesNames;
                    obj.aesNames = cell(obj.channels,1);
                    for ch = 1:obj.channels
                        obj.aesNames{ch} = temp;
                    end
                end

                if (hasSmpl)
                    if (isempty(obj.dimension))
                        tempTiff = double(imread(strcat(maskPath,filesep,'sample',filesep,obj.smplNames{1},'.tif'),1));
                        obj.dimension = [size(tempTiff,1),size(tempTiff,2)];
                    end

                    temp = obj.smplNames;
                    obj.smplNames = cell(obj.channels,1);
                    for ch = 1:obj.channels
                        obj.smplNames{ch} = temp;
                    end
                end
            end

            %% Process Trace data
            % check which traces are present and initialize matrices
            if (useDlg)
                interruptDialog.Title = "Reading trace data...";
                interruptDialog.Cancelable = "on";
                interruptDialog.Indeterminate = "off";
                interruptDialog.Value = 0;
            end
            binPath = strcat(path,filesep,obj.basenames{1},filesep,'bin_files');
            hasAesTraces = false;
            hasSmplTraces = false;
            hasBgTraces = false;
            hasMtnTraces = false;
            hasOverlapTraces = false;
            obj.traceZeroed = false;
            if (exist(strcat(binPath,filesep,'aes_roi_traces'),'dir'))
                obj.aesTraces = cell(numFiles,1);
                hasAesTraces = true;
            end

            if (exist(strcat(binPath,filesep,'smpl_roi_traces'),'dir'))
                obj.smplTraces = cell(numFiles,1);
                hasSmplTraces = true;
            end

            if (exist(strcat(binPath,filesep,obj.basenames{1},'_displacement.bin'),'file'))
                obj.mtnTraces = cell(numFiles,1);
                hasMtnTraces = true;
            end

            if (obj.channels==1)
                if (exist(strcat(binPath,filesep,obj.basenames{1},'_bg_mean.bin'),'file'))
                    obj.bgTraces = cell(numFiles,1);
                    hasBgTraces = true;
                end
            else
                if (exist(strcat(binPath,filesep,obj.basenames{1},'_bg_mean_ch1.bin'),'file'))
                    obj.bgTraces = cell(numFiles,1);
                    hasBgTraces = true;
                end
            end

            if (obj.splitChannels)
                if (exist(strcat(binPath,filesep,obj.basenames{1},'_in_bounds_ch1.bin'),'file'))
                    obj.overlapTraces = cell(numFiles,1);
                    hasOverlapTraces = true;
                end
            else
                if (exist(strcat(binPath,filesep,obj.basenames{1},'_in_bounds.bin'),'file'))
                    obj.overlapTraces = cell(numFiles,1);
                    hasOverlapTraces = true;
                end
            end

            numTotTraces = hasMtnTraces+hasBgTraces*obj.channels+hasOverlapTraces*((obj.channels-1)*obj.splitChannels+1);
            if (hasAesTraces)
                for ch=1:obj.channels
                    numTotTraces = numTotTraces + length(obj.aesNames{ch});
                end
            end

            if (hasSmplTraces)
                for ch=1:obj.channels
                    numTotTraces = numTotTraces + length(obj.smplNames{ch});
                end
            end
            numTotTraces = numTotTraces * numFiles;
            numFinTraces = 0;

            if (hasSmplTraces)
                if (obj.channels==1)
                    if (exist(strcat(binPath,filesep,'smpl_roi_traces',filesep,obj.basenames{1},'_',obj.smplNames{ch}{1},'_zeroed.bin'),'file'))
                        obj.traceZeroed = true;
                    end
                else
                    ch=1;
                    while(isempty(obj.smplNames{ch}))
                        ch=ch+1;
                    end
                    if (exist(strcat(binPath,filesep,'smpl_roi_traces',filesep,obj.basenames{1},'_',obj.smplNames{ch}{1},'_zeroed_ch1.bin'),'file'))
                        obj.traceZeroed = true;
                    end
                end
            elseif (hasAesTraces)
                if (obj.channels==1)
                    if (exist(strcat(binPath,filesep,'aes_roi_traces',filesep,obj.basenames{1},'_',obj.aesNames{ch}{1},'_zeroed.bin'),'file'))
                        obj.traceZeroed = true;
                    end
                else
                    ch=1;
                    while(isempty(obj.smplNames{ch}))
                        ch=ch+1;
                    end
                    if (exist(strcat(binPath,filesep,'aes_roi_traces',filesep,obj.basenames{1},'_',obj.aesNames{ch}{1},'_zeroed_ch1.bin'),'file'))
                        obj.traceZeroed = true;
                    end
                end
            end

            % read in data
            for fl=1:numFiles
                if (useDlg)
                    interruptDialog.Message = strcat('File (',num2str(fl),'/',num2str(numFiles),')');
                end

                binPath = strcat(path,filesep,obj.basenames{fl},filesep,'bin_files');
                if (hasBgTraces)
                    if (obj.channels == 1)
                        obj.bgTraces{fl} = cell(1,1);
                        obj.bgTraces{fl}{1} = AESFile.readFullFile(strcat(binPath,filesep,obj.basenames{fl},'_bg_mean.bin'));
                        numFinTraces = numFinTraces+1;
                    else
                        tempName = strcat(binPath,filesep,obj.basenames{fl},'_bg_mean_ch');
                        obj.bgTraces{fl} = cell(obj.channels,1);
                        for ch=1:obj.channels
                            obj.bgTraces{fl}{ch} = AESFile.readFullFile(strcat(tempName,num2str(ch),'.bin'));
                            numFinTraces = numFinTraces+1;
                        end
                    end
                    if (useDlg)
                        interruptDialog.Value = numFinTraces/numTotTraces;
                        if (interruptDialog.CancelRequested)
                            return;
                        end
                    end
                end

                if (hasAesTraces)
                    obj.aesTraces{fl} = cell(obj.channels,1);
                    for ch = 1:obj.channels
                        numTraces = length(obj.aesNames{ch});
                        obj.aesTraces{fl}{ch} = cell(numTraces,1);
                        for ii=1:numTraces
                            fname = strcat(binPath,filesep,'aes_roi_traces',filesep,obj.basenames{fl},'_',obj.aesNames{ch}{ii});
                            if (obj.traceZeroed)
                                fname = strcat(fname,'_zeroed');
                            end

                            if (obj.channels==1)
                                fname = strcat(fname,'.bin');
                            else
                                fname = strcat(fname,'_ch',num2str(ch),'.bin');
                            end
                            temp = AESFile.readFullFile(fname);

                            if (hasBgTraces && obj.traceZeroed)
                                obj.aesTraces{fl}{ch}{ii} = [(mean(temp,2)-obj.bgTraces{fl}{ch}(:,1)),std(temp,0,2)];
                            else
                                obj.aesTraces{fl}{ch}{ii} = [mean(temp,2),std(temp,0,2)];
                            end
                            numFinTraces = numFinTraces+1;
                            if (useDlg)
                                interruptDialog.Value = numFinTraces/numTotTraces;
                                if (interruptDialog.CancelRequested)
                                    return;
                                end
                            end
                        end
                    end
                end
    
                if (hasSmplTraces)
                    obj.smplTraces{fl} = cell(obj.channels,1);
                    for ch = 1:obj.channels
                        numTraces = length(obj.smplNames{ch});
                        obj.smplTraces{fl}{ch} = cell(numTraces,1);
                        for ii=1:numTraces
                            fname = strcat(binPath,filesep,'smpl_roi_traces',filesep,obj.basenames{fl},'_',obj.smplNames{ch}{ii});
                            if (obj.traceZeroed)
                                fname = strcat(fname,'_zeroed');
                            end

                            if (obj.channels==1)
                                fname = strcat(fname,'.bin');
                            else
                                fname = strcat(fname,'_ch',num2str(ch),'.bin');
                            end
                            temp = AESFile.readFullFile(fname);

                            if (hasBgTraces&&obj.traceZeroed)
                                obj.smplTraces{fl}{ch}{ii} = [(mean(temp,2)-obj.bgTraces{fl}{ch}(:,1)),std(temp,0,2)];
                            else
                                obj.smplTraces{fl}{ch}{ii} = [mean(temp,2),std(temp,0,2)];
                            end
                            numFinTraces = numFinTraces+1;
                            if (useDlg)
                                interruptDialog.Value = numFinTraces/numTotTraces;
                                if (interruptDialog.CancelRequested)
                                    return;
                                end
                            end
                        end
                    end
                end
    
                if (hasMtnTraces)
                    obj.mtnTraces{fl} = AESFile.readFullFile(strcat(binPath,filesep,obj.basenames{fl},'_displacement.bin'));
                    numFinTraces = numFinTraces+1;
                    if (useDlg)
                        interruptDialog.Value = numFinTraces/numTotTraces;
                        if (interruptDialog.CancelRequested)
                            return;
                        end
                    end
                end
    
                if (hasOverlapTraces)
                    if (obj.splitChannels)
                        obj.overlapTraces{fl} = cell(obj.channels,1);
                        for ch=1:obj.channels
                            obj.overlapTraces{fl}{ch} = AESFile.readFullFile(fullfile(binPath,strcat(obj.basenames{fl},'_in_bounds_ch',num2str(ch),'.bin')));
                            numFinTraces = numFinTraces+1;
                        end
                    else
                        obj.overlapTraces{fl} = cell(obj.channels,1);
                        temp = AESFile.readFullFile(strcat(binPath,filesep,obj.basenames{fl},'_in_bounds.bin'));
                        for ch=1:obj.channels
                            obj.overlapTraces{fl}{ch} = temp;
                            numFinTraces = numFinTraces+1;
                        end
                    end
                    if (useDlg)
                        interruptDialog.Value = numFinTraces/numTotTraces;
                        if (interruptDialog.CancelRequested)
                            return;
                        end
                    end
                end
            end
            obj.traceZeroed = obj.traceZeroed && ~hasBgTraces;

            %% get video information
            % read in projection
            if (useDlg)
                interruptDialog.Message = '';
                interruptDialog.Title = "Reading in masks...";
                interruptDialog.Cancelable = "off";
                interruptDialog.Indeterminate = "on";
            end
            tempName = strcat(path,filesep,obj.basenames{1},filesep,obj.basenames{1});
            projAppend = [];
            if (exist(strcat(tempName,'_motion_corrected_projection.tif'),'file'))
                projAppend = '_motion_corrected_projection.tif';
            elseif (exist(strcat(tempName,'_projection.tif'),'file'))
                projAppend = '_projection.tif';
            end

            vidAppend = [];
            if (exist(strcat(tempName,'_motion_corrected.tif'),'file'))
                vidAppend = '_motion_corrected.tif';
            elseif (exist(strcat(tempName,'_compressed.tif'),'file'))
                vidAppend = '_compressed.tif';
            elseif (exist(strcat(tempName,'_processed.tif'),'file'))
                vidAppend = '_processed.tif';
            end

            if (~isempty(vidAppend))
                obj.vidFiles = cell(numFiles,1);
                for fl=1:numFiles                        
                    obj.vidFiles{fl} = strcat(path,filesep,obj.basenames{fl},filesep,obj.basenames{fl},vidAppend);
                end
            end

            if (isempty(projAppend))
                if (~isempty(vidAppend))
                    if (isempty(obj.dimension))
                        tempTiff = double(imread(strcat(tempName,vidAppend),1));
                        obj.dimension = [size(tempTiff,1),size(tempTiff,2)];
                    end

                    obj.projImages = zeros(obj.dimension(1),obj.dimension(2),obj.channels,length(obj.basenames));
                    params.bidirectional = false;
                    params.channels = obj.channels;
                    stackManager = StackTracker(params);
                    for fl=1:numFiles
                        obj.projImages(:,:,:,fl) = stackManager.buildProjection(obj.vidFiles{fl},3000,interruptDialog);
                    end
                else
                    if (isempty(obj.dimension))
                        disp('Not enough information.');
                        return;
                    end
                end
            else
                if (isempty(obj.dimension))
                    tempTiff = double(imread(strcat(tempName,projAppend),1));
                    obj.dimension = [size(tempTiff,1),size(tempTiff,2)];
                end
                obj.projImages = zeros(obj.dimension(1),obj.dimension(2),obj.channels,length(obj.basenames));
                for fl=1:length(obj.basenames)
                    tempFile = strcat(path,filesep,obj.basenames{fl},filesep,obj.basenames{fl},projAppend);
                    for ch=1:obj.channels
                        tempTiff = double(imread(tempFile,ch));
                        minVal = min(tempTiff,[],'all');
                        obj.projImages(:,:,ch,fl) = (tempTiff-minVal) / (max(tempTiff,[],'all')-minVal);
                    end
                end
            end

            %% read in masks
            obj.maskHandler = MaskHandler(obj.splitChannels,obj.channels,1,obj.dimension);
            if (obj.splitChannels)
                if (hasAES)
                    for ch=1:obj.channels
                        obj.maskHandler.setChannel(ch);
                        for ii=1:length(obj.aesNames{ch})
                            tempTiff = double(imread(strcat(maskPath,filesep,'channel_',num2str(ch),filesep,'AES',filesep,obj.aesNames{ch}{ii},'.tif'),1));
                            obj.maskHandler.addRoiAes(obj.aesNames{ch}{ii},(tempTiff>0));
                        end
                    end
                end

                if (hasSmpl)
                    for ch=1:obj.channels
                        obj.maskHandler.setChannel(ch);
                        for ii=1:length(obj.smplNames{ch})
                            tempTiff = double(imread(strcat(maskPath,filesep,'channel_',num2str(ch),filesep,'sample',filesep,obj.smplNames{ch}{ii},'.tif'),1));
                            obj.maskHandler.addRoiSmpl(obj.smplNames{ch}{ii},(tempTiff>0));
                        end
                    end
                end

                for ch=1:obj.channels
                    tempName = strcat(maskPath,filesep,'channel_',num2str(ch),filesep,'exp_mask.tif');
                    obj.maskHandler.setChannel(ch);
                    if (exist(tempName,'file'))
                        tempTiff = double(imread(tempName,1));
                        obj.maskHandler.setExpMask(tempTiff>0);
                    else
                        temp = obj.maskHandler.getRoiAes();
                        obj.maskHandler.setExpMask((sum(temp,3)>0));
                    end
                end
            else
                if (hasAES)
                    for ii=1:length(obj.aesNames{1})
                        tempTiff = double(imread(strcat(maskPath,filesep,'AES',filesep,obj.aesNames{1}{ii},'.tif'),1));
                        obj.maskHandler.addRoiAes(obj.aesNames{1}{ii},(tempTiff>0));
                    end
                end

                if (hasSmpl)
                    for ii=1:length(obj.smplNames{1})
                        tempTiff = double(imread(strcat(maskPath,filesep,'sample',filesep,obj.smplNames{1}{ii},'.tif'),1));
                        obj.maskHandler.addRoiSmpl(obj.smplNames{1}{ii},(tempTiff>0));
                    end
                end

                tempName = strcat(maskPath,filesep,'exp_mask.tif');
                if (exist(tempName,'file'))
                    tempTiff = double(imread(tempName,1));
                    obj.maskHandler.setExpMask((tempTiff>0));
                else
                    temp = obj.maskHandler.getRoiAes();
                    obj.maskHandler.setExpMask((sum(temp,3)>0));
                end
            end

            if (hasAesTraces)
                obj.aesNumPx = cell(obj.channels,1);
                for ch=1:obj.channels
                    numTraces = length(obj.aesNames{ch});
                    if (numTraces > 0)
                        obj.maskHandler.setChannel(ch);
                        tempMasks = obj.maskHandler.getRoiAes();
                        obj.aesNumPx{ch} = zeros(numTraces,1);
                        for ii=1:numTraces
                            obj.aesNumPx{ch}(ii) = sum(tempMasks(:,:,ii),'all');
                        end
                    end
                end
            end

            if (hasSmplTraces)
                obj.smplNumPx = cell(obj.channels,1);
                for ch=1:obj.channels
                    numTraces = length(obj.smplNames{ch});
                    if (numTraces > 0)
                        obj.maskHandler.setChannel(ch);
                        tempMasks = obj.maskHandler.getRoiSmpl();
                        obj.smplNumPx{ch} = zeros(numTraces,1);
                        for ii=1:numTraces
                            obj.smplNumPx{ch}(ii) = sum(tempMasks(:,:,ii),'all');
                        end
                    end
                end
            end

            %% get number of frames
            if (useDlg)
                interruptDialog.Title = "Finishing up...";
            end
            obj.frames = zeros(numFiles,1);
            if (~isempty(obj.bgTraces))
                for fl=1:numFiles
                    obj.frames(fl) = size(obj.bgTraces{fl}{1},1);
                end
            elseif (~isempty(obj.mtnTraces))
                for fl=1:numFiles
                    obj.frames(fl) = size(obj.mtnTraces{fl},1);
                end
            elseif (~isempty(obj.overlapTraces))
                for fl=1:numFiles
                    obj.frames(fl) = size(obj.overlapTraces{fl}{1},1);
                end
            elseif (hasAesTraces)
                ch=1;
                while(isempty(obj.aesTraces{1}{ch}))
                    ch=ch+1;
                end
                for fl=1:numFiles
                    obj.frames(fl) = size(obj.aesTraces{fl}{ch}{1},1);
                end
            elseif (hasSmplTraces)
                ch=1;
                while(isempty(obj.smplTraces{1}{ch}))
                    ch=ch+1;
                end
                for fl=1:numFiles
                    obj.frames(fl) = size(obj.smplTraces{fl}{ch}{1},1);
                end
            else
                disp('No meaningful trace data found.');
                return;
            end

            obj.status = true;
        end
    end

    methods (Static)
        % tries to find number of channels based on naming conventions of .bin files
        function ch = checkBinChannels(basename)
            ch = 0;
            if (exist(strcat(basename,'_ch1.bin'),'file'))
                ch = 1;
                while (exist(strcat(basename,'_ch',num2str(ch+1),'.bin'),'file'))
                    ch = ch + 1;
                end
            elseif (exist(strcat(basename,'.bin'),'file'))
                ch = 1;
            end
        end

        % get list of filenames in a folder with a given extension
        function names = getFileList(path,ext)
            temp = dir(fullfile(path,ext));
            names = cell(length(temp),1);
            for ii = 1:length(temp)
                temp2 = split(temp(ii).name,'.');
                names{ii} = temp2{1};
            end
        end

        % get ordered list of sample names from roi index file
        function names = readRoiIndex(fname,splitChannels,channels)
            fid = fopen(fname);
            tline = fgetl(fid);
            if (splitChannels)
                names = cell(channels,1);
                tline = fgetl(fid);
                while(ischar(tline))
                    name = split(tline,',');
                    ch = str2num(name{1});
                    names{ch} = [names{ch}, name(3)];
                    tline = fgetl(fid);
                end
            else
                names = {};
                tline = fgetl(fid);
                while (ischar(tline))
                    name = split(tline,',');
                    names = [names,name(2)];
                    tline = fgetl(fid);
                end
            end
        end
    end
end