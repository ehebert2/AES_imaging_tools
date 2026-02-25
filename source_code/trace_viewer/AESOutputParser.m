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
        flBgTraces
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

            fname = fullfile(path,'settings.json');
            if (~isfile(fname))
                return;
            end

            fid = fopen(fname);
            raw = fread(fid,inf);
            str = char(raw');
            fclose(fid);
            aesSettings = jsondecode(str);

            maskPath = fullfile(path,'masks');
            if (~isfolder(maskPath))
                return;
            end

            obj.basenames = aesSettings.Videos;
            numFiles = length(obj.basenames);

            if (~exist(strcat(path,filesep,obj.basenames{1},filesep,'bin_files'),'dir'))
                disp('No bin folder found.');
                return;
            end

            %% unpack json file
            hasAesTraces = aesSettings.Trace_Extraction.aes;
            hasSmplTraces = aesSettings.Trace_Extraction.sample;
            hasBgTraces = aesSettings.Trace_Extraction.electrical_bg;
            hasFlBgTraces = aesSettings.Trace_Extraction.fluorescent_bg;
            hasMtnTraces = aesSettings.Trace_Extraction.displacement;
            hasOverlapTraces = aesSettings.Trace_Extraction.in_bounds;
            temp = aesSettings.Trace_Extraction.background_subtracted;
            if (strcmp(temp,'electrical background'))
                obj.traceZeroed = Background.Electric;
            elseif (strcmp(temp,'fluorescent background'))
                obj.traceZeroed = Background.Fluorescent;
            else
                obj.traceZeroed = Background.None;
            end
            
            obj.channels = aesSettings.Stack_Properties.channels;
            obj.dimension = [aesSettings.Stack_Properties.dimensions.y,aesSettings.Stack_Properties.dimensions.x];
            
            % volume not currently implemented
            if (aesSettings.Stack_Properties.slices>1)
                return;
            end

            obj.aesNames = cell(obj.channels,1);
            obj.smplNames = cell(obj.channels,1);
            for ch = 1:obj.channels
                tempNames = aesSettings.ROIs.aes_names.(strcat('slice_1_channel_',num2str(ch)));
                obj.aesNames{ch} = split(tempNames,',');
                tempNames = aesSettings.ROIs.smpl_names.(strcat('slice_1_channel_',num2str(ch)));
                obj.smplNames{ch} = split(tempNames,',');
            end

            %% Process Trace data
            % check which traces are present and initialize matrices
            if (useDlg)
                interruptDialog.Title = "Reading trace data...";
                interruptDialog.Cancelable = "on";
                interruptDialog.Indeterminate = "off";
                interruptDialog.Value = 0;
            end
            
            numTotTraces = hasMtnTraces+(hasBgTraces+hasFlBgTraces+hasOverlapTraces)*obj.channels;
            if (hasAesTraces)
                obj.aesTraces = cell(numFiles,1);
                for ch=1:obj.channels
                    numTotTraces = numTotTraces + length(obj.aesNames{ch});
                end
            end

            if (hasSmplTraces)
                obj.smplTraces = cell(numFiles,1);
                for ch=1:obj.channels
                    numTotTraces = numTotTraces + length(obj.smplNames{ch});
                end
            end

            if (hasMtnTraces)
                obj.mtnTraces = cell(numFiles,1);
            end

            if (hasBgTraces)
                obj.bgTraces = cell(numFiles,1);
            end

            if (hasOverlapTraces)
                obj.overlapTraces = cell(numFiles,1);
            end
            numTotTraces = numTotTraces * numFiles;
            numFinTraces = 0;

            % read in data
            multiChannel = obj.channels>1;
            for fl=1:numFiles
                if (useDlg)
                    interruptDialog.Message = strcat('File (',num2str(fl),'/',num2str(numFiles),')');
                end

                binPath = strcat(path,filesep,obj.basenames{fl},filesep,'bin_files');
                if (hasBgTraces)
                    obj.bgTraces{fl} = cell(obj.channels,1);
                    tempName = strcat(binPath,filesep,obj.basenames{fl},'_bg_mean');
                    for ch=1:obj.channels
                        fname = tempName;
                        if (multiChannel)
                            fname = strcat(fname,'_ch',num2str(ch));
                        end
                        fname = strcat(fname,'.bin');
                        obj.bgTraces{fl}{ch} = AESFile.readFullFile(fname);
                        numFinTraces = numFinTraces+1;
                    end

                    if (useDlg)
                        interruptDialog.Value = numFinTraces/numTotTraces;
                        if (interruptDialog.CancelRequested)
                            return;
                        end
                    end
                end

                if (hasFlBgTraces)
                    obj.flBgTraces{fl} = cell(obj.channels,1);
                    tempName = strcat(binPath,filesep,obj.basenames{fl},'_flbg_mean');
                    for ch=1:obj.channels
                        fname = tempName;
                        if (multiChannel)
                            fname = strcat(fname,'_ch',num2str(ch));
                        end
                        fname = strcat(fname,'.bin');
                        obj.flBgTraces{fl}{ch} = AESFile.readFullFile(fname);
                        numFinTraces = numFinTraces+1;
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
                            if (obj.traceZeroed~=Background.None)
                                fname = strcat(fname,'_zeroed');
                            end

                            if (multiChannel)
                                fname = strcat(fname,'_ch',num2str(ch));
                            end
                            fname = strcat(fname,'.bin');
                            temp = AESFile.readFullFile(fname);
                            obj.aesTraces{fl}{ch}{ii} = [mean(temp,2),std(temp,0,2)];
                            if ((obj.traceZeroed==Background.Fluorescent)&&hasFlBgTraces)
                                obj.aesTraces{fl}{ch}{ii}(:,1) = obj.aesTraces{fl}{ch}{ii}(:,1) + obj.flBgTraces{fl}{ch}(:,1);
                            elseif ((obj.traceZeroed==Background.Electric)&&hasBgTraces)
                                obj.aesTraces{fl}{ch}{ii}(:,1) = obj.aesTraces{fl}{ch}{ii}(:,1) + obj.bgTraces{fl}{ch}(:,1);
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
                            if (obj.traceZeroed~=Background.None)
                                fname = strcat(fname,'_zeroed');
                            end

                            if (multiChannel)
                                fname = strcat(fname,'_ch',num2str(ch));
                            end
                            fname = strcat(fname,'.bin');
                            temp = AESFile.readFullFile(fname);
                            obj.smplTraces{fl}{ch}{ii} = [mean(temp,2),std(temp,0,2)];
                            if ((obj.traceZeroed==Background.Fluorescent)&&hasFlBgTraces)
                                obj.smplTraces{fl}{ch}{ii}(:,1) = obj.smplTraces{fl}{ch}{ii}(:,1) + obj.flBgTraces{fl}{ch}(:,1);
                            elseif ((obj.traceZeroed==Background.Electric)&&hasBgTraces)
                                obj.smplTraces{fl}{ch}{ii}(:,1) = obj.smplTraces{fl}{ch}{ii}(:,1) + obj.bgTraces{fl}{ch}(:,1);
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
                    obj.overlapTraces{fl} = cell(obj.channels,1);
                    for ch=1:obj.channels
                        fname = strcat(binPath,filesep,obj.basenames{fl},'_in_bounds');
                        if (multiChannel)
                            fname = strcat(fname,'_ch',num2str(ch));
                        end
                        fname = strcat(fname,'.bin');
                        obj.overlapTraces{fl}{ch} = AESFile.readFullFile(fname);
                        numFinTraces = numFinTraces+1;
                    end

                    if (useDlg)
                        interruptDialog.Value = numFinTraces/numTotTraces;
                        if (interruptDialog.CancelRequested)
                            return;
                        end
                    end
                end
            end
            if (((obj.traceZeroed==Background.Fluorescent)&&hasFlBgTraces)||((obj.traceZeroed==Background.Electric)&&hasBgTraces))
                obj.traceZeroed = Background.None;
            end

            %% get video information
            % read in projection
            if (useDlg)
                interruptDialog.Message = '';
                interruptDialog.Title = "Reading in masks...";
                interruptDialog.Cancelable = "off";
                interruptDialog.Indeterminate = "on";
            end

            if (aesSettings.Video_Properties.saved_video)
                changed = false;
                if (aesSettings.Motion_Correction.motion_corrected)
                    changed = true;
                    vidAppend = '_motion_corrected';
                elseif (aesSettings.Video_Properties.compressed)
                    changed = true;
                    vidAppend = '_compressed';
                end

                temp = aesSettings.Video_Properties.background_subtracted;
                if (~strcmp(temp,'none'))
                    changed = true;
                    vidAppend = strcat(vidAppend,'_zeroed');
                end

                if (~changed)
                    vidAppend = 'processed';
                end

                vidAppend = strcat(vidAppend,'.tif');
                obj.vidFiles = cell(numFiles,1);
                for fl=1:numFiles
                    obj.vidFiles{fl} = strcat(path,filesep,obj.basenames{fl},filesep,obj.basenames{fl},vidAppend);
                end
            end

            obj.projImages = zeros(obj.dimension(1),obj.dimension(2),obj.channels,length(obj.basenames));
            if (aesSettings.Video_Properties.saved_projection)
                if (aesSettings.Motion_Correction.motion_corrected)
                    projAppend = '_motion_corrected_projection.tif';
                else
                    projAppend = '_projection.tif';
                end

                for fl=1:length(obj.basenames)
                    tempFile = strcat(path,filesep,obj.basenames{fl},filesep,obj.basenames{fl},projAppend);
                    for ch=1:obj.channels
                        tempTiff = double(imread(tempFile,ch));
                        minVal = min(tempTiff,[],'all');
                        obj.projImages(:,:,ch,fl) = (tempTiff-minVal) / (max(tempTiff,[],'all')-minVal);
                    end
                end
            elseif (aesSettings.Video_Properties.saved_video)
                params.bidirectional = false;
                params.channels = obj.channels;
                stackManager = StackTracker(params);
                for fl=1:numFiles
                    obj.projImages(:,:,:,fl) = stackManager.buildProjection(obj.vidFiles{fl},3000,interruptDialog);
                    for ch=1:obj.channels
                        tempTiff = obj.projImages(:,:,ch,fl);
                        minVal = min(tempTiff,[],'all');
                        obj.projImages(:,:,ch,fl) = (tempTiff-minVal) / (max(tempTiff,[],'all')-minVal);
                    end
                end
            end

            %% read in masks
            obj.maskHandler = MaskHandler(obj.channels,1,obj.dimension);
            for ch=1:obj.channels
                obj.maskHandler.setChannel(ch);
                for ii=1:length(obj.aesNames{ch})
                    fname = strcat(maskPath,filesep,'AES');
                    if (multiChannel)
                        fname = strcat(fname,filesep,'channel_',num2str(ch));
                    end
                    fname = strcat(fname,filesep,obj.aesNames{ch}{ii},'.tif');
                    tempTiff = double(imread(fname,1));
                    obj.maskHandler.addRoiAes(obj.aesNames{ch}{ii},(tempTiff>0));
                end

                for ii=1:length(obj.smplNames{ch})
                    fname = strcat(maskPath,filesep,'sample');
                    if (multiChannel)
                        fname = strcat(fname,filesep,'channel_',num2str(ch));
                    end
                    fname = strcat(fname,filesep,obj.smplNames{ch}{ii},'.tif');
                    tempTiff = double(imread(fname,1));
                    obj.maskHandler.addRoiSmpl(obj.smplNames{ch}{ii},(tempTiff>0));
                end

                fname = strcat(maskPath,filesep,'exposure',filesep,'exp_mask');
                if (multiChannel)
                    fname = strcat(fname,'_ch',num2str(ch));
                end
                fname = strcat(fname,'.tif');
                tempTiff = double(imread(fname,1));
                obj.maskHandler.setExpMask((tempTiff>0));
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