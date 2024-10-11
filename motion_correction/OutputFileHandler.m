%%%%%%%%%%%%%%%%%%%
% handles file output when processing videos
%%%%%%%%%%%%%%%%%%%

classdef OutputFileHandler < handle
    properties
        params
        channels
        aesMasks
        smplMasks
        aesNames
        smplNames
        fullMask
        compMask
        numAes
        numSmpl
        overlapMask
        profileNorm

        tout
        firstFrame
        tagstruct
        fBgPrfl
        fAes
        fSmpl
        fBgm
        fMtn
        fOver
        projImages
        projFname
        numFrames

        image
        images
        xDisp
        yDisp

        getBG
        writeVid
    end

    methods
        function obj = OutputFileHandler(params)
            obj.params = validateParams(params);
            obj.channels = obj.params.channels;
            obj.aesMasks = obj.params.roiAes;
            obj.aesNames = obj.params.aesNames;
            obj.smplMasks = obj.params.roiSmpl;
            obj.smplNames = obj.params.smplNames;
            obj.fullMask = obj.params.expMask > 0;
            obj.numAes = obj.params.numRoiAes;
            obj.numSmpl = obj.params.numRoiSmpl;
            obj.compMask = int16(obj.fullMask);
        end

        % open video (false if something goes wrong)
        function status = open(obj,outPath,inFilename,numFrames)
            obj.numFrames = numFrames;
            basename = split(inFilename,'.');
            basename = basename{1};
            [status, msg] = mkdir(outPath,basename);
            if (~status)
                opts.WindowStyle = 'non-modal';
                waitfor(errordlg(msg,'Folder Write Issue',opts));
                return;
            end
            outPath = fullfile(outPath,strcat('\',basename));

            % store useful booleans so I don't have to recalc every frame
            multiChannel = (obj.channels > 1);
            obj.getBG = obj.params.bgTrace || obj.params.zeroVid;
            obj.writeVid = obj.params.mtnVid || obj.params.compress;
            obj.images = int16(zeros(obj.params.dim(1),obj.params.dim(2),obj.channels));

            hasTraces = obj.params.aesTrace + obj.params.smplTrace + obj.params.bgTrace + obj.params.mtnTrace;
            if (obj.params.splitChannels)
                hasTraces = hasTraces + (sum(obj.params.mtnOverlap) > 0);
            else
                hasTraces = hasTraces + obj.params.mtnOverlap;
            end
            hasTraces = hasTraces > 0;
            
            % determine output folder for trace data
            if (hasTraces)
                [status, msg] = mkdir(outPath,'bin_files');
                if (~status)
                    opts.WindowStyle = 'non-modal';
                    waitfor(errordlg(msg,'Folder Write Issue',opts));
                    return;
                end
                binPath = fullfile(outPath,'bin_files');
            end

            % build tagstruct for writing tiffs
            obj.tagstruct.ImageLength = obj.params.dim(1);
            obj.tagstruct.ImageWidth = obj.params.dim(2);
            obj.tagstruct.RowsPerStrip = obj.params.dim(1);
            obj.tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
            obj.tagstruct.BitsPerSample = 16;
            obj.tagstruct.SamplesPerPixel = 1;
            obj.tagstruct.SampleFormat = Tiff.SampleFormat.Int;
            obj.tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
            if (obj.params.compress)
                obj.tagstruct.Compression = Tiff.Compression.PackBits;
            else
                obj.tagstruct.Compression = Tiff.Compression.None;
            end
            
            % open output tiff stream
            if (obj.params.compress || obj.params.mtnVid)
                if (obj.params.mtnVid)
                    obj.tout = Tiff(fullfile(outPath,strcat(basename,'_motion_corrected.tif')),'w8');
                else
                    obj.tout = Tiff(fullfile(outPath,strcat(basename,'_compressed.tif')),'w8');
                end
                obj.firstFrame = true;
            end

            % open bg profile stream
            if (obj.params.bgProfile)
                obj.profileNorm = zeros(obj.params.dim(1),obj.channels);
                for ch=1:obj.channels
                    obj.profileNorm(:,ch) = sum((1-obj.fullMask(:,:,ch)),2);
                end

                if (prod(obj.profileNorm,'all') == 0)
                        disp('Warning: lack background samples across some rows of image. Turning off profile zeroing');
                        obj.params.bgProfile = false;
                        obj.params.bgTrace = true;
                else
                    obj.fBgPrfl = cell(1,obj.channels);
                    if (multiChannel)
                        for ii = 1:obj.channels
                            obj.fBgPrfl{ii} = AESFile.getWriter(fullfile(binPath,strcat(basename,'_bg_profile_ch',num2str(ii),'.bin')),obj.params.dim(1),numFrames,16,true);
                        end
                    else
                        obj.fBgPrfl{1} = AESFile.getWriter(fullfile(binPath,strcat(basename,'_bg_profile.bin')),obj.params.dim(1),numFrames,16,true);
                    end
                end
            end

            % open aes file streams
            if (obj.params.aesTrace)
                [status, msg] = mkdir(binPath,'aes_roi_traces');
                if (~status)
                    opts.WindowStyle = 'non-modal';
                    waitfor(errordlg(msg,'Folder Write Issue',opts));
                    return;
                end
                aesPath = fullfile(binPath,'aes_roi_traces');

                obj.fAes = cell(obj.channels);
                if (multiChannel)
                    for ch=1:obj.channels
                        obj.fAes{ch} = cell(obj.numAes(ch),1);
                        for ii=1:obj.numAes(ch)
                            obj.fAes{ch}{ii} = AESFile.getWriter(fullfile(aesPath,strcat(basename,'_',obj.aesNames{ch}{ii},'_ch',num2str(ch),'.bin')),sum(obj.aesMasks{ch}(:,:,ii),'all'),numFrames,16,true);
                        end
                    end
                else
                    obj.fAes{1} = cell(obj.numAes(1),1);
                    for ii=1:obj.numAes(1)
                        obj.fAes{1}{ii} = AESFile.getWriter(fullfile(aesPath,strcat(basename,'_',obj.aesNames{1}{ii},'.bin')),sum(obj.aesMasks{1}(:,:,ii),'all'),numFrames,16,true);
                    end
                end
            end

            % open sample file streams
            if (obj.params.smplTrace)
                [status, msg] = mkdir(binPath,'smpl_roi_traces');
                if (~status)
                    opts.WindowStyle = 'non-modal';
                    waitfor(errordlg(msg,'Folder Write Issue',opts));
                    return;
                end
                smplPath = fullfile(binPath,'smpl_roi_traces');

                obj.fSmpl = cell(obj.channels);
                if (multiChannel)
                    for ch=1:obj.channels
                        obj.fSmpl{ch} = cell(obj.numSmpl(ch),1);
                        for ii=1:obj.numSmpl(ch)
                            obj.fSmpl{ch}{ii} = AESFile.getWriter(fullfile(smplPath,strcat(basename,'_',obj.smplNames{ch}{ii},'_ch',num2str(ch),'.bin')),sum(obj.smplMasks{ch}(:,:,ii),'all'),numFrames,16,true);
                        end
                    end
                else
                    obj.fSmpl{1} = cell(obj.numSmpl(1),1);
                    for ii=1:obj.numSmpl(1)
                        obj.fSmpl{1}{ii} = AESFile.getWriter(fullfile(smplPath,strcat(basename,'_',obj.smplNames{1}{ii},'.bin')),sum(obj.smplMasks{1}(:,:,ii),'all'),numFrames,16,true);
                    end
                end
            end

            % open background mean/std filestreams
            if (obj.params.bgTrace)
                obj.fBgm = cell(1,obj.channels);
                if (multiChannel)
                    for ii = 1:obj.channels
                        obj.fBgm{ii} = AESFile.getWriter(fullfile(binPath,strcat(basename,'_bg_mean_ch',num2str(ii),'.bin')),2,numFrames,AESFile.DOUBLE,true);
                    end
                else
                    obj.fBgm{1} = AESFile.getWriter(fullfile(binPath,strcat(basename,'_bg_mean.bin')),2,numFrames,AESFile.DOUBLE,true);
                end
            end

            % generate overlap check matrices
            if (obj.params.splitChannels)
                if (sum(obj.params.mtnOverlap)>0)
                    obj.overlapMask = cell(obj.channels,1);
                    obj.fOver = cell(obj.channels,1);
                    numRoi = sum(obj.numSmpl(obj.params.mtnOverlap));
                    tCh = zeros(numRoi,1);
                    tInd = zeros(numRoi,1);
                    tName = cell(numRoi,1);
                    ind = 0;
                    for ch = 1:obj.channels
                        if (obj.params.mtnOverlap(ch))
                            obj.overlapMask{ch} = zeros((2*obj.params.dim(1)-1),(2*obj.params.dim(2)-1),obj.numSmpl(ch));
                            totalAes = squeeze(sum(obj.aesMasks{ch},3)) > 0;
                            totalAes = totalAes+0.0;
                            for ii=1:obj.numSmpl(ch)
                                obj.overlapMask{ch}(:,:,ii) = (xcorr2(totalAes,(obj.smplMasks{ch}(:,:,ii)+0.0))/sum(squeeze(obj.smplMasks{ch}(:,:,ii)),'all')) >= 1;
                                ind = ind+1;
                                tCh(ind) = ch;
                                tInd(ind) = ii;
                                tName{ind} = obj.smplNames{ch}{ii};
                            end
                            obj.fOver{ch} = AESFile.getWriter(fullfile(binPath,strcat(basename,'_in_bounds_ch',num2str(ch),'.bin')),obj.numSmpl(ch),numFrames,1,false);
                        end
                    end
                    T = table(tCh,tInd,categorical(tName),'VariableNames',{'Channel','Index','ROI Title'});
                    writetable(T,fullfile(binPath,strcat(basename,'_roi_index.txt')));
                end
            else
                if (obj.params.mtnOverlap)
                    obj.overlapMask = zeros((2*obj.params.dim(1)-1),(2*obj.params.dim(2)-1),obj.numSmpl(1));
                    totalAes = squeeze(sum(obj.aesMasks{1},3)) > 0;
                    totalAes = totalAes+0.0;
                    for ii=1:obj.numSmpl(1)
                        obj.overlapMask(:,:,ii) = (xcorr2(totalAes,(obj.smplMasks{1}(:,:,ii)+0.0))/sum(squeeze(obj.smplMasks{1}(:,:,ii)),'all')) >= 1;
                    end
                    T = table((1:obj.numSmpl(1))',categorical(obj.smplNames{1}'),'VariableNames',{'Index','ROI Title'});
                    writetable(T,fullfile(binPath,strcat(basename,'_roi_index.txt')));
                    obj.fOver = AESFile.getWriter(fullfile(binPath,strcat(basename,'_in_bounds.bin')),obj.numSmpl(1),numFrames,1,false);
                end
            end

            % open file streams if saving motion registration displacement
            if (obj.params.mtnTrace)
                obj.fMtn = AESFile.getWriter(fullfile(binPath,strcat(basename,'_displacement.bin')),2,numFrames,16,true);
            end

            if (obj.params.proj)
                obj.projImages = zeros(size(obj.images));
                if (obj.params.mtn)
                    obj.projFname = fullfile(outPath,strcat(basename,'_motion_corrected_projection.tif'));
                else
                    obj.projFname = fullfile(outPath,strcat(basename,'_projection.tif'));
                end
            end

            status = true;
            obj.xDisp = 0;
            obj.yDisp = 0;
        end

        % close file output streams
        function close(obj)
            if (obj.params.aesTrace)
                for ch = 1:obj.channels
                    for ii = 1:obj.numAes(ch)
                        fclose(obj.fAes{ch}{ii});
                    end
                end
            end

            if(obj.params.smplTrace)
                for ch = 1:obj.channels
                    for ii = 1:obj.numSmpl(ch)
                        fclose(obj.fSmpl{ch}{ii});
                    end
                end
            end

            if (obj.params.bgTrace)
                for ii = 1:obj.channels
                    fclose(obj.fBgm{ii});
                end
            end

            if (obj.params.bgProfile)
                for ii = 1:obj.channels
                    fclose(obj.fBgPrfl{ii});
                end
            end

            if (obj.params.mtn)
                if (obj.params.mtnTrace)
                    fclose(obj.fMtn);
                end

                if (obj.params.splitChannels)
                    for ch=1:obj.channels
                        if (obj.params.mtnOverlap(ch))
                            fclose(obj.fOver{ch});
                        end
                    end
                else
                    if (obj.params.mtnOverlap)
                        fclose(obj.fOver);
                    end
                end
            end

            if (obj.params.compress || obj.params.mtnVid)
                close(obj.tout);
            end

            if (obj.params.proj)
                toutProj = Tiff(obj.projFname,'w8');
                setTag(toutProj,obj.tagstruct);
                write(toutProj,int16(obj.projImages(:,:,1)));
                if (obj.channels > 1)
                    for ch=2:obj.channels
                        toutProj.writeDirectory();
                        setTag(toutProj,obj.tagstruct);
                        write(toutProj,int16(obj.projImages(:,:,ch)));
                    end
                end
                close(toutProj);
            end
        end

        % set current image to process
        function setImage(obj,image,channel)
            obj.images(:,:,channel) = image;
        end

        % set current images to process
        function setFrame(obj,frame)
            obj.images = frame;
        end

        % set displacement for motion tracking
        function setDisp(obj,xDisp,yDisp)
            obj.xDisp = xDisp;
            obj.yDisp = yDisp;
        end

        % write data to file
        function processFrame(obj)
            if (obj.params.mtnTrace)
                fwrite(obj.fMtn,[obj.xDisp,obj.yDisp],'int16','b');
            end

            if (obj.params.splitChannels)
                for ch=1:obj.channels
                    if (obj.params.mtnOverlap(ch))
                            fwrite(obj.fOver{ch},obj.overlapMask{ch}((obj.params.dim(1)-obj.yDisp),(obj.params.dim(2)-obj.xDisp),:),'ubit1','b');
                    end
                end
            else
                if (obj.params.mtnOverlap)
                    fwrite(obj.fOver,obj.overlapMask((obj.params.dim(1)-obj.yDisp),(obj.params.dim(2)-obj.xDisp),:),'ubit1','b');
                end
            end

            for ch = 1:obj.channels
                obj.image = obj.images(:,:,ch);

                if (obj.getBG)
                    bg = double(obj.image(~obj.fullMask(:,:,ch)));
                    bgMean = mean(bg);
                    if (obj.params.bgTrace)
                        fwrite(obj.fBgm{ch},[bgMean,std(bg)],'double','b');
                    end
                    bgMean = int16(bgMean);
                end
    
                if (obj.params.bgProfile)
                    profile = int16(sum((obj.image .* (1-obj.compMask(:,:,ch))),2) ./ obj.profileNorm(:,ch));
                    obj.image = obj.image - profile;
                    fwrite(obj.fBgPrfl{ch},profile,'int16','b');
                end
    
                if (obj.params.aesTrace)
                    for ii = 1:obj.numAes(ch)
                        fwrite(obj.fAes{ch}{ii},obj.image(obj.aesMasks{ch}(:,:,ii)),'int16','b');
                    end
                end
    
                if (obj.params.compress)
                    obj.image = obj.image .* obj.compMask(:,:,ch);
                    if (obj.params.zeroVid && ~obj.params.bgProfile)
                        obj.image = obj.image + (1 - obj.compMask(:,:,ch)) * bgMean;
                    end
                end
    
                if (obj.params.mtn)
                    obj.image = circshift(obj.image,[obj.yDisp,obj.xDisp]);
                end
    
                if (obj.params.smplTrace) 
                    for ii = 1:obj.numSmpl(ch)
                        fwrite(obj.fSmpl{ch}{ii},obj.image(obj.smplMasks{ch}(:,:,ii)),'int16','b');
                    end
                end
    
                if (obj.writeVid)
                    if (obj.params.zeroVid)
                        if (~obj.params.bgProfile)
                            obj.image = obj.image - bgMean;
                        end
                    elseif (obj.params.bgProfile)
                        obj.image = obj.image + circshift(profile,1,obj.yDisp);
                    end
        
                    if (obj.firstFrame)
                        obj.firstFrame = false;
                    else
                        obj.tout.writeDirectory();
                    end
                            
                    setTag(obj.tout,obj.tagstruct);
                    write(obj.tout,obj.image);
                end

                if (obj.params.proj)
                    obj.projImages(:,:,ch) = obj.projImages(:,:,ch) + double(obj.image)/obj.numFrames;
                end
            end
        end
    end
end