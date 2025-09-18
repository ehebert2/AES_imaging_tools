%%%%%%%%%%%%%%%%%%%
% handles file output when processing videos
%%%%%%%%%%%%%%%%%%%

classdef OutputFileHandler < handle
    properties
        params
        channels
        slices
        isVolume
        reslice
        aesMasks
        smplMasks
        aesNames
        smplNames
        bgMask
        fullMask16
        numAes
        numSmpl
        overlapMasks
        totAesMasks
        totSmplMasks
        totAesMasks16
        flMask
        imageDbl
        bgPx

        tout
        firstFrame
        tagstruct
        fBgPrfl
        fAes
        fSmpl
        fBgm
        fBgmFl
        fMtn
        fOver
        projImages
        projFname
        numFrames

        checkOverlap
        mtnOverlap
        checkAcc
        maxAcc2
        frame
        dxm1
        dxm2
        dym1
        dym2
        accPass

        image
        images
        xDisp
        yDisp

        getBG
        getFL
    end

    methods
        function obj = OutputFileHandler(params)
            obj.params = validateParams(params);
            obj.reslice = obj.params.reslice;
            if (obj.reslice)
                obj.channels = params.channelsOut;
                obj.slices = params.slicesOut;
            else
                obj.channels = obj.params.channels;
                obj.slices = params.slices;
            end
            obj.isVolume = params.volume && (obj.slices>1);
            obj.aesMasks = obj.params.roiAes;
            obj.aesNames = obj.params.aesNames;
            obj.smplMasks = obj.params.roiSmpl;
            obj.smplNames = obj.params.smplNames;
            obj.bgMask = obj.params.expMask <= 0;
            obj.bgPx = sum(obj.bgMask,'all');
            obj.numAes = obj.params.numRoiAes;
            obj.numSmpl = obj.params.numRoiSmpl;
            obj.fullMask16 = int16(~obj.bgMask);

            obj.getFL = obj.params.flTrace || obj.params.fillBG || obj.params.zeroTraceFl || obj.params.zeroVidFl;
            if (obj.getFL)
                obj.totAesMasks = zeros(size(obj.bgMask));
                obj.totSmplMasks = zeros(size(obj.bgMask));
                kernel = ones(3,3);
                for sl=1:obj.slices
                    for ch=1:obj.channels
                        obj.totAesMasks(:,:,sl,ch) = sum(obj.aesMasks{sl,ch},3);
                        if (~isempty(obj.smplMasks{sl,ch}))
                            obj.totSmplMasks(:,:,sl,ch) = conv2(sum(obj.smplMasks{sl,ch},3),kernel,'same');
                        end
                    end
                end
                obj.totAesMasks = obj.totAesMasks > 0;
                obj.totSmplMasks = obj.totSmplMasks > 0;
                obj.totAesMasks16 = int16(obj.totAesMasks);
            end
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
            multiSlice = (obj.slices > 1);
            obj.getBG = obj.params.bgTrace || obj.params.zeroVid || obj.params.zeroTrace;
            obj.images = int16(zeros(obj.params.dim(1),obj.params.dim(2),obj.slices,obj.channels));

            obj.checkOverlap = obj.params.mtnOverlap;
            obj.checkAcc = obj.params.checkAcc;
            obj.maxAcc2 = obj.params.maxAcc^2;
            if (obj.checkOverlap)
                obj.mtnOverlap = (zeros(obj.slices,obj.channels)>0);
                for sl=1:obj.slices
                    for ch=1:obj.channels
                        if ((obj.numAes(sl,ch)>0)&&(obj.numSmpl(sl,ch)>0))
                            obj.mtnOverlap(sl,ch)=true;
                        end
                    end
                end
                obj.checkOverlap = (sum(obj.mtnOverlap,'all')>0);
            end

            hasTraces = obj.params.aesTrace + obj.params.smplTrace + obj.params.bgTrace + obj.params.mtnTrace;
            hasTraces = hasTraces + obj.checkOverlap;
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
            if (obj.params.saveVid)
                vidName = fullfile(outPath,strcat(basename));
                changed = false;
                if (obj.params.mtn)
                    vidName = strcat(vidName,'_motion_corrected');
                    changed = true;
                elseif (obj.params.compress)
                    vidName = strcat(vidName,'_compressed');
                    changed = true;
                end

                if (obj.params.zeroVid || obj.params.zeroVidFl)
                    vidName = strcat(vidName,'_zeroed');
                    changed = true;
                end

                if (~changed)
                    vidName = strcat(vidName,'_processed');
                end
                vidName = strcat(vidName,'.tif');
                obj.tout = Tiff(vidName,'w8');
                obj.firstFrame = true;
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

                obj.fAes = cell(obj.slices,obj.channels);
                for sl=1:obj.slices
                    for ch=1:obj.channels
                        obj.fAes{sl,ch} = cell(obj.numAes(sl,ch),1);
                        for ii=1:obj.numAes(sl,ch)
                            aesName = strcat(basename,'_',obj.aesNames{sl,ch}{ii});
                            if (obj.params.zeroTrace)
                                aesName = strcat(aesName,'_zeroed');
                            end
                            
                            if (multiChannel)
                                aesName = strcat(aesName,'_ch',num2str(ch));
                            end
    
                            if (multiSlice)
                                aesName = strcat(aesName,'_sl',num2str(sl));
                            end
                            aesName = strcat(aesName,'.bin');
                            obj.fAes{sl,ch}{ii} = AESFile.getWriter(fullfile(aesPath,aesName),sum(obj.aesMasks{sl,ch}(:,:,ii),'all'),numFrames,16,true);
                        end
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

                obj.fSmpl = cell(obj.slices,obj.channels);
                for ch=1:obj.channels
                    for sl=1:obj.slices
                        obj.fSmpl{sl,ch} = cell(obj.numSmpl(sl,ch),1);
                        for ii=1:obj.numSmpl(sl,ch)
                            smplName = strcat(basename,'_',obj.smplNames{sl,ch}{ii});
                            if (obj.params.zeroTrace)
                                smplName = strcat(smplName,'_zeroed');
                            end

                            if (multiChannel)
                                smplName = strcat(smplName,'_ch',num2str(ch));
                            end

                            if (multiSlice)
                                smplName = strcat(smplName,'_sl',num2str(sl));
                            end
                            smplName = strcat(smplName,'.bin');
                            obj.fSmpl{sl,ch}{ii} = AESFile.getWriter(fullfile(smplPath,smplName),sum(obj.smplMasks{sl,ch}(:,:,ii),'all'),numFrames,16,true);
                        end
                    end
                end
            end

            % open background mean/std filestreams
            if (obj.params.bgTrace)
                obj.fBgm = cell(obj.slices,obj.channels);
                for sl=1:obj.slices
                    for ch=1:obj.channels
                        bgName = strcat(basename,'_bg_mean');
                        if (multiChannel)
                            bgName = strcat(bgName,'_ch',num2str(ch));
                        end

                        if (multiSlice)
                            bgName = strcat(bgName,'_sl',num2str(sl));
                        end
                        bgName = strcat(bgName,'.bin');
                        obj.fBgm{sl,ch} = AESFile.getWriter(fullfile(binPath,bgName),2,numFrames,AESFile.DOUBLE,true);
                    end
                end
            end

            if (obj.params.flTrace)
                obj.fBgmFl = cell(obj.slices,obj.channels);
                for sl=1:obj.slices
                    for ch=1:obj.channels
                        flName = strcat(basename,'_flbg_mean');
                        if (multiChannel)
                            flName = strcat(flName,'_ch',num2str(ch));
                        end

                        if (multiSlice)
                            flName = strcat(flName,'_sl',num2str(sl));
                        end
                        flName = strcat(flName,'.bin');
                        obj.fBgmFl{sl,ch} = AESFile.getWriter(fullfile(binPath,flName),2,numFrames,AESFile.DOUBLE,true);
                    end
                end
            end

            % generate overlap check matrices
            if (obj.checkOverlap)
                obj.overlapMasks = cell(obj.slices,obj.channels);
                obj.fOver = cell(obj.slices,obj.channels);
                numRoi = sum(obj.numSmpl.*obj.mtnOverlap,'all');
                tCh = zeros(numRoi,1);
                tSl = zeros(numRoi,1);
                tInd = zeros(numRoi,1);
                tName = cell(numRoi,1);
                ind=0;
                for sl=1:obj.slices
                    for ch=1:obj.channels
                        if (obj.mtnOverlap(sl,ch))
                            obj.overlapMasks{sl,ch} = zeros((2*obj.params.dim(1)-1),(2*obj.params.dim(2)-1),obj.numSmpl(sl,ch));
                            totalAes = squeeze(sum(obj.aesMasks{sl,ch},3))>0;
                            totalAes = totalAes+0.0;
                            for ii=1:obj.numSmpl(sl,ch)
                                obj.overlapMasks{sl,ch}(:,:,ii) = (xcorr2(totalAes,(obj.smplMasks{sl,ch}(:,:,ii)+0.0))/sum(squeeze(obj.smplMasks{sl,ch}(:,:,ii)),'all')) >= 1;
                                ind = ind+1;
                                tCh(ind) = ch;
                                tSl(ind) = sl;
                                tInd(ind) = ii;
                                tName{ind} = obj.smplNames{sl,ch}{ii};
                            end
                            overName = strcat(basename,'_in_bounds');
                            if (multiChannel)
                                overName = strcat(overName,'_ch',num2str(ch));
                            end

                            if (multiSlice)
                                overName = strcat(overName,'_sl',num2str(sl));
                            end
                            overName = strcat(overName,'.bin');
                            obj.fOver{sl,ch} = AESFile.getWriter(fullfile(binPath,overName),obj.numSmpl(sl,ch),numFrames,1,false);
                        end
                    end
                end

                if (multiSlice)
                    if (multiChannel)
                        T = table(tCh,tSl,tInd,categorical(tName),'VariableNames',{'Channel','Slice','Index','ROI Title'});
                    else
                        T = table(tSl,tInd,categorical(tName),'VariableNames',{'Slice','Index','ROI Title'});
                    end
                else
                    if (multiChannel)
                        T = table(tCh,tInd,categorical(tName),'VariableNames',{'Channel','Index','ROI Title'});
                    else
                        T = table(tInd,categorical(tName),'VariableNames',{'Index','ROI Title'});
                    end
                end
                writetable(T,fullfile(binPath,strcat(basename,'_roi_index.txt')));

                if (obj.checkAcc)
                    obj.frame = 1;
                    obj.dxm1 = 0;
                    obj.dxm2 = 0;
                    obj.dym1 = 0;
                    obj.dym2 = 0;
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
                for sl=1:obj.slices
                    for ch = 1:obj.channels
                        for ii = 1:obj.numAes(sl,ch)
                            fclose(obj.fAes{sl,ch}{ii});
                        end
                    end
                end
            end

            if(obj.params.smplTrace)
                for sl=1:obj.slices
                    for ch = 1:obj.channels
                        for ii = 1:obj.numSmpl(sl,ch)
                            fclose(obj.fSmpl{sl,ch}{ii});
                        end
                    end
                end
            end

            if (obj.params.bgTrace)
                for sl=1:obj.slices
                   for ch=1:obj.channels
                        fclose(obj.fBgm{sl,ch});
                    end
                end
            end

            if (obj.params.flTrace)
                for sl=1:obj.slices
                    for ch=1:obj.channels
                        fclose(obj.fBgmFl{sl,ch});
                    end
                end
            end

            if (obj.params.mtn)
                if (obj.params.mtnTrace)
                    fclose(obj.fMtn);
                end

                if (obj.checkOverlap)
                    for sl=1:obj.slices
                        for ch=1:obj.channels
                            if (obj.mtnOverlap(sl,ch))
                                fclose(obj.fOver{sl,ch});
                            end
                        end
                    end
                end
            end

            if (obj.params.saveVid)
                close(obj.tout);
            end

            if (obj.params.proj)
                toutProj = Tiff(obj.projFname,'w8');
                setTag(toutProj,obj.tagstruct);
                write(toutProj,int16(obj.projImages(:,:,1,1)));
                for sl=1:obj.slices
                    for ch=1:obj.channels
                        if ((sl>1)||(ch>1))
                            toutProj.writeDirectory();
                            setTag(toutProj,obj.tagstruct);
                            write(toutProj,int16(obj.projImages(:,:,sl,ch)));
                        end
                    end
                end
                close(toutProj);
            end
        end

        % set current image to process
        function setImage(obj,image,slice,channel)
            obj.images(:,:,slice,channel) = image;
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

            if (obj.checkOverlap)
                if (obj.checkAcc)
                    if (~obj.params.bidirectional || (obj.params.bidirectional&&(mod(obj.frame,2)==1)))
                        obj.accPass = ((obj.xDisp+obj.dxm2-2*obj.dxm1)^2+(obj.yDisp+obj.dym2-2*obj.dym1)^2)<=obj.maxAcc2;
                        obj.dxm2 = obj.dxm1;
                        obj.dym2 = obj.dym1;
                        obj.dxm1 = obj.xDisp;
                        obj.dym1 = obj.yDisp;
                    end

                    if (obj.frame==1)
                        obj.dxm1 = obj.xDisp;
                        obj.dym1 = obj.yDisp;
                    elseif (obj.frame==obj.numFrames)
                        for sl=1:obj.slices
                            for ch=1:obj.channels
                                if (obj.mtnOverlap(sl,ch))
                                    if (obj.accPass)
                                        fwrite(obj.fOver{sl,ch},obj.overlapMasks{sl,ch}((obj.params.dim(1)-obj.dym1),(obj.params.dim(2)-obj.dxm1),:),'ubit1','b');
                                        fwrite(obj.fOver{sl,ch},obj.overlapMasks{sl,ch}((obj.params.dim(1)-obj.yDisp),(obj.params.dim(2)-obj.xDisp),:),'ubit1','b');
                                    else
                                        fwrite(obj.fOver{sl,ch},zeros(obj.numSmpl(sl,ch))>0,'ubit1','b');
                                        fwrite(obj.fOver{sl,ch},zeros(obj.numSmpl(sl,ch))>0,'ubit1','b');
                                    end
                                end
                            end
                        end
                    else
                        for sl=1:obj.slices
                            for ch=1:obj.channels
                                if (obj.mtnOverlap(sl,ch))
                                    if (obj.accPass)
                                        fwrite(obj.fOver{sl,ch},obj.overlapMasks{sl,ch}((obj.params.dim(1)-obj.dym1),(obj.params.dim(2)-obj.dxm1),:),'ubit1','b');
                                    else
                                        fwrite(obj.fOver{sl,ch},zeros(obj.numSmpl(sl,ch))>0,'ubit1','b');
                                    end
                                end
                            end
                        end
                    end
                    obj.frame = obj.frame+1;
                else
                    for sl=1:obj.slices
                        for ch=1:obj.channels
                            if (obj.mtnOverlap(sl,ch))
                                fwrite(obj.fOver{sl,ch},obj.overlapMasks{sl,ch}((obj.params.dim(1)-obj.yDisp),(obj.params.dim(2)-obj.xDisp),:),'ubit1','b');
                            end
                        end
                    end
                end
            end

            for sl=1:obj.slices
                for ch = 1:obj.channels
                    obj.image = obj.images(:,:,sl,ch);
    
                    if (obj.getBG)
                        obj.imageDbl = double(obj.image);
                        bgMean = sum(obj.imageDbl.*obj.bgMask,'all')/obj.bgPx;
                        if (obj.params.bgTrace)
                            fwrite(obj.fBgm{sl,ch},[bgMean,sqrt(sum(((obj.imageDbl-bgMean).^2).*obj.bgMask,'all')/obj.bgPx)],'double','b');
                        end
                        bgMean = int16(bgMean);
                    end
        
                    if (obj.getFL)
                        obj.flMask = (obj.totAesMasks(:,:,sl,ch) - circshift(obj.totSmplMasks(:,:,sl,ch),[-obj.yDisp,-obj.xDisp]))>0;
                        flPx = sum(obj.flMask,'all');
                        if (flPx==0)
                            flMean = int16(0);
                            if (obj.params.flTrace)
                                fwrite(obj.fBgmFl{sl,ch},[0,0],'double','b');
                            end
                        else
                            if (~obj.getBG)
                                obj.imageDbl = double(obj.image);
                            end
                            flMean = sum(obj.imageDbl.*obj.flMask,'all')/flPx;
                            if (obj.params.flTrace)
                                fwrite(obj.fBgmFl{sl,ch},[flMean,sqrt(sum(((obj.imageDbl-flMean).^2).*obj.flMask,'all')/flPx)],'double','b');
                            end
                            flMean = int16(flMean);
                        end
                    end

                    if (obj.params.aesTrace)
                        if (obj.params.zeroTrace)
                            for ii = 1:obj.numAes(sl,ch)
                                fwrite(obj.fAes{sl,ch}{ii},(obj.image(obj.aesMasks{sl,ch}(:,:,ii))-bgMean),'int16','b');
                            end
                        elseif (obj.params.zeroTraceFl)
                            for ii=1:obj.numAes(sl,ch)
                                write(obj.fAes{sl,ch}{ii},(obj.image(obj.aesMasks{sl,ch}(:,:,ii))-flMean),'int16','b');
                            end
                        else
                            for ii = 1:obj.numAes(sl,ch)
                                fwrite(obj.fAes{sl,ch}{ii},obj.image(obj.aesMasks{sl,ch}(:,:,ii)),'int16','b');
                            end
                        end
                    end

                    if (obj.params.zeroVid)
                        obj.image = obj.image - bgMean;
                    elseif (obj.params.zeroVidFl)
                        obj.image = obj.image - flMean;
                    end
                    
                    if (obj.params.fillBG)
                        if (obj.params.zeroVidFl)
                            obj.image = obj.image .* obj.totAesMasks16(:,:,sl,ch);
                        elseif (obj.params.zeroVid)
                            obj.image = obj.image.*obj.totAesMasks16(:,:,sl,ch) + (flMean-bgMean) .* (1-obj.totAesMasks16(:,:,sl,ch));
                        else
                            obj.image = obj.image.*obj.totAesMasks16(:,:,sl,ch) + flMean .* (1-obj.totAesMasks16(:,:,sl,ch));
                        end
                    elseif (obj.params.compress)
                        obj.image = obj.image .* obj.fullMask16(:,:,sl,ch);
                    end
        
                    if (obj.params.mtn)
                        obj.image = circshift(obj.image,[obj.yDisp,obj.xDisp]);
                    end
        
                    if (obj.params.smplTrace)
                        traceDisp = int16(0);
                        if (obj.params.zeroVid)
                            traceDisp = -bgMean;
                        elseif (obj.params.zeroVidFl)
                            traceDisp = -flMean;
                        end

                        if (obj.params.zeroTrace)
                            traceDisp = traceDisp + bgMean;
                        elseif (obj.params.zeroTraceFl)
                            traceDisp = traceDisp + flMean;
                        end

                        for ii=1:obj.numSmpl(sl,ch)
                            fwrite(obj.fSmpl{sl,ch}{ii},(obj.image(obj.smplMasks{sl,ch}(:,:,ii))-traceDisp),'int16','b');
                        end
                    end
        
                    if (obj.params.saveVid)
                        if (obj.firstFrame)
                            obj.firstFrame = false;
                        else
                            obj.tout.writeDirectory();
                        end
                                
                        setTag(obj.tout,obj.tagstruct);
                        write(obj.tout,obj.image);
                    end
    
                    if (obj.params.proj)
                        obj.projImages(:,:,sl,ch) = obj.projImages(:,:,sl,ch) + double(obj.image)/obj.numFrames;
                    end
                end
            end
        end
    end
end