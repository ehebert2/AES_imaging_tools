%%%%%%%%%%%%%%%%%%%%
% object for reading in and reorganizing frames from file and rigid lateral 
% motion registration
%%%%%%%%%%%%%%%%%%%%

classdef MotionTracker < handle
    properties (Access = private)
        params
        channelsIn
        channelsOut
        mainChannel
        isVolume
        slicesIn
        slicesOut
        reslice
        reslicer

        tempBlur
        spaceBlur
        tempBlurWght
        spaceBlurWght

        aesMasks
        totMasks
        refImages
        fftRefImage
        corrX
        corrY
        corrBase
        adj_corr
        roiX
        roiY
        maxWidth
        maxHeight
        norm
        roiNumel
        refSub
        numRois
        maxNumRois

        bidirectional
        yBidirectional
        zBidirectional
        bidiShift        
        frameParity

        tin
        filename
        frames
        index
        resetLimit
        currentOffset
        buffer
        mtnBuffer
        smoothBuf
        updateMtnBfr
        fftXCorr
        intlWndw
        images
        firstFrame
        
        multiPass
        pass
        dispLog

        % stored to avoid reallocating memory
        useGPU = false;
        tempImages
        tempSl
        tempCh
        tempRois
        A
    end

    properties (GetAccess = public)
        h
        htot
    end

    properties (SetAccess = public)
        frameLimit
    end

    methods
        function obj = MotionTracker(params)
            obj.frameLimit = [];
            obj.updateMtnBfr = true;
            obj.params = validateParams(params);
            if (obj.params.splitChannels)
                obj.aesMasks = obj.params.roiAes(:,obj.params.mainChannel);
            else
                obj.aesMasks = obj.params.roiAes(:,1);
            end
            obj.mainChannel = obj.params.mainChannel;
            obj.bidirectional = obj.params.bidirectional;
            obj.yBidirectional = obj.params.yBidirectional;
            obj.zBidirectional = obj.params.zBidirectional;
            obj.bidiShift = obj.params.bidiShift;
            obj.isVolume = obj.params.volume;
            obj.reslice = obj.params.reslice;
            if (obj.reslice)
                obj.slicesIn = obj.params.reslicer.slicesIn;
                obj.slicesOut = obj.params.reslicer.slicesOut;
                obj.channelsIn = obj.params.reslicer.channelsIn;
                obj.channelsOut = obj.params.reslicer.channelsOut;
                obj.reslicer = obj.params.reslicer;
            else
                obj.slicesIn = obj.params.slices;
                obj.slicesOut = obj.params.slices;
                obj.channelsIn = obj.params.channels;
                obj.channelsOut = obj.params.channels;
            end
            obj.spaceBlur = obj.params.sptlWndw;
            obj.tempBlur = obj.params.mtnWndw;
            obj.fftXCorr = obj.params.fft;
            obj.intlWndw = obj.params.mtnIntlWndw;
            obj.multiPass = (obj.params.passes > 1);
            obj.filename = [];
            obj.pass = 0;
            obj.frames = 0;
            obj.totMasks = (zeros(size(obj.aesMasks{1},1),size(obj.aesMasks{1},2),obj.slicesOut)>0);
            obj.numRois = zeros(obj.slicesOut,1);
            for sl=1:obj.slicesOut
                obj.totMasks(:,:,sl) = (sum(obj.aesMasks{sl},3)>0);
                if (isempty(obj.aesMasks{sl}))
                    obj.numRois(sl) = 0;
                else
                    obj.numRois(sl) = size(obj.aesMasks{sl},3);
                end
            end

            % get rid of overlapping pixels in masks
            for sl=1:obj.slicesOut
                if (obj.numRois(sl) > 0)
                    keep = ones(1,1,obj.numRois(sl));
                    overlapMask = zeros(size(obj.totMasks,1),size(obj.totMasks,2));
                    for ii = 1:obj.numRois(sl)
                        obj.aesMasks{sl}(:,:,ii) = obj.aesMasks{sl}(:,:,ii) - obj.aesMasks{sl}(:,:,ii) .* overlapMask;
                        overlapMask = overlapMask + obj.aesMasks{sl}(:,:,ii);
                        if (sum(obj.aesMasks{sl}(:,:,ii),'all')==0)
                            keep(ii)=0;
                        end
                    end
                    obj.aesMasks{sl} = obj.aesMasks{sl}(:,:,(keep>0));
                    obj.numRois(sl) = sum(keep);
                end
            end
        end

        % initial setup before running through a full video
        function setup(obj,filename)
            obj.currentOffset = ceil((obj.bidirectional+1)*obj.tempBlur/2);

            % build weights for temporal smoothing of frames
            temp = linspace(-1,1,obj.tempBlur);
            temp = exp(-2*temp.^2);
            temp = temp / sum(temp);
            obj.tempBlurWght = reshape(temp,1,1,[]);

            if (isempty(obj.tin))
                obj.filename = filename;
                obj.tin = Tiff(filename,'r');
                obj.countFrames();
                obj.dispLog = zeros(obj.frames,2);
            end

            % setup for bidirectional data
            tempTiff = read(obj.tin);
            obj.htot = size(tempTiff,1);
            if (obj.isVolume || (~obj.yBidirectional))
                obj.h = obj.htot;
            else
                obj.h = floor(obj.htot/2);
            end

            % calculte variables used for accessing motion corrected sub
            % regions
            obj.maxNumRois = max(obj.numRois);
            obj.roiNumel = zeros(obj.slicesOut,obj.maxNumRois);
            for sl=1:obj.slicesOut
                for ii=1:obj.numRois(sl)
                    obj.roiNumel(sl,ii) = sum(obj.aesMasks{sl}(:,:,ii),'all');
                end
            end
            obj.images = int16(zeros(size(obj.totMasks,1),size(obj.totMasks,2),obj.slicesOut,obj.channelsOut));
            obj.tempImages = zeros(size(obj.totMasks,1),size(obj.totMasks,2),obj.slicesOut);
            obj.smoothBuf = zeros(size(obj.totMasks,1),size(obj.totMasks,2),obj.slicesOut);

            if (obj.fftXCorr)
                obj.corrBase = zeros((2*size(obj.totMasks,1)-1),(2*size(obj.totMasks,2)-1));
                obj.maxHeight = size(obj.totMasks,1);
                obj.maxWidth = size(obj.totMasks,2);

                if (obj.spaceBlur > 0)
                    tempX = obj.maxWidth*linspace(-0.5,0.5,obj.maxWidth);
                    tempX = repmat(tempX,obj.maxHeight,1);
                    tempY = obj.maxHeight*linspace(-0.5,0.5,obj.maxHeight);
                    tempY = repmat(tempY',1,obj.maxWidth);
                    blur = exp(-2*(tempX.^2+tempY.^2)/obj.spaceBlur^2);
                    blur = blur / sum(blur,'all');
                    obj.spaceBlurWght = obj.corrBase;
                    obj.spaceBlurWght(1:obj.maxHeight,1:obj.maxWidth) = blur;
                    obj.spaceBlurWght = fft2(obj.spaceBlurWght);
                end
            else
                if (obj.spaceBlur > 0)
                    tempWidth = ceil(obj.spaceBlur);
                    tempX = repmat((-tempWidth):tempWidth,(2*tempWidth+1),1);
                    tempY = tempX';
                    obj.spaceBlurWght = exp(-2*(tempX.^2+tempY.^2)/(obj.spaceBlur^2));
                    obj.spaceBlurWght = obj.spaceBlurWght / sum(obj.spaceBlurWght,'all');
                end
                obj.roiX = zeros(obj.slicesOut,obj.maxNumRois,2);
                obj.roiY = zeros(obj.slicesOut,obj.maxNumRois,2);
                for sl=1:obj.slicesOut
                    for ii=1:obj.numRois(sl)
                        proj = sum(squeeze(obj.aesMasks{sl}(:,:,ii)),1);
                        obj.roiX(sl,ii,1) = 1;
                        while(proj(obj.roiX(sl,ii,1))==0)
                            obj.roiX(sl,ii,1) = obj.roiX(sl,ii,1) + 1;
                        end
    
                        obj.roiX(sl,ii,2) = length(proj);
                        while(proj(obj.roiX(sl,ii,2))==0)
                            obj.roiX(sl,ii,2) = obj.roiX(sl,ii,2) - 1;
                        end
    
                        proj = sum(squeeze(obj.aesMasks{sl}(:,:,ii)),2);
                        obj.roiY(sl,ii,1) = 1;
                        while(proj(obj.roiY(sl,ii,1))==0)
                            obj.roiY(sl,ii,1) = obj.roiY(sl,ii,1) + 1;
                        end
    
                        obj.roiY(sl,ii,2) = length(proj);
                        while(proj(obj.roiY(sl,ii,2))==0)
                            obj.roiY(sl,ii,2) = obj.roiY(sl,ii,2) - 1;
                        end
                    end
                end
                roiWidth = obj.roiX(:,:,2) - obj.roiX(:,:,1) + 1;
                roiHeight = obj.roiY(:,:,2) - obj.roiY(:,:,1) + 1;
                obj.maxWidth = max(roiWidth,[],'all');
                obj.maxHeight = max(roiHeight,[],'all');
                obj.norm = zeros(2*obj.maxHeight-1,2*obj.maxWidth-1);
                obj.corrBase = zeros(2*obj.maxHeight-1,2*obj.maxWidth-1);
                obj.corrY = zeros(obj.slicesOut,obj.maxNumRois,2);
                obj.corrX = zeros(obj.slicesOut,obj.maxNumRois,2);
                obj.corrY(:,:,1) = obj.maxHeight - roiHeight + 1;
                obj.corrY(:,:,2) = obj.maxHeight + roiHeight - 1;
                obj.corrX(:,:,1) = obj.maxWidth - roiWidth + 1;
                obj.corrX(:,:,2) = obj.maxWidth + roiWidth - 1;
                for sl=1:obj.slicesOut
                    for ii=1:obj.numRois(sl)
                        temp = obj.corrBase;
                        temp(obj.corrY(sl,ii,1):obj.corrY(sl,ii,2),obj.corrX(sl,ii,1):obj.corrX(sl,ii,2)) = ones((2*roiHeight(sl,ii)-1),(2*roiWidth(sl,ii)-1));
                        obj.norm = obj.norm + temp;
                    end
                end
            end

            % initialize some more variables for buffer and ref image
            obj.fftRefImage = zeros(size(obj.corrBase,1),size(obj.corrBase,2),obj.slicesOut);
            obj.adj_corr = 0*obj.corrBase;
            obj.intlWndw = min(obj.intlWndw,floor(obj.frames/(obj.bidirectional+1)));
            obj.refSub = cell(obj.slicesOut,obj.maxNumRois);
            obj.resetLimit = round(4000*(obj.bidirectional+1)/(obj.slicesOut*obj.channelsOut));
            obj.frameParity = true;
            obj.buffer = cell(obj.slicesOut,obj.channelsOut);
            obj.mtnBuffer = cell(obj.slicesOut,1);
            obj.tempRois = cell(obj.slicesOut,obj.maxNumRois);
            obj.A = cell(obj.slicesOut,obj.maxNumRois);
            for sl=1:obj.slicesOut
                obj.mtnBuffer{sl} = zeros(size(obj.totMasks,1),size(obj.totMasks,2),obj.tempBlur);
                if (~obj.fftXCorr)
                    for ii=1:obj.numRois(sl)
                        obj.refSub{sl,ii} = zeros((obj.roiY(sl,ii,2)-obj.roiY(sl,ii,1)+1),(obj.roiX(sl,ii,2)-obj.roiX(sl,ii,1)+1));
                        obj.tempRois{sl,ii} = obj.refSub{sl,ii};
                        obj.A{sl,ii} = zeros((2*size(obj.refSub{sl,ii},1)-1),(2*size(obj.refSub{sl,ii},2)-1));
                    end
                end

                for ch = 1:obj.channelsOut
                    obj.buffer{sl,ch} = int16(zeros(size(obj.totMasks,1),size(obj.totMasks,2),(obj.tempBlur*(obj.bidirectional+1))));
                end
            end
            
            % considered using for gpu acceleration
            % not currently implemented because slow for small images
            % may be worth looking into for proper 3d correction
            if (obj.useGPU)
                obj.totMasks = gpuArray(obj.totMasks);
                obj.smoothBuf = gpuArray(obj.smoothBuf);
                obj.tempBlurWght = gpuArray(obj.tempBlurWght);
                obj.corrBase = gpuArray(obj.corrBase);
                obj.fftRefImage = gpuArray(obj.fftRefImage);
                obj.spaceBlurWght = gpuArray(obj.spaceBlurWght);
                obj.adj_corr = gpuArray(obj.adj_corr);
                for sl=1:obj.slicesOut
                    obj.mtnBuffer{sl} = gpuArray(obj.mtnBuffer{sl});
                    obj.aesMasks{sl} = gpuArray(obj.aesMasks{sl});
                    if (~obj.fftXCorr)
                        for ii=1:obj.numRois(sl)
                            obj.refSub{sl,ii} = gpuArray(obj.refSub{sl,ii});
                            obj.tempRois{sl,ii} = gpuArray(obj.tempRois{sl,ii});
                            obj.A{sl,ii} = gpuArray(obj.A{sl,ii});
                        end
                    end
                end
            end

            % build reference image and do additional motion correction
            % passes to clean up
            obj.buildRef();
            if (obj.params.intlPasses > 1)
                for p = 2:obj.params.intlPasses
                    obj.buildBuffer(1);
                    for fr = 1:(obj.intlWndw*(obj.bidirectional+1))
                        obj.readFrame();
                        obj.calculateDisp();
                    end
                    obj.buildRef();
                end
            end
        end

        % open new file or increment passes if doing multi-pass correction
        % complete setup before starting main motion correction loop
        function open(obj,filename)            
            if (~strcmp(obj.filename,filename))
                if (~isempty(obj.tin))
                    close(obj.tin);
                    obj.tin = [];
                end

                obj.tin = Tiff(filename,'r');
                obj.countFrames();
                obj.pass = 1;
                obj.dispLog = zeros(obj.frames,2);
                obj.filename = filename;
            else
                if (isempty(obj.tin))
                    obj.tin = Tiff(filename,'r');
                end
                obj.pass = obj.pass + 1;
            end

            obj.buildBuffer(1);
        end

        % read in one frame (all channels)
        function readFrame(obj)
            obj.index = obj.index + 1;
            obj.frameParity = (mod(obj.index,2) == 1);

            if (obj.index > (obj.frames - obj.currentOffset))
                if (obj.yBidirectional && ~obj.isVolume)
                    if (obj.frameParity)
                        for ch = 1:obj.channelsOut
                            obj.buffer{1,ch} = circshift(obj.buffer{1,ch},2,3);
                        end
                    end
                else
                    for sl=1:obj.slicesOut
                        for ch = 1:obj.channelsOut
                            obj.buffer{sl,ch} = circshift(obj.buffer{sl,ch},1,3);
                        end
                    end
                end

                if (obj.updateMtnBfr && ((~obj.bidirectional) || obj.frameParity))
                    for sl=1:obj.slicesOut
                        obj.mtnBuffer{sl} = circshift(obj.mtnBuffer{sl},1,3);
                        obj.mtnBuffer{sl}(:,:,1) = 0 * obj.mtnBuffer{sl}(:,:,1);
                    end
                end
            else
                if (mod(obj.index,obj.resetLimit) == 0)
                    tempDirIn = obj.tin.currentDirectory();
                    close(obj.tin);
                    clear obj.tin;
                    obj.tin = Tiff(obj.filename,'r');
                    obj.tin.setDirectory(tempDirIn);
                end
    
                if (obj.isVolume)
                    for sl=1:obj.slicesIn
                        for ch=1:obj.channelsIn
                            if (obj.firstFrame)
                                obj.firstFrame = false;
                            else
                                obj.tin.nextDirectory();
                            end

                            obj.tempCh = ch;
                            if (obj.zBidirectional && (~obj.frameParity))
                                obj.tempSl = obj.slicesIn-sl+1;
                            else
                                obj.tempSl = sl;
                            end

                            if (obj.reslice)
                                coord = obj.reslicer.getCoord(obj.tempSl,obj.tempCh);
                                obj.tempSl = coord(1);
                                obj.tempCh = coord(2);
                            end

                            obj.buffer{obj.tempSl,obj.tempCh} = circshift(obj.buffer{obj.tempSl,obj.tempCh},1,3);
                            if (obj.yBidirectional)
                                if (obj.frameParity)
                                    obj.buffer{obj.tempSl,obj.tempCh}(:,:,1) = circshift(read(obj.tin),obj.bidiShift,1);
                                else
                                    obj.buffer{obj.tempSl,obj.tempCh}(:,:,1) = flipud(read(obj.tin));
                                end
                            else
                                obj.buffer{obj.tempSl,obj.tempCh}(:,:,1) = read(obj.tin);
                            end
                        end
                    end
                else
                    if (obj.yBidirectional)
                        if (obj.frameParity)
                            for ch = 1:obj.channelsIn
                                if (obj.firstFrame)
                                    obj.firstFrame = false;
                                else
                                    obj.tin.nextDirectory();
                                end
                                obj.buffer{1,ch} = circshift(obj.buffer{1,ch},2,3);
                                temp = read(obj.tin);
                                obj.buffer{1,ch}(:,:,2) = flipud(temp((obj.htot-obj.h+1):obj.htot,:));
                                obj.buffer{1,ch}(:,:,1) = circshift(temp(1:obj.h,:),obj.bidiShift,1);
                            end
                        end
                    else
                        for ch = 1:obj.channelsOut
                            if (obj.firstFrame)
                                obj.firstFrame = false;
                            else
                                obj.tin.nextDirectory();
                            end
                            obj.buffer{1,ch} = circshift(obj.buffer{1,ch},1,3);
                            obj.buffer{1,ch}(:,:,1) = read(obj.tin);
                        end
                    end
                end

                if (obj.updateMtnBfr && ((~obj.bidirectional) || obj.frameParity))
                    for sl=1:obj.slicesOut
                        obj.mtnBuffer{sl} = circshift(obj.mtnBuffer{sl},1,3);
                        obj.mtnBuffer{sl}(:,:,1) = obj.processROIs(double(obj.buffer{sl,obj.mainChannel}(:,:,1)),sl,(obj.currentOffset+obj.index-1));
                    end
                end
            end
        end

        % main motion registration function
        function [xDisp,yDisp] = calculateDisp(obj)
            if ((~obj.bidirectional) || obj.frameParity)
                obj.adj_corr(:,:) = 0*obj.adj_corr;
                for sl=1:obj.slicesOut
                    obj.smoothBuf(:,:,sl) = sum(obj.mtnBuffer{sl} .* obj.tempBlurWght,3);
                    if (obj.fftXCorr)
                        obj.corrBase(1:obj.maxHeight,1:obj.maxWidth) = obj.smoothBuf(end:-1:1,end:-1:1,sl);
                        obj.adj_corr(:,:) = obj.adj_corr + real(ifft2(obj.fftRefImage(:,:,sl).*fft2(obj.corrBase)));
                    else
                        for jj = 1:obj.numRois(sl)
                            obj.tempRois{sl,jj}(:,:) = obj.smoothBuf(obj.roiY(sl,jj,1):obj.roiY(sl,jj,2),obj.roiX(sl,jj,1):obj.roiX(sl,jj,2),sl);
                            if (obj.spaceBlur > 0)
                                obj.tempRois{sl,jj}(:,:) = conv2(obj.tempRois{sl,jj},obj.spaceBlurWght,'same');
                            end
                            obj.A{sl,jj}(:,:) = xcorr2(obj.refSub{sl,jj},obj.tempRois{sl,jj});
                            obj.adj_corr(obj.corrY(sl,jj,1):obj.corrY(sl,jj,2),obj.corrX(sl,jj,1):obj.corrX(sl,jj,2)) = obj.adj_corr(obj.corrY(sl,jj,1):obj.corrY(sl,jj,2),obj.corrX(sl,jj,1):obj.corrX(sl,jj,2)) + obj.A{sl,jj};
                        end
                    end
                end
                [M1,temp_y] = max(obj.adj_corr);
                [~,temp_x] = max(M1);
                yDisp = temp_y(temp_x) - obj.maxHeight + obj.dispLog(obj.index,1);
                xDisp = temp_x - obj.maxWidth + obj.dispLog(obj.index,2);
                obj.dispLog(obj.index,:) = [yDisp,xDisp];
                if (obj.bidirectional)
                    obj.dispLog((obj.index+1),:) = [yDisp,xDisp];
                end
            end
        end

        function close(obj)
            close(obj.tin);
            obj.tin = [];
        end

        % Not currently usable. Needs to be rewritten
        function proj = buildProjection(obj,filename,frames,dlg)
            obj.tin = Tiff(filename,'r');
            if (frames == 0)
                obj.countFrames();
            else
                obj.frames = 1;
                while(~obj.tin.lastDirectory() && (obj.frames < obj.channelsOut*(frames+obj.params.mtnWndw)))
                    obj.frames = obj.frames + 1;
                end
                obj.frames = round(obj.frames/obj.channelsOut);
            end
            obj.dispLog = zeros((1+obj.bidirectional)*obj.frames,2);
            obj.params.mtnIntlWndw = min(obj.params.mtnIntlWndw,obj.frames-obj.params.mtnWndw);
            tempTiff = read(obj.tin);
            if (isempty(obj.aesMasks))
                obj.fftXCorr = true;
                obj.htot = size(tempTiff,1);
                if (obj.bidirectional)
                    obj.h = floor(obj.htot/2);
                else
                    obj.h = obj.htot;
                end
                obj.aesMasks = ones(obj.h,size(tempTiff,2),1);
                obj.totMasks = ones(obj.h,size(tempTiff,2));
            end
            obj.params.IntlPasses = 1;
            obj.setup(filename);
            obj.buildBuffer(1);
            proj = zeros(obj.h,size(tempTiff,2),obj.slicesOut,(obj.channelsOut*(obj.bidirectional+1)));
            notify_count = round(obj.frames/20);
            for fr = 1:obj.frames
                if (mod(fr,notify_count)==0)
                    dlg.Value = fr/obj.frames;
                end

                obj.readFrame();
                obj.calculateDisp();
                for ch=1:obj.channelsOut
                    if (obj.bidirectional)
                        proj(:,:,(2*ch-1)) = proj(:,:,(2*ch-1)) + double(obj.buffer{ch}(:,:,obj.currentOffset))/obj.frames;
                        proj(:,:,(2*ch)) = proj(:,:,(2*ch)) + double(obj.buffer{ch}(:,:,(obj.currentOffset+1)))/obj.frames;
                    else
                        proj(:,:,ch) = proj(:,:,ch) + double(obj.buffer{ch}(:,:,obj.currentOffset))/obj.frames;
                    end
                end
                obj.index = obj.index + obj.bidirectional;
            end
            close(obj.tin);
            obj.tin = [];
        end

        %% Getter and Setter Functions
        function image = getImage(obj,slice,channel)
            if (obj.yBidirectional && (~obj.isVolume))
                if (obj.frameParity)
                    image = obj.buffer{slice,channel}(:,:,obj.currentOffset);
                else
                    image = obj.buffer{slice,channel}(:,:,obj.currentOffset+1);
                end
            else
                image = obj.buffer{slice,channel}(:,:,obj.currentOffset);
            end
        end

        function images = getFrame(obj)
            if ((~obj.isVolume) && obj.yBidirectional && (~obj.frameParity))
                for ch=1:obj.channelsOut
                    obj.images(:,:,1,ch) = obj.buffer{ch}(:,:,obj.currentOffset+1);
                end
            else
                for sl=1:obj.slicesOut
                    for ch=1:obj.channelsOut
                        obj.images(:,:,sl,ch) = obj.buffer{sl,ch}(:,:,obj.currentOffset);
                    end
                end
            end
            images = obj.images;
        end

        function images = getReference(obj)
            if (obj.fftXCorr && (obj.spaceBlur > 0))
                for sl=1:obj.slicesOut
                    obj.corrBase(:,:) = 0*obj.corrBase;
                    obj.corrBase(1:obj.maxHeight,1:obj.maxWidth) = obj.refImages(:,:,sl);
                    obj.corrBase(:,:) = real(ifft2(fft2(obj.corrBase).*abs(obj.spaceBlurWght)));
                    obj.tempImages(:,:,sl) = obj.corrBase(1:obj.maxHeight,1:obj.maxWidth);
                    obj.corrBase(:,:) = 0*obj.corrBase;
                end
                images = obj.tempImages;
            else
                images = obj.refImages;
            end
        end

        function images = getCompImage(obj)
            if (obj.fftXCorr && (obj.spaceBlur > 0))
                for sl=1:obj.slicesOut
                    obj.corrBase(:,:) = 0*obj.corrBase;
                    obj.corrBase(1:obj.maxHeight,1:obj.maxWidth) = obj.smoothBuf(:,:,sl);
                    obj.corrBase(:,:) = real(ifft2(fft2(obj.corrBase).*abs(obj.spaceBlurWght)));
                    obj.tempImages(:,:,sl) = obj.corrBase(1:obj.maxHeight,1:obj.maxWidth);
                    obj.corrBase(:,:) = 0*obj.corrBase;
                end
                images = obj.tempImages;
            else
                images = obj.smoothBuf;
            end
        end

        function [xDisp,yDisp] = getDisp(obj)
            xDisp = obj.dispLog(obj.index,2);
            yDisp = obj.dispLog(obj.index,1);
        end

        function frames = getNumFrames(obj)
            frames = obj.frames;
        end

        function setIndex(obj,frame)
            obj.buildBuffer(frame);
        end
    end

    methods (Access = private)
        %% helper functions
        function countFrames(obj)
            obj.tin.setDirectory(1);
            obj.frames = 1;
            if (isempty(obj.frameLimit))
                while (~obj.tin.lastDirectory())
                    obj.frames = obj.frames + 1;
                    obj.tin.nextDirectory();
                end
            else
                while (~obj.tin.lastDirectory && (obj.frames < obj.frameLimit))
                    obj.frames = obj.frames+1;
                    obj.tin.nextDirectory();
                end
            end

            if (obj.isVolume)
                obj.frames = obj.frames / (obj.channelsOut*obj.slicesOut);
            else
                obj.frames = obj.frames * (obj.yBidirectional + 1) / obj.channelsOut;
            end
            obj.tin.setDirectory(1);
        end

        % prepare image for motion registration
        function image = processROIs(obj,image,slice,frame)
            image = image .* obj.totMasks(:,:,slice);
            for ii = 1:obj.numRois(slice)
                mV = sum((image.*obj.aesMasks{slice}(:,:,ii)),'all')/obj.roiNumel(slice,ii);
                image = image - mV * obj.aesMasks{slice}(:,:,ii);
            end
            image = circshift(image,obj.dispLog(frame,:));
        end

        % fill buffer ahead of current frame
        function buildBuffer(obj, frame)
            frame = max(frame,1);
            targetFrame = (ceil(frame/(1+obj.bidirectional))-1)*(1+obj.bidirectional);
            startIndex = targetFrame - (obj.params.mtnWndw-1)*(1+obj.bidirectional);
            numFramesToRead = (obj.params.mtnWndw-1)*(1+obj.bidirectional);
            startFrame = startIndex+obj.currentOffset;
            if (startFrame<1)
                numFramesToRead = numFramesToRead - (1-startFrame);
                startIndex = startIndex + (1-startFrame);
                startFrame = 1;
            end

            if (obj.yBidirectional && ~obj.isVolume)
                startDir = (startFrame-1)*obj.channelsOut/2;
            else
                startDir = (startFrame-1)*obj.channelsOut*obj.slicesOut;
            end

            obj.firstFrame = (startDir==0);
            startDir = startDir+obj.firstFrame;
            obj.tin.setDirectory(startDir);
            obj.index = startIndex;

            for sl=1:obj.slicesOut
                obj.mtnBuffer{sl} = 0*obj.mtnBuffer{sl};
                for ch = 1:obj.channelsOut
                    obj.buffer{sl,ch} = int16(0*obj.buffer{sl,ch});
                end
            end

            for ii=1:numFramesToRead
                obj.readFrame();
            end
        end

        % build initial reference image
        function buildRef(obj)
            obj.updateMtnBfr = false;
            obj.refImages = zeros(size(obj.totMasks,1),size(obj.totMasks,2),obj.slicesOut);
            obj.tin.setDirectory(1);
            obj.index = 0;
            obj.firstFrame = true;
            for ii = 1:(obj.intlWndw*(1+obj.bidirectional))
                obj.readFrame();
                if (obj.frameParity)
                    for sl=1:obj.slicesOut
                        obj.refImages(:,:,sl) = obj.refImages(:,:,sl) + obj.processROIs(double(obj.buffer{sl,obj.mainChannel}(:,:,1)),sl,obj.index);
                    end
                end
            end
            obj.refImages = obj.refImages / obj.intlWndw;

            obj.adj_corr = 0*obj.adj_corr;
            if (obj.fftXCorr)
                obj.fftRefImage = 0*obj.fftRefImage;
                for sl=1:obj.slicesOut
                    obj.fftRefImage(1:obj.maxHeight,1:obj.maxWidth,sl) = obj.refImages(:,:,sl);
                    obj.fftRefImage(:,:,sl) = fft2(obj.fftRefImage(:,:,sl));
                    if (obj.spaceBlur > 0)
                        obj.fftRefImage(:,:,sl) = obj.fftRefImage(:,:,sl) .* obj.spaceBlurWght .* conj(obj.spaceBlurWght);
                    end
                end
            else
                for sl=1:obj.slicesOut
                    if (obj.spaceBlur > 0)                
                        obj.refImages(:,:,sl) = conv2(obj.refImages(:,:,sl),obj.spaceBlurWght,'same');
                    end
        
                    for ii = 1:obj.numRois(sl)
                        obj.refSub{sl,ii}(:,:) = obj.refImages(obj.roiY(sl,ii,1):obj.roiY(sl,ii,2),obj.roiX(sl,ii,1):obj.roiX(sl,ii,2),sl);
                    end
                end
            end
            obj.updateMtnBfr = true;
        end
    end
end