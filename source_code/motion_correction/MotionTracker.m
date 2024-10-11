%%%%%%%%%%%%%%%%%%%%
% object for reading in and reorganizing frames from file and rigid lateral 
% motion registration
%%%%%%%%%%%%%%%%%%%%

classdef MotionTracker < handle
    properties (Access = private)
        params
        channels
        mainChannel
        tempBlur
        spaceBlur
        tempBlurWght
        spaceBlurWght

        aesMask
        totMask
        refImage
        fftRefImage
        corrX
        corrY
        corrBase
        roiX
        roiY
        maxWidth
        maxHeight
        norm
        roiNumel
        roiSub
        refSub
        numRois

        bidirectional
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
        fftXCorr
        intlWndw
        images
        firstFrame
        
        multiPass
        pass
        dispLog
    end

    properties (GetAccess = public)
        h
        htot
    end

    methods
        function obj = MotionTracker(params)
            obj.params = validateParams(params);
            if (obj.params.splitChannels)
                obj.aesMask = obj.params.roiAes{obj.params.mainChannel};
            else
                obj.aesMask = obj.params.roiAes{1};
            end
            obj.channels = obj.params.channels;
            obj.mainChannel = obj.params.mainChannel;
            obj.bidirectional = obj.params.bidirectional;
            obj.bidiShift = obj.params.bidiShift;
            obj.spaceBlur = obj.params.sptlWndw;
            obj.tempBlur = obj.params.mtnWndw;
            obj.fftXCorr = obj.params.fft;
            obj.intlWndw = obj.params.mtnIntlWndw;
            obj.multiPass = (obj.params.passes > 1);
            obj.filename = [];
            obj.totMask = (sum(obj.aesMask,3)>0);
            obj.pass = 0;
            obj.frames = 0;
            if (isempty(obj.aesMask))
                obj.numRois = 0;
            else
                obj.numRois = size(obj.aesMask,3);
            end

            % get rid of overlapping pixels in masks
            if (obj.numRois > 0)
                keep = ones(1,1,obj.numRois);
                overlapMask = zeros(size(obj.totMask));
                for ii = 1:obj.numRois
                    obj.aesMask(:,:,ii) = obj.aesMask(:,:,ii) - obj.aesMask(:,:,ii) .* overlapMask;
                    overlapMask = overlapMask + obj.aesMask(:,:,ii);
                    if (sum(obj.aesMask(:,:,ii),'all')==0)
                        keep(ii)=0;
                    end
                end
                obj.aesMask = obj.aesMask(:,:,(keep>0));
                obj.numRois = sum(keep);
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
            if (obj.bidirectional)
                obj.h = floor(obj.htot/2);
            else
                obj.h = obj.htot;
            end

            % calculte variables used for accessing motion corrected sub
            % regions
            obj.roiNumel = zeros(obj.numRois,1);
            for ii=1:obj.numRois
                obj.roiNumel(ii) = sum(obj.aesMask(:,:,ii),'all');
            end
            obj.images = int16(zeros(size(obj.totMask,1),size(obj.totMask,2),obj.channels));

            if (obj.fftXCorr)
                obj.corrBase = zeros((2*size(obj.totMask)-1));
                obj.maxHeight = size(obj.totMask,1);
                obj.maxWidth = size(obj.totMask,2);

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
                    kernX = repmat(-obj.spaceBlur:obj.spaceBlur,(2*obj.spaceBlur+1),1)/obj.spaceBlur;
                    kernY = repmat((-obj.spaceBlur:obj.spaceBlur)',1,(2*obj.spaceBlur+1))/obj.spaceBlur;
                    obj.spaceBlurWght = exp(-2*(kernX.^2+kernY.^2));
                    obj.spaceBlurWght = obj.spaceBlurWght / sum(obj.spaceBlurWght,'all');
                end
                obj.roiX = zeros(obj.numRois,2);
                obj.roiY = zeros(obj.numRois,2);
                obj.roiSub = cell(obj.numRois,1);
                for ii=1:obj.numRois
                    proj = sum(squeeze(obj.aesMask(:,:,ii)),1);
                    obj.roiX(ii,1) = 1;
                    while(proj(obj.roiX(ii,1))==0)
                        obj.roiX(ii,1) = obj.roiX(ii,1) + 1;
                    end
    
                    obj.roiX(ii,2) = length(proj);
                    while(proj(obj.roiX(ii,2))==0)
                        obj.roiX(ii,2) = obj.roiX(ii,2) - 1;
                    end
    
                    proj = sum(squeeze(obj.aesMask(:,:,ii)),2);
                    obj.roiY(ii,1) = 1;
                    while(proj(obj.roiY(ii,1))==0)
                        obj.roiY(ii,1) = obj.roiY(ii,1) + 1;
                    end
    
                    obj.roiY(ii,2) = length(proj);
                    while(proj(obj.roiY(ii,2))==0)
                        obj.roiY(ii,2) = obj.roiY(ii,2) - 1;
                    end
    
                    obj.roiSub{ii} = squeeze(obj.aesMask(obj.roiY(ii,1):obj.roiY(ii,2),obj.roiX(ii,1):obj.roiX(ii,2),ii));
                end
    
                roiWidth = obj.roiX(:,2) - obj.roiX(:,1) + 1;
                roiHeight = obj.roiY(:,2) - obj.roiY(:,1) + 1;
                obj.maxWidth = max(roiWidth);
                obj.maxHeight = max(roiHeight);
                obj.norm = zeros(2*obj.maxHeight-1,2*obj.maxWidth-1);
                obj.corrBase = zeros(2*obj.maxHeight-1,2*obj.maxWidth-1);
                obj.corrY = zeros(length(roiHeight),2);
                obj.corrX = zeros(length(roiWidth),2);
                obj.corrY(:,1) = obj.maxHeight - roiHeight + 1;
                obj.corrY(:,2) = obj.maxHeight + roiHeight - 1;
                obj.corrX(:,1) = obj.maxWidth - roiWidth + 1;
                obj.corrX(:,2) = obj.maxWidth + roiWidth - 1;
                for ii = 1:obj.numRois
                    temp = obj.corrBase;
                    temp(obj.corrY(ii,1):obj.corrY(ii,2),obj.corrX(ii,1):obj.corrX(ii,2)) = ones((2*roiHeight(ii)-1),(2*roiWidth(ii)-1));
                    obj.norm = obj.norm + temp;
                end
            end

            % build reference image
            obj.intlWndw = min(obj.intlWndw,floor(obj.frames/(obj.bidirectional+1)));
            obj.refSub = cell(obj.numRois,1);
            obj.resetLimit = round(5000*(obj.bidirectional+1)/obj.channels);
            obj.frameParity = true;
            obj.buildRef();
            % do additional motion correction passes to create clearer
            % reference image
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
            if (obj.index > (obj.frames - obj.currentOffset))
                if (obj.bidirectional)
                    obj.frameParity = (mod(obj.index,2) == 1);
                    if (obj.frameParity)
                        for ch = 1:obj.channels
                            obj.buffer{ch} = circshift(obj.buffer{ch},2,3);
                        end
                    end
                else
                    for ch = 1:obj.channels
                        obj.buffer{ch} = circshift(obj.buffer{ch},1,3);
                    end
                end
                obj.mtnBuffer = circshift(obj.mtnBuffer,1,3);
                obj.mtnBuffer(:,:,1) = 0 * obj.mtnBuffer(:,:,1);
            else
                if (mod(obj.index,obj.resetLimit) == 0)
                    tempDirIn = obj.tin.currentDirectory();
                    close(obj.tin);
                    clear obj.tin;
                    obj.tin = Tiff(obj.filename,'r');
                    obj.tin.setDirectory(tempDirIn);
                end
    
                if (obj.bidirectional)
                    obj.frameParity = (mod(obj.index,2) == 1);
                    if(obj.frameParity)
                        for ch = 1:obj.channels
                            if (obj.firstFrame)
                                obj.firstFrame = false;
                            else
                                obj.tin.nextDirectory();
                            end
                            obj.buffer{ch} = circshift(obj.buffer{ch},2,3);
                            temp = read(obj.tin);
                            obj.buffer{ch}(:,:,2) = flipud(temp((obj.htot-obj.h+1):obj.htot,:));
                            obj.buffer{ch}(:,:,1) = circshift(temp(1:obj.h,:),obj.bidiShift,1);
                        end
                        obj.mtnBuffer = circshift(obj.mtnBuffer,1,3);
                        obj.mtnBuffer(:,:,1) = obj.processROIs(double(obj.buffer{obj.mainChannel}(:,:,1)),(obj.currentOffset+obj.index-1));
                    end
                else
                    for ch = 1:obj.channels
                        if (obj.firstFrame)
                            obj.firstFrame = false;
                        else
                            obj.tin.nextDirectory();
                        end
                        obj.buffer{ch} = circshift(obj.buffer{ch},1,3);
                        obj.buffer{ch}(:,:,1) = read(obj.tin);
                    end
                    obj.mtnBuffer = circshift(obj.mtnBuffer,1,3);
                    obj.mtnBuffer(:,:,1) = obj.processROIs(double(obj.buffer{obj.mainChannel}(:,:,1)),(obj.currentOffset+obj.index-1));
                end
            end
        end

        % main motion registration function
        function [xDisp,yDisp] = calculateDisp(obj)
            if (obj.frameParity)
                obj.smoothBuf = sum(obj.mtnBuffer .* obj.tempBlurWght,3);
                if (obj.fftXCorr)
                    obj.corrBase(1:obj.maxHeight,1:obj.maxWidth) = obj.smoothBuf(end:-1:1,end:-1:1);
                    adj_corr = real(ifft2(obj.fftRefImage.*fft2(obj.corrBase)));
                else
                    if (obj.spaceBlur > 0)
                        obj.smoothBuf = conv2(obj.smoothBuf,obj.spaceBlurWght,'same');
                    end
    
                    adj_corr = obj.corrBase;
                    for jj = 1:obj.numRois
                        tempRoi = obj.smoothBuf(obj.roiY(jj,1):obj.roiY(jj,2),obj.roiX(jj,1):obj.roiX(jj,2));
                        if (obj.spaceBlur > 0)
                            tempRoi = conv2(tempRoi,obj.spaceBlurWght,'same');
                        end
                        A = xcorr2(obj.refSub{jj},tempRoi);
                        adj_corr(obj.corrY(jj,1):obj.corrY(jj,2),obj.corrX(jj,1):obj.corrX(jj,2)) = adj_corr(obj.corrY(jj,1):obj.corrY(jj,2),obj.corrX(jj,1):obj.corrX(jj,2)) + A;
                    end
                end
                [M1,temp_y] = max(adj_corr);
                [M2,temp_x] = max(M1);
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

        function proj = buildProjection(obj,filename,frames,dlg)
            obj.tin = Tiff(filename,'r');
            if (frames == 0)
                obj.countFrames();
            else
                obj.frames = 1;
                while(~obj.tin.lastDirectory() && (obj.frames < obj.channels*(frames+obj.params.mtnWndw)))
                    obj.frames = obj.frames + 1;
                end
                obj.frames = round(obj.frames/obj.channels);
            end
            obj.dispLog = zeros((1+obj.bidirectional)*obj.frames,2);
            obj.params.mtnIntlWndw = min(obj.params.mtnIntlWndw,obj.frames-obj.params.mtnWndw);
            tempTiff = read(obj.tin);
            if (isempty(obj.aesMask))
                obj.fftXCorr = true;
                obj.htot = size(tempTiff,1);
                if (obj.bidirectional)
                    obj.h = floor(obj.htot/2);
                else
                    obj.h = obj.htot;
                end
                obj.aesMask = ones(obj.h,size(tempTiff,2),1);
                obj.totMask = ones(obj.h,size(tempTiff,2));
            end
            obj.params.IntlPasses = 1;
            obj.setup(filename);
            obj.buildBuffer(1);
            proj = zeros(obj.h,size(tempTiff,2),(obj.channels*(obj.bidirectional+1)));
            notify_count = round(obj.frames/20);
            for fr = 1:obj.frames
                if (mod(fr,notify_count)==0)
                    dlg.Value = fr/obj.frames;
                end

                obj.readFrame();
                obj.calculateDisp();
                for ch=1:obj.channels
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
        function image = getImage(obj,channel)
            if (obj.bidirectional)
                if (obj.frameParity)
                    image = obj.buffer{channel}(:,:,obj.currentOffset);
                else
                    image = obj.buffer{channel}(:,:,obj.currentOffset+1);
                end
            else
                image = obj.buffer{channel}(:,:,obj.currentOffset);
            end
        end

        function images = getFrame(obj)
            if (obj.bidirectional && ~obj.frameParity)
                for ch=1:obj.channels
                    obj.images(:,:,ch) = obj.buffer{ch}(:,:,obj.currentOffset+1);
                end
            else
                for ch=1:obj.channels
                    obj.images(:,:,ch) = obj.buffer{ch}(:,:,obj.currentOffset);
                end
            end
            images = obj.images;
        end

        function image = getReference(obj)
            image = obj.refImage;
        end

        function image = getCompImage(obj)
            image = obj.smoothBuf;
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
            while (~obj.tin.lastDirectory())
                obj.frames = obj.frames + 1;
                obj.tin.nextDirectory();
            end
            obj.frames = obj.frames * (obj.bidirectional + 1) / obj.channels;
            obj.tin.setDirectory(1);
        end

        % prepare image for motion registration
        function image = processROIs(obj,image,frame)
            image = image .* obj.totMask;
            for ii = 1:obj.numRois
                mV = sum((image.*obj.aesMask(:,:,ii)),'all')/obj.roiNumel(ii);
                image = image - mV * obj.aesMask(:,:,ii);
            end
            image = image .* obj.totMask;
            image = circshift(image,obj.dispLog(frame,:));
        end

        % fill buffer ahead of current frame
        function buildBuffer(obj, frame)
            frame = max(frame,1);
            realFrame = ceil(frame/(1+obj.bidirectional));
            obj.index = (realFrame-1)*(1+obj.bidirectional);
            startFrameDisp = min(ceil(obj.params.mtnWndw/2),realFrame)-1;
            obj.tin.setDirectory((realFrame-startFrameDisp-1)*obj.channels+1);

            obj.buffer = cell(obj.channels,1);
            obj.mtnBuffer = zeros(size(obj.refImage,1),size(obj.refImage,2),obj.tempBlur);
            for ii = 1:obj.channels
                obj.buffer{ii} = int16(zeros(size(obj.refImage,1),size(obj.refImage,2),(obj.tempBlur*(obj.bidirectional+1))));
            end

            obj.firstFrame = true;
            if (obj.bidirectional)
                for ii = (-startFrameDisp+1):(ceil(obj.tempBlur/2)-1)
                    for jj = 1:obj.channels
                        if (obj.firstFrame)
                            obj.firstFrame = false;
                        else
                            obj.tin.nextDirectory();
                        end
                        temp = read(obj.tin);
                        obj.buffer{jj}(:,:,obj.currentOffset-2*ii) = circshift(temp(1:obj.h,:),obj.bidiShift,1);
                        obj.buffer{jj}(:,:,obj.currentOffset-2*ii+1) = flipud(temp((obj.htot-obj.h+1):obj.htot,:));
                    end
                    obj.mtnBuffer(:,:,(ceil(obj.tempBlur/2)-ii)) = obj.processROIs(double(obj.buffer{obj.mainChannel}(:,:,obj.currentOffset-2*ii)),(obj.index+2*ii));
                end
            else
                for ii = (-startFrameDisp+1):(obj.currentOffset-1)
                    for jj = 1:obj.channels
                        if (obj.firstFrame)
                            obj.firstFrame = false;
                        else
                            obj.tin.nextDirectory();
                        end
                        obj.buffer{jj}(:,:,obj.currentOffset-ii) = read(obj.tin);
                    end
                    obj.mtnBuffer(:,:,obj.currentOffset-ii) = obj.processROIs(double(obj.buffer{obj.mainChannel}(:,:,obj.currentOffset-ii)),(obj.index+ii));
                end
            end
        end

        % build initial reference image
        function buildRef(obj)
            obj.refImage = zeros(size(obj.totMask,1),size(obj.totMask,2));
            for ii = 1:obj.intlWndw
                obj.tin.setDirectory((ii-1)*obj.channels+obj.mainChannel);
                tempTiff = double(read(obj.tin));
                if (obj.bidirectional)
                    tempTiff = circshift(tempTiff(1:obj.h,:),obj.bidiShift,1);
                end
                tempTiff = obj.processROIs(tempTiff,ii*(1+obj.bidirectional));
                obj.refImage = obj.refImage + tempTiff;
            end
            obj.refImage = obj.refImage / obj.intlWndw;

            if (obj.fftXCorr)
                obj.fftRefImage = obj.corrBase;
                obj.fftRefImage(1:obj.maxHeight,1:obj.maxWidth) = obj.refImage;
                obj.fftRefImage = fft2(obj.fftRefImage);
                if (obj.spaceBlur > 0)
                    obj.fftRefImage = obj.fftRefImage .* obj.spaceBlurWght.^2;
                end
            else
                if (obj.spaceBlur > 0)                
                    obj.refImage = conv2(obj.refImage,obj.spaceBlurWght,'same');
                end
    
                obj.refSub = cell(obj.numRois,1);
                for ii = 1:obj.numRois
                    obj.refSub{ii} = obj.refImage(obj.roiY(ii,1):obj.roiY(ii,2),obj.roiX(ii,1):obj.roiX(ii,2));
                end
            end
        end
    end
end