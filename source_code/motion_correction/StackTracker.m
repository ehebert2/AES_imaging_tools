%%%%%%%%%%%%%%%%%%%
% Object used to read in frames, reorganize them, and read them out in
% sequence.
%%%%%%%%%%%%%%%%%%%

classdef StackTracker < handle
    properties
        channelsIn
        channelsOut
        bidirectional
        yBidirectional
        zBidirectional
        bidiShift
        isVolume
        slicesIn
        slicesOut

        filename
        tin
        frames
        index
        frameParity
        resetLimit
        buffer
        images
        firstFrame

        resliceMap
        reslice
        tempSl
        tempCh
    end

    properties (GetAccess = public)
        h
        htot
        w
    end

    properties (SetAccess = public)
        frameLimit
    end

    methods
        function obj = StackTracker(params)
            params = validateParams(params);
            obj.channelsIn = params.channels;
            obj.bidirectional = params.bidirectional;
            obj.yBidirectional = params.yBidirectional;
            obj.zBidirectional = params.zBidirectional;
            obj.bidiShift = params.bidiShift;
            obj.frameLimit = [];
            obj.slicesIn = params.slices;
            obj.isVolume = params.volume;
            obj.resetLimit = round(5000*(1+(1-obj.isVolume)*obj.bidirectional)/(obj.channelsIn*obj.slicesIn));
            obj.reslice = params.reslice;
            if (obj.reslice)
                obj.slicesOut = params.slicesOut;
                obj.channelsOut = params.channelsOut;
                obj.resliceMap = params.resliceMap;
            else
                obj.slicesOut = obj.slicesIn;
                obj.channelsOut = obj.channelsIn;
            end
        end

        % open tiff reader and setup image parameters
        function open(obj,filename)
            % assumes image dimensions are the same between files
            if (isempty(obj.filename))
                obj.tin = Tiff(filename,'r');
                tempTiff = read(obj.tin);
                obj.w = size(tempTiff,2);
                if (obj.isVolume)
                    obj.h = size(tempTiff,1);
                    obj.buffer = int16(zeros(obj.h,obj.w,obj.slicesOut,obj.channelsOut));
                    if (~obj.bidirectional)
                        obj.frameParity = true;
                    end
                elseif (obj.yBidirectional)
                    obj.htot = size(tempTiff,1);
                    obj.h = floor(obj.htot/2);
                    obj.buffer = int16(zeros(obj.h,size(tempTiff,2),1,(2*obj.channelsOut)));
                else
                    obj.buffer = int16(zeros(size(tempTiff,1),size(tempTiff,2),1,obj.channelsOut));
                    obj.h = size(tempTiff,1);
                end
                obj.frames = 0;
            elseif (~strcmp(obj.filename,filename))
                if (~isempty(obj.tin))
                    close(obj.tin);
                    obj.tin = [];
                end
                obj.tin = Tiff(filename,'r');
                obj.frames = 0;
            else
                if (isempty(obj.tin))
                    obj.tin = Tiff(filename,'r');
                else
                    obj.tin.setDirectory(1);
                end
            end
            obj.filename = filename;
            obj.index = 0;
            obj.firstFrame = true;
        end

        % read in frame and organize channels, reorient bidirectional data
        % if relevant
        function readFrame(obj)
            obj.index = obj.index+1;
            obj.frameParity = (mod(obj.index,2)==1);
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
                            temp = obj.resliceMap(obj.tempSl,obj.tempCh,1);
                            obj.tempCh = obj.resliceMap(obj.tempSl,obj.tempCh,2);
                            obj.tempSl = temp;
                        end

                        if (obj.yBidirectional)
                            if (obj.frameParity)
                                obj.buffer(:,:,obj.tempSl,obj.tempCh) = circshift(read(obj.tin),obj.bidiShift,1);
                            else
                                obj.buffer(:,:,obj.tempSl,obj.tempCh) = flipud(read(obj.tin));
                            end
                        else
                            obj.buffer(:,:,obj.tempSl,obj.tempCh) = read(obj.tin);
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
                            tempTiff = read(obj.tin);
                            obj.buffer(:,:,1,(2*ch-1)) = circshift(tempTiff(1:obj.h,:),obj.bidiShift,1);
                            obj.buffer(:,:,1,(2*ch)) = flipud(tempTiff((obj.htot-obj.h+1):obj.htot,:));
                        end
                    end
                else
                    for ch = 1:obj.channelsIn
                        if (obj.firstFrame)
                            obj.firstFrame = false;
                        else
                            obj.tin.nextDirectory();
                        end
                        obj.buffer(:,:,1,ch) = read(obj.tin);
                    end
                end
            end
        end

        % close reader
        function close(obj)
            close(obj.tin);
            obj.tin = [];
        end

        % build projection of first several frames (used for generating 
        % reference images)
        function proj = buildProjection(obj,filename,frames,dlg)
            obj.open(filename);
            proj = zeros(obj.h,obj.w,obj.slicesIn,(obj.channelsIn*(1+obj.bidirectional)));
            dlg.Indeterminate = "on";
            dlg.Title = "Counting Frames...";
            if (frames==0)
                obj.getNumFrames();
            else
                obj.frames = 1;
                tiffsPerFrame = obj.channelsIn*obj.slicesIn/(1+obj.yBidirectional*(~obj.isVolume));
                upperLim = frames*tiffsPerFrame;
                while(~obj.tin.lastDirectory && obj.frames < upperLim)
                    obj.frames = obj.frames + 1;
                    obj.tin.nextDirectory();
                end
                obj.frames = obj.frames/tiffsPerFrame;
            end
            dlg.Indeterminate = 'off';
            dlg.Title = 'Loading Image...';

            obj.tin.setDirectory(1);
            notify_count = round(obj.frames/20);
            obj.index = 0;
            obj.firstFrame = true;
            for fr=1:obj.frames
                if (mod(fr,notify_count)==0)
                    dlg.Message = ['Frame (',num2str(fr),'/',num2str(obj.frames),')'];
                    dlg.Value = fr/obj.frames;
                end

                obj.readFrame();
                if (obj.isVolume)
                    for sl=1:obj.slicesIn
                        for ch=1:obj.channelsIn
                            if (obj.bidirectional)
                                proj(:,:,sl,2*ch-obj.frameParity) = proj(:,:,sl,2*ch-obj.frameParity) + double(obj.buffer(:,:,sl,ch))*2/obj.frames;
                            else
                                proj(:,:,sl,ch) = proj(:,:,sl,ch) + double(obj.buffer(:,:,sl,ch))/obj.frames;
                            end
                        end
                    end
                else
                    if (obj.bidirectional)
                        if (obj.frameParity)
                            for ch=1:obj.channelsIn
                                proj(:,:,1,(2*ch-1))=proj(:,:,(2*ch-1))+double(obj.buffer(:,:,1,(2*ch-1)))*2/obj.frames;
                                proj(:,:,1,(2*ch))=proj(:,:,(2*ch))+double(obj.buffer(:,:,1,(2*ch)))*2/obj.frames;
                            end
                        end
                    else
                        for ch=1:obj.channelsIn
                            proj(:,:,1,ch) = proj(:,:,ch)+double(obj.buffer(:,:,1,ch))/obj.frames;
                        end
                    end
                end
            end
            close(obj.tin);
            obj.tin = [];
        end

        % get image stored from current frame index
        function image = getImage(obj,slice,channel)
            if ((~obj.isVolume) && obj.yBidirectional)
                if (obj.frameParity)
                    image = obj.buffer(:,:,1,(2*channel-1));
                else
                    image = obj.buffer(:,:,1,(2*channel));
                end
            else
                image = obj.buffer(:,:,slice,channel);
            end
        end

        % get all images stored from current frame index
        function images = getFrame(obj)
            if (obj.yBidirectional && (~obj.isVolume))
                if (obj.frameParity)
                    images = obj.buffer(:,:,1,(2*(1:obj.channelsIn)-1));
                else
                    images = obj.buffer(:,:,1,(2*(1:obj.channelsIn)));
                end
            else
                images = obj.buffer;
            end
        end

        % calculate number of frames if not known
        function frames = getNumFrames(obj)
            if (obj.frames < 1)
                obj.countFrames();
            end

            frames = obj.frames;
        end

        % set place in video
        function setIndex(obj, frame)
            if (obj.isVolume)
                realFrame = frame;
                obj.index = realFrame-1;
            else
                realFrame = ceil(frame/(1+obj.bidirectional));
                obj.index = (realFrame-1)*(1+obj.bidirectional);
            end

            obj.firstFrame = (obj.index <= 0);
            if (obj.firstFrame)
                obj.tin.setDirectory(1);           
            else
                obj.tin.setDirectory((realFrame-1)*obj.channelsIn*obj.slicesIn);
            end
        end
    end

    methods (Access = private)
        function countFrames(obj)
            obj.tin.setDirectory(1);
            obj.frames = 1;
            if (isempty(obj.frameLimit))
                while (~obj.tin.lastDirectory())
                    obj.frames = obj.frames + 1;
                    obj.tin.nextDirectory();
                end
            else
                while (~obj.tin.lastDirectory() && (obj.frames < obj.frameLimit))
                    obj.frames = obj.frames + 1;
                    obj.tin.nextDirectory();
                end
            end
            if (obj.isVolume)
                obj.frames = floor(obj.frames / (obj.channelsIn * obj.slicesIn));
            else
                obj.frames = floor(obj.frames * (obj.yBidirectional + 1) / obj.channelsIn);
            end
            obj.tin.setDirectory(1);
        end
    end
end