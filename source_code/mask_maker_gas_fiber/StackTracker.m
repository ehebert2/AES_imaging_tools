%%%%%%%%%%%%%%%%%%%
% Object used to read in frames, reorganize them, and read them out in
% sequence.
%%%%%%%%%%%%%%%%%%%

classdef StackTracker < handle
    properties
        channels
        bidirectional
        bidiShift
        volume
        slices

        filename
        tin
        frames
        index
        frameParity
        resetLimit
        buffer
        images
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
            obj.channels = params.channels;
            obj.bidirectional = params.bidirectional;
            obj.bidiShift = params.bidiShift;
            obj.frameLimit = [];
            obj.slices = params.slices;
            obj.volume = params.volume;
            obj.resetLimit = round(5000*(1+(1-obj.volume)*obj.bidirectional)/(obj.channels*obj.slices));
        end

        % open tiff reader and setup image parameters
        function open(obj,filename)
            % assumes image dimensions are the same between files
            if (isempty(obj.filename))
                obj.tin = Tiff(filename,'r');
                tempTiff = read(obj.tin);
                obj.w = size(tempTiff,2);
                if (obj.volume)
                    obj.h = size(tempTiff,1);
                    obj.buffer = int16(zeros(obj.h,obj.w,obj.slices,obj.channels));
                    if (~obj.bidirectional)
                        obj.frameParity = true;
                    end
                elseif (obj.bidirectional)
                    obj.htot = size(tempTiff,1);
                    obj.h = floor(obj.htot/2);
                    obj.buffer = int16(zeros(obj.h,size(tempTiff,2),(2*obj.channels)));
                else
                    obj.buffer = int16(zeros(size(tempTiff,1),size(tempTiff,2),obj.channels));
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
        end

        % read in frame and organize channels, reorient bidirectional data
        % if relevant
        function readFrame(obj)
            obj.index = obj.index+1;
            if (mod(obj.index,obj.resetLimit) == 0)
                tempDirIn = obj.tin.currentDirectory();
                close(obj.tin);
                clear obj.tin;
                obj.tin = Tiff(obj.filename,'r');
                obj.tin.setDirectory(tempDirIn);
            end

            if (obj.volume)
                obj.frameParity = (mod(obj.index,2)==1);
                    
                if (obj.frameParity)
                    for sl = 1:obj.slices
                        for ch = 1:obj.channels
                            if (obj.index > 1 || ch > 1 || sl > 1)
                                obj.tin.nextDirectory();
                            end
                            obj.buffer(:,:,sl,ch) = read(obj.tin);
                        end
                    end
                else
                    for sl = 1:obj.slices
                        for ch = 1:obj.channels
                            obj.tin.nextDirectory();
                            if (obj.bidirectional)
                                obj.buffer(:,:,obj.slices-sl+1,ch) = flipud(read(obj.tin));
                            else
                                obj.buffer(:,:,obj.slices-sl+1,ch) = read(obj.tin);
                            end
                        end
                    end
                end
            else
                if (~obj.bidirectional)
                    for ch = 1:obj.channels
                        if (obj.index > 1 || ch > 1)
                            obj.tin.nextDirectory();
                        end
                        obj.buffer(:,:,ch) = read(obj.tin);
                    end
                else
                    obj.frameParity = (mod(obj.index,2)==1);
                    if (obj.frameParity)
                        for ch = 1:obj.channels
                            if (obj.index > 1 || ch > 1)
                                obj.tin.nextDirectory();
                            end
                            tempTiff = read(obj.tin);
                            obj.buffer(:,:,(2*ch-1)) = circshift(tempTiff(1:obj.h,:),obj.bidiShift,1);
                            obj.buffer(:,:,(2*ch)) = flipud(tempTiff((obj.htot-obj.h+1):obj.htot,:));
                        end
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
            if (obj.volume)
                proj = zeros(obj.h,obj.w,obj.slices,(obj.channels*(1+obj.bidirectional)));
            else
                proj = zeros(obj.h,obj.w,(obj.channels*(1+obj.bidirectional)));
            end
            dlg.Indeterminate = "on";
            dlg.Title = "Counting Frames...";
            if (frames==0)
                obj.getNumFrames();
                if (~obj.volume)
                    obj.frames = obj.frames/(1+obj.bidirectional);
                end
            else
                obj.frames = 1;
                upperLim = frames*obj.channels*obj.slices*(1+obj.bidirectional*obj.volume);
                while(~obj.tin.lastDirectory && obj.frames < upperLim)
                    obj.frames = obj.frames + 1;
                    obj.tin.nextDirectory();
                end
                obj.frames = obj.frames/(obj.channels*obj.slices);
            end
            dlg.Indeterminate = 'off';
            dlg.Title = 'Loading Image...';

            obj.tin.setDirectory(1);
            notify_count = round(obj.frames/20);
            obj.index = 0;
            for fr=1:obj.frames
                if (mod(fr,notify_count)==0)
                    dlg.Message = ['Frame (',num2str(fr),'/',num2str(obj.frames),')'];
                    dlg.Value = fr/obj.frames;
                end

                obj.readFrame();
                if (obj.volume)
                    for sl=1:obj.slices
                        for ch=1:obj.channels
                            if (obj.bidirectional)
                                proj(:,:,sl,2*ch-obj.frameParity) = proj(:,:,sl,2*ch-obj.frameParity) + double(obj.buffer(:,:,sl,ch))*2/obj.frames;
                            else
                                proj(:,:,sl,ch) = proj(:,:,sl,ch) + double(obj.buffer(:,:,sl,ch))/obj.frames;
                            end
                        end
                    end
                else
                    if (obj.bidirectional)
                        obj.index = obj.index+1;
                        for ch=1:obj.channels
                            proj(:,:,(2*ch-1))=proj(:,:,(2*ch-1))+double(obj.buffer(:,:,(2*ch-1)))/obj.frames;
                            proj(:,:,(2*ch))=proj(:,:,(2*ch))+double(obj.buffer(:,:,(2*ch)))/obj.frames;
                        end
                    else
                        for ch=1:obj.channels
                            proj(:,:,ch) = proj(:,:,ch)+double(obj.buffer(:,:,ch))/obj.frames;
                        end
                    end
                end
            end
            close(obj.tin);
            obj.tin = [];
        end

        % get image stored from current frame index
        function image = getImage(obj,channel)
            if (obj.bidirectional)
                if (obj.frameParity)
                    image = obj.buffer(:,:,(2*channel-1));
                else
                    image = obj.buffer(:,:,(2*channel));
                end
            else
                image = obj.buffer(:,:,channel);
            end
        end

        function image = getVolumeImage(obj,channel,slice)
            image = obj.buffer(:,:,slice,channel);
        end

        % get all images stored from current frame index
        function images = getFrame(obj)
            if (obj.bidirectional && ~obj.volume)
                if (obj.frameParity)
                    images = obj.buffer(:,:,(2*(1:obj.channels)-1));
                else
                    images = obj.buffer(:,:,(2*(1:obj.channels)));
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
            if (obj.volume)
                realFrame = frame;
                obj.index = realFrame-1;
            else
                realFrame = ceil(frame/(1+obj.bidirectional));
                obj.index = (realFrame-1)*(1+obj.bidirectional);
            end

            if (obj.index > 0)
                obj.tin.setDirectory((realFrame-1)*obj.channels*obj.slices);
            else
                obj.tin.setDirectory(1);
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
            if (obj.volume)
                obj.frames = floor(obj.frames / (obj.channels * obj.slices));
            else
                obj.frames = floor(obj.frames * (obj.bidirectional + 1) / obj.channels);
            end
            obj.tin.setDirectory(1);
        end
    end
end