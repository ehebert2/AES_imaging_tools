classdef StackTracker < handle
    properties
        channels
        bidirectional
        bidiShift
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

    methods
        function obj = StackTracker(params)
            params = validateParams(params);
            obj.channels = params.channels;
            obj.bidirectional = params.bidirectional;
            obj.bidiShift = params.bidiShift;
        end

        function open(obj,filename)
            obj.filename = filename;
            obj.tin = Tiff(filename,'r');
            tempTiff = read(obj.tin);
            if (obj.bidirectional)
                obj.htot = size(tempTiff,1);
                obj.h = floor(obj.htot/2);
                obj.buffer = int16(zeros(obj.h,size(tempTiff,2),(2*obj.channels)));
            else
                obj.buffer = int16(zeros(size(tempTiff,1),size(tempTiff,2),obj.channels));
                obj.h = size(tempTiff,1);
            end
            obj.w = size(tempTiff,2);
            obj.index = 0;
            obj.frames = 0;
        end

        function readFrame(obj)
            obj.index = obj.index+1;
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

        function close(obj)
            close(obj.tin);
        end

        function proj = buildProjection(obj,filename,frames,dlg)
            obj.open(filename);
            proj = zeros(obj.h,obj.w,(obj.channels*(1+obj.bidirectional)));
            if (frames==0)
                obj.getNumFrames();
                obj.frames = obj.frames/(1+obj.bidirectional);
            else
                obj.frames = 1;
                while(~obj.tin.lastDirectory && obj.frames < (frames*obj.channels))
                    obj.frames = obj.frames + 1;
                    obj.tin.nextDirectory();
                end
                obj.frames = obj.frames/obj.channels;
            end

            obj.tin.setDirectory(1);
            notify_count = round(obj.frames/20);
            obj.index = 0;
            for fr=1:obj.frames
                if (mod(fr,notify_count)==0)
                    dlg.Value = fr/obj.frames;
                end

                obj.readFrame();
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
            close(obj.tin);
            obj.tin = [];
        end

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

        function images = getFrame(obj)
            if (obj.bidirectional)
                if (obj.frameParity)
                    images = obj.buffer(:,:,(2*(1:obj.channels)-1));
                else
                    images = obj.buffer(:,:,(2*(1:obj.channels)));
                end
            else
                images = obj.buffer;
            end
        end

        function frames = getNumFrames(obj)
            if (obj.frames < 1)
                obj.countFrames();
            end

            frames = obj.frames;
        end

        function setIndex(obj, frame)
            realFrame = ceil(frame/(1+obj.bidirectional));
            obj.index = (realFrame-1)*(1+obj.bidirectional);
            if (obj.index > 0)
                obj.tin.setDirectory((realFrame-1)*obj.channels);
            else
                obj.tin.setDirectory(1);
            end
        end
    end

    methods (Access = private)
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
    end
end