classdef MaskFiller < handle
    properties (GetAccess = private, SetAccess = private)
        minWidth
        imSwitches
        fillerMasks

        imDim
        channels
        channel
        sepChannels
        slices
        slice

        baseMask
        line
        lastEdge
        polar
        width
        baseDspLeft
        baseDspRight
        dspLeft
        dspRight
        bndLeft
        bndRight
        diff

        threshMaskHandler
        shapeMaskHandler
    end

    properties (SetAccess = public, GetAccess = public)
        visible
    end

    methods
        function obj = MaskFiller(params,threshMaskHandler, shapeMaskHandler)
            obj.threshMaskHandler = threshMaskHandler;
            obj.shapeMaskHandler = shapeMaskHandler;

            obj.imDim = [size(obj.threshMaskHandler.fullMask,1),size(obj.threshMaskHandler.fullMask,2)];
            obj.channels = params.channels*params.splitChannels + (1-params.splitChannels);
            obj.sepChannels = params.splitChannels;
            obj.channel = 1;
            obj.slices = params.slices*params.volume + (1-params.volume);
            obj.slice = 1;

            obj.fillerMasks = zeros(obj.imDim(1),obj.imDim(2),obj.slices,obj.channels);
            obj.imSwitches = zeros(obj.slices,obj.channels);
            obj.minWidth = 1;
            obj.visible = true;
        end

        function setChannel(obj,channel)
            if (obj.sepChannels)
                obj.channel = channel;
            end
        end

        function setSlice(obj,slice)
            obj.slice = slice;
        end

        function setMinWidth(obj,minWidth)
            if (obj.minWidth~=minWidth)
                obj.minWidth = minWidth;
                obj.baseDspLeft = floor(minWidth/2);
                obj.baseDspRight = ceil(minWidth/2)-1;
                obj.updateMaskTot();
            end
        end

        function mask = getMask(obj)
            mask = (obj.fillerMasks(:,:,obj.slice,obj.channel)*obj.visible);
        end

        function mask = getChMask(obj,sl,ch)
            mask = obj.fillerMasks(:,:,sl,ch);
        end

        function sw = getSwIm(obj)
            sw = obj.imSwitches(obj.slice,obj.channel);
        end

        function sw = getSwCh(obj)
            sw = sum(obj.imSwitches(:,obj.channel));
        end

        function sw = getSwTot(obj)
            sw = sum(obj.imSwitches,'all');
        end

        function updateMask(obj)
            obj.fill(obj.slice,obj.channel);
        end

        function updateMaskCh(obj)
            for sl=1:obj.slices
                obj.fill(sl,obj.channel);
            end
        end

        function updateMaskTot(obj)
            for sl=1:obj.slices
                for ch=1:obj.channels
                    obj.fill(sl,ch);
                end
            end
        end

        function fill(obj,sl,ch)
            obj.baseMask = (obj.threshMaskHandler.getFullMask(sl,ch) + obj.shapeMaskHandler.getChMask(sl,ch))>0;
            obj.fillerMasks(:,:,sl,ch) = 0*obj.fillerMasks(:,:,sl,ch);
            if ((obj.minWidth<2) || (sum(obj.baseMask,'all')==0))
                edges = (obj.baseMask(:,2:end) > obj.baseMask(:,1:(end-1)));
                obj.imSwitches(sl,ch) = sum(edges,'all') + sum(obj.baseMask(:,1));
                return;
            end

            obj.imSwitches(sl,ch) = 0;
            for ln = 1:obj.imDim(1)
                obj.line = obj.baseMask(ln,:);
                if (sum(obj.line)==0)
                    continue;
                end
                obj.polar = obj.line(1);
                if (obj.polar)
                    obj.lastEdge = 1;
                end

                px = 1;
                while (px<obj.imDim(2))
                    px = px+1;
                    if (obj.line(px) && ~obj.polar)
                        obj.polar = true;
                        obj.lastEdge = px;
                    elseif (obj.polar && ~obj.line(px))
                        obj.polar = false;
                        obj.width = px-obj.lastEdge;
                        if (obj.width < obj.minWidth)
                            obj.diff = obj.minWidth-obj.width;
                            obj.bndLeft = obj.lastEdge-obj.baseDspLeft+ceil(obj.width/2);
                            obj.bndRight = px+obj.baseDspRight-floor(obj.width/2);
                            if (obj.bndLeft<1)
                                obj.bndRight=obj.bndRight-obj.bndLeft+1;
                                obj.bndLeft = 1;
                            elseif (obj.bndRight > obj.imDim(2))
                                obj.bndLeft = obj.bndLeft - (obj.bndRight-obj.imDim(2));
                                obj.bndRight = obj.imDim(2);
                            end

                            obj.line(obj.bndLeft:obj.bndRight) = true;
                            px = obj.bndRight;
                        end
                    end
                end

                if (obj.polar)
                    obj.width = px-obj.lastEdge+1;
                    if (obj.width < obj.minWidth)
                        obj.bndLeft = obj.imDim(2)-obj.minWidth+1;
                        obj.line(obj.bndLeft:obj.imDim(2)) = true;
                    end
                end

                px = 1;
                while(~obj.line(px))
                    px = px+1;
                end
                obj.polar = true;
                while (px < obj.imDim(2))
                    px = px+1;
                    if (obj.line(px)&&~obj.polar)
                        obj.polar = true;
                        obj.width = px-obj.lastEdge;
                        if (obj.width<obj.minWidth)
                            obj.line(obj.lastEdge:(px-1)) = true;
                            obj.imSwitches(sl,ch) = obj.imSwitches(sl,ch)-1;
                        end
                    elseif (obj.polar&&~obj.line(px))
                        obj.polar = false;
                        obj.lastEdge = px;
                        obj.imSwitches(sl,ch) = obj.imSwitches(sl,ch)+1;
                    end
                end

                obj.imSwitches(sl,ch) = obj.imSwitches(sl,ch)+obj.line(px);
                obj.fillerMasks(ln,:,sl,ch) = obj.line;
            end

            obj.fillerMasks(:,:,sl,ch) = obj.fillerMasks(:,:,sl,ch) - obj.baseMask;
        end
    end
end