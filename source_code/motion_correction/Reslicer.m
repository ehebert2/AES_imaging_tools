classdef Reslicer < handle
    properties (GetAccess = public, SetAccess = private)
        map
        imDim
        invMap

        channelsIn
        channelsOut
        slicesIn
        slicesOut
        
        mappedIms
        idMap
        mapActive
    end

    methods
        function obj = Reslicer(imDim,slices,channels)
            obj.imDim = imDim;
            obj.idMap = zeros(slices,channels,2);
            obj.channelsIn = channels;
            obj.slicesIn = slices;
            for sl = 1:slices
                for ch = 1:channels
                    obj.idMap(sl,ch,:) = [sl,ch];
                end
            end
            obj.resetMap();
        end

        function newArr = mapCells(obj,oldArr)
            newArr = cell(obj.slicesOut,obj.channelsOut);
            for sl=1:obj.slicesIn
                for ch=1:obj.channelsIn
                    newArr{obj.map(sl,ch,1),obj.map(sl,ch,2)} = oldArr{sl,ch};
                end
            end
        end

        function newIms = mapImages(obj,ims)
            newIms = obj.mappedIms;
            for sl=1:obj.slicesIn
                for ch=1:obj.channelsIn
                    newIms(:,:,obj.map(sl,ch,1),obj.map(sl,ch,2)) = ims(:,:,sl,ch);
                end
            end
        end

        function B = mapMatrix(obj,A)
            B = zeros(obj.slicesOut,obj.channelsOut);
            for sl=1:obj.slicesIn
                for ch=1:obj.channelsIn
                    B(obj.map(sl,ch,1),obj.map(sl,ch,2)) = A(sl,ch);
                end
            end
        end

        function coord = getCoord(obj,sl,ch)
            coord = obj.map(sl,ch,:);
        end

        function coord = getInvCoord(obj,sl,ch)
            coord = obj.invMap(sl,ch,:);
        end

        function setMap(obj,map)
            if ((size(map,1)~=obj.slicesIn)||(size(map,2)~=obj.channelsIn))
                return;
            end

            obj.map = map;
            obj.slicesOut = max(obj.map(:,:,1),[],'all');
            obj.channelsOut = max(obj.map(:,:,2),[],'all');
            obj.invMap = zeros(obj.slicesOut,obj.channelsOut,2);
            for sl2=1:obj.slicesOut
                for ch2=1:obj.channelsOut
                    for sl1=1:obj.slicesIn
                        for ch1=1:obj.channelsIn
                            if ((obj.map(sl1,ch1,1) == sl2) && (obj.map(sl1,ch1,2) == ch2))
                                obj.invMap(sl2,ch2,:) = [sl1,ch1];
                            end
                        end
                    end
                end
            end
            obj.mappedIms = int16(zeros(obj.imDim(1),obj.imDim(2),obj.slicesOut,obj.channelsOut));
            obj.mapActive = true;
        end
        
        function setInverseMap(obj,invMap)
            if ((max(invMap(:,:,1),[],'all')~=obj.slicesIn)||(max(invMap(:,:,2),[],'all')~=obj.channelsIn))
                return;
            end

            obj.slicesOut = size(invMap,1);
            obj.channelsOut = size(invMap,2);
            obj.invMap = invMap;
            obj.map = zeros(obj.slicesIn,obj.channelsIn,2);
            for sl1 = 1:obj.slicesIn
                for ch1 = 1:obj.channelsIn
                    for sl2 = 1:obj.slicesOut
                        for ch2 = 1:obj.channelsOut
                            if ((invMap(sl2,ch2,1)==sl1)&&(invMap(sl2,ch2,2)==ch1))
                                obj.map(sl1,ch1,:) = [sl2,ch2];
                            end
                        end
                    end
                end
            end
            obj.mappedIms = int16(zeros(obj.imDim(1),obj.imDim(2),obj.slicesOut,obj.channelsOut));
            obj.mapActive = true;
        end

        function resetMap(obj)
            obj.map = obj.idMap;
            obj.invMap = obj.idMap;
            obj.slicesOut = obj.slicesIn;
            obj.channelsOut = obj.channelsIn;
            obj.mappedIms = int16(zeros(obj.imDim(1),obj.imDim(2),obj.slicesOut,obj.channelsOut));
            obj.mapActive = false;
        end

        function params = attachParams(obj, params)
            params.reslice = obj.mapActive;
            params.resliceMap = obj.map;
            params.slicesOut = obj.slicesOut;
            params.channelsOut = obj.channelsOut;
        end
    end
end