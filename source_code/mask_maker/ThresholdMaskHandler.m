classdef ThresholdMaskHandler < handle
    properties (GetAccess = private, SetAccess = private)
        callingApp
        
        mainChannel
        channels
        slices
        yShift
        bidirectional
        sepChannels
        maskGenerators

        chngAllSlices
        blockOverlap
        maxHole
        compMasks
    end

    properties (GetAccess = public, SetAccess = private)
        baseMask
        fullMask
        selectedMask
        excludedMask
        channel
        slice
    end

    methods
        function obj = ThresholdMaskHandler(callingApp,params,images)
            obj.callingApp = callingApp;
            obj.sepChannels = params.splitChannels;
            if (obj.sepChannels)
                obj.channels = params.channels;
            else
                obj.channels = 1;
                obj.mainChannel = params.mainChannel;
            end

            if (params.volume)
                obj.slices = params.slices;
            else
                obj.slices = 1;
            end

            obj.channel = 1;
            obj.slice = 1;
            obj.baseMask = zeros(size(images,1),size(images,2));
            obj.fullMask = zeros(size(obj.baseMask));
            obj.selectedMask = zeros(size(obj.baseMask));
            obj.excludedMask = zeros(size(obj.baseMask));
            obj.compMasks = zeros(size(obj.baseMask,1),size(obj.baseMask,2),obj.slices);

            obj.chngAllSlices = false;
            obj.blockOverlap = false;
            obj.bidirectional = params.bidirectional;
            if (obj.bidirectional)
                obj.yShift = params.bidiShift;
            else
                obj.yShift = 0;
            end

            obj.maskGenerators = cell(obj.slices,obj.channels);
            for sl=1:obj.slices
                if (obj.sepChannels)
                    for ch=1:obj.channels
                        obj.maskGenerators{sl,ch} = MaskGenerator(images(:,:,sl,(1+obj.bidirectional)*ch-obj.bidirectional));
                    end
                else
                    obj.maskGenerators{sl,obj.mainChannel} = MaskGenerator(images(:,:,sl,(1+obj.bidirectional)*obj.mainChannel-obj.bidirectional));
                end
            end
        end

        %% Setters
        function setChannel(obj, channel)
            if (obj.sepChannels)
                obj.deselectROI();
                obj.maskGenerators{obj.slice,obj.channel}.deselectROI();
                obj.channel = channel;
                obj.updateMasks();
            end
        end

        function setSlice(obj, slice)
            obj.deselectROI();
            obj.slice = slice;
            obj.updateMasks();
        end

        function setYShift(obj,yShift)
            obj.yShift = yShift;
            obj.updateMasks();
        end

        function setBorder(obj,border)
            if (obj.chngAllSlices)
                for sl=1:obj.slices
                    obj.maskGenerators{sl,obj.channel}.setBorder(border);
                    obj.maskGenerators{sl,obj.channel}.buildBorder();
                end
                obj.callingApp.updateChOcc();
            else
                obj.maskGenerators{obj.slice,obj.channel}.setBorder(border);
                obj.maskGenerators{obj.slice,obj.channel}.buildBorder();
                obj.callingApp.updateOcc();
            end

            if (obj.bidirectional)
                obj.fullMask = circshift(obj.maskGenerators{obj.slice,obj.channel}.fullMask,obj.yShift,1);
            else
                obj.fullMask = obj.maskGenerators{obj.slice,obj.channel}.fullMask;
            end
        end

        function setMaxHole(obj,maxHole)
            obj.deselectROI();
            for sl=1:obj.slices
                for ch=1:obj.channels
                    obj.maskGenerators{sl,ch}.maxHole = maxHole;
                    obj.maskGenerators{sl,ch}.buildBaseMask();
                end
            end
            
            if (obj.blockOverlap)
                for ch=1:obj.channels
                    obj.removeMaskOverlap(ch);
                end
            else
                for sl=1:obj.slices
                    for ch=1:obj.channels
                        obj.maskGenerators{sl,ch}.buildBorder();
                    end
                end
            end

            obj.maxHole = maxHole;
            obj.updateMasks();
            obj.callingApp.updateTotOcc();
        end

        function setImages(obj,images)
            if (obj.sepChannels)
                for sl=1:obj.slices
                    for ch=1:obj.channels
                        obj.maskGenerators{sl,ch}.setImage(images(:,:,sl,(1+obj.bidirectional)*ch-obj.bidirectional));
                    end
                end
            else
                for sl=1:obj.slices
                    obj.maskGenerators{sl,1}.setImage(images(:,:,sl,(1+obj.bidirectional)*obj.mainChannel-obj.bidirectional));
                end
            end

            if (obj.blockOverlap)
                for ch=1:obj.channels
                    obj.removeMaskOverlap(ch);
                end
            else
                for sl=1:obj.slices
                    for ch=1:obj.channels
                        obj.maskGenerators{sl,ch}.buildBorder();
                    end
                end
            end

            obj.updateMasks();
            obj.callingApp.updateTotOcc();
        end

        function setBlockOverlap(obj,blockOverlap)
            if (obj.blockOverlap==blockOverlap)
                return;
            else
                obj.blockOverlap = blockOverlap;
                if (obj.blockOverlap)
                    for ch=1:obj.channels
                        obj.removeMaskOverlap(ch);
                    end
                else
                    for sl=1:obj.slices
                        for ch=1:obj.channels
                            obj.maskGenerators{sl,ch}.quickFilter();
                            obj.maskGenerators{sl,ch}.buildBorder();
                        end
                    end
                end
            end
            obj.updateMasks();
            obj.callingApp.updateTotOcc();
        end

        function setChngAllSlices(obj,chngAllSlices)
            obj.chngAllSlices = chngAllSlices;
        end

        %% Getters
        function border = getBorder(obj)
            border = obj.maskGenerators{obj.slice,obj.channel}.border;
        end

        function threshold = getThreshold(obj)
            threshold = obj.maskGenerators{obj.slice,obj.channel}.threshold;
        end

        function gen = getMaskGenerator(obj)
            gen = obj.maskGenerators{obj.slice,obj.channel};
        end

        function mask = getFullMask(obj,sl,ch)
            mask = obj.maskGenerators{sl,ch}.fullMask;
        end

        function value = getMaxHole(obj)
            value = obj.maskGenerators{obj.slice,obj.channel}.maxHole;
        end

        function value = getChngAllSlices(obj)
            value = obj.chngAllSlices;
        end

        function value = getBlockOverlap(obj)
            value = obj.blockOverlap;
        end

        %% threshold and min feature size modifying functions
        function chngThresholdFast(obj, threshold)
            if (obj.bidirectional)
                obj.baseMask = circshift(obj.maskGenerators{obj.slice,obj.channel}.image > threshold,obj.yShift,1);
            else
                obj.baseMask = obj.maskGenerators{obj.slice,obj.channel}.image > threshold;
            end
            obj.excludedMask = 0*obj.excludedMask;
            obj.selectedMask = 0*obj.selectedMask;
            obj.fullMask = obj.baseMask;
        end

        function chngThreshold(obj, thresholdType, threshold)
            if (obj.chngAllSlices)
                for sl=1:obj.slices
                    obj.maskGenerators{sl,obj.channel}.thresholdType = thresholdType;
                    switch (thresholdType)
                        case obj.maskGenerators{1,1}.GLOBAL
                            obj.maskGenerators{sl,obj.channel}.threshold = threshold;
                        case obj.maskGenerators{1,1}.ADAPTIVE
                            obj.maskGenerators{sl,obj.channel}.sensitivity = threshold;
                    end
                    obj.maskGenerators{sl,obj.channel}.buildBaseMask();
                end

                if (obj.blockOverlap)
                    obj.removeMaskOverlap(obj.channel);
                else
                    for sl=1:obj.slices
                        obj.maskGenerators{sl,obj.channel}.buildBorder();
                    end
                end
                obj.callingApp.updateChOcc();
            else
                obj.maskGenerators{obj.slice,obj.channel}.thresholdType = thresholdType;
                switch (thresholdType)
                    case obj.maskGenerators{1,1}.GLOBAL
                        obj.maskGenerators{obj.slice,obj.channel}.threshold = threshold;
                    case obj.maskGenerators{1,1}.ADAPTIVE
                        obj.maskGenerators{obj.slice,obj.channel}.sensitivity = threshold;
                end
                obj.maskGenerators{obj.slice,obj.channel}.buildBaseMask();

                if (obj.blockOverlap)
                    obj.removeMaskOverlap(obj.channel);
                else
                    obj.maskGenerators{obj.slice,obj.channel}.buildBorder();
                end
                obj.callingApp.updateOcc();
            end
            obj.maskGenerators{obj.slice,obj.channel}.deselectROI();
            obj.updateMasks();
        end

        function chngSizeFast(obj, minFeature)
            obj.maskGenerators{obj.slice,obj.channel}.minFeature = minFeature;
            obj.maskGenerators{obj.slice,obj.channel}.quickFilter();
            if (obj.bidirectional)
                obj.baseMask = circshift(obj.maskGenerators{obj.slice,obj.channel}.baseMask,obj.yShift,1);
                obj.excludedMask = circshift(obj.maskGenerators{obj.slice,obj.channel}.excludedMask,obj.yShift,1);
            else
                obj.baseMask = obj.maskGenerators{obj.slice,obj.channel}.baseMask;
                obj.excludedMask = obj.maskGenerators{obj.slice,obj.channel}.excludedMask;
            end
            obj.fullMask = obj.baseMask;
            obj.selectedMask = 0*obj.selectedMask;
        end

        function chngSize(obj,minFeature)
            if (obj.chngAllSlices)
                for sl=1:obj.slices
                    obj.maskGenerators{sl,obj.channel}.minFeature = minFeature;
                    obj.maskGenerators{sl,obj.channel}.quickFilter();
                end

                if (obj.blockOverlap)
                    obj.removeMaskOverlap(obj.channel);
                else
                    for sl=1:obj.slices
                        obj.maskGenerators{sl,obj.channel}.buildBorder();
                    end
                end
                obj.callingApp.updateChOcc();
            else
                obj.maskGenerators{obj.slice,obj.channel}.minFeature = minFeature;
                obj.maskGenerators{obj.slice,obj.channel}.quickFilter();
                if (obj.blockOverlap)
                    obj.removeMaskOverlap(obj.channel);
                else
                    obj.maskGenerators{obj.slice,obj.channel}.buildBorder();
                end
                obj.callingApp.updateOcc();
            end

            obj.maskGenerators{obj.slice,obj.channel}.deselectROI();
            obj.updateMasks();
        end

        %% Other functions
        function removeMaskOverlap(obj,ch)
            for sl=1:obj.slices
                obj.compMasks(:,:,sl) = obj.maskGenerators{sl,ch}.filterMask .* obj.maskGenerators{sl,ch}.image;
            end

            maxVals = repmat(max(obj.compMasks,[],3),1,1,obj.slices);
            obj.compMasks = (obj.compMasks == maxVals).*(obj.compMasks>0);

            for sl=1:obj.slices
                obj.maskGenerators{sl,ch}.baseMask = obj.maskGenerators{sl,ch}.filterMask .* obj.compMasks(:,:,sl);
                obj.maskGenerators{sl,ch}.excludedMask = obj.maskGenerators{sl,ch}.filterMask + obj.maskGenerators{sl,ch}.excludedFilterMask - obj.maskGenerators{sl,ch}.baseMask;
                obj.maskGenerators{sl,ch}.buildBorder();
            end
        end

        function updateMasks(obj)
            if (obj.bidirectional)
                obj.baseMask = circshift(obj.maskGenerators{obj.slice,obj.channel}.baseMask,obj.yShift,1);
                obj.fullMask = circshift(obj.maskGenerators{obj.slice,obj.channel}.fullMask,obj.yShift,1);
                obj.excludedMask = circshift(obj.maskGenerators{obj.slice,obj.channel}.excludedMask,obj.yShift,1);
                obj.selectedMask = circshift(obj.maskGenerators{obj.slice,obj.channel}.selectedMask,obj.yShift,1);
            else
                obj.baseMask = obj.maskGenerators{obj.slice,obj.channel}.baseMask;
                obj.fullMask = obj.maskGenerators{obj.slice,obj.channel}.fullMask;
                obj.excludedMask = obj.maskGenerators{obj.slice,obj.channel}.excludedMask;
                obj.selectedMask = obj.maskGenerators{obj.slice,obj.channel}.selectedMask;
            end
        end

        function hit = imgClk(obj,x,y,hold)
            obj.maskGenerators{obj.slice,obj.channel}.clckDetect(x,(y-obj.bidirectional*obj.yShift),hold);
            hit = obj.maskGenerators{obj.slice,obj.channel}.roiSelected;
            if (hit)
                if (obj.bidirectional)
                    obj.selectedMask = circshift(obj.maskGenerators{obj.slice,obj.channel}.selectedMask,obj.yShift,1);
                else
                    obj.selectedMask = obj.maskGenerators{obj.slice,obj.channel}.selectedMask;
                end
            else
                obj.selectedMask = 0*obj.selectedMask;
            end
        end

        function deselectROI(obj)
            obj.maskGenerators{obj.slice,obj.channel}.deselectROI();
            obj.selectedMask = 0*obj.selectedMask;
        end

        function deleteROI(obj)
            obj.maskGenerators{obj.slice,obj.channel}.deleteROI();
            obj.updateMasks();
        end

        function printROIMasks(obj,basename,path,ind)
            for ch=1:obj.channels
                if (obj.channels>1)
                    folderTitle = strcat('channel_',num2str(ch));
                    chPath = fullfile(path,folderTitle);
                else
                    chPath = path;
                end

                for sl=1:obj.slices
                    if (obj.slices>1)
                        folderTitle = strcat('slice_',num2str(sl));
                        slPath = fullfile(chPath,folderTitle);
                    else
                        slPath = chPath;
                    end

                    obj.maskGenerators{sl,ch}.printBlobs(basename,slPath,ind(sl,ch));
                end
            end
        end

        function saveROIs(obj,fname)
            realCh = obj.channels*obj.sepChannels+1-obj.sepChannels;
            threshMaskParams.chngAllSlices = obj.chngAllSlices;
            threshMaskParams.blockOverlap  = obj.blockOverlap;
            threshMaskParams.maskData = cell(obj.slices,realCh);
            obj.maskGenerators{obj.slice,obj.channel}.deselectROI();
            for sl=1:obj.slices
                for ch=1:realCh
                    threshMaskParams.maskData{sl,ch} = obj.maskGenerators{sl,ch}.save();
                end
            end
            save(fname,"threshMaskParams",'-append');
        end

        function loadROIs(obj,data)
            if (~isfield(data,'threshMaskParams'))
                return;
            end

            realCh = obj.channels*obj.sepChannels+1-obj.sepChannels;
            obj.maskGenerators{obj.slice,obj.channel}.deselectROI();
            params = data.threshMaskParams;
            for sl=1:obj.slices
                for ch=1:realCh
                    obj.maskGenerators{sl,ch}.load(params.maskData{sl,ch});
                end
            end
            obj.maxHole = obj.maskGenerators{1,1}.maxHole;
            obj.chngAllSlices = params.chngAllSlices;
            obj.blockOverlap = params.blockOverlap;
            if (obj.blockOverlap)
                for ch=1:realCh
                    obj.removeMaskOverlap(ch);
                end
            end

            for sl=1:obj.slices
                for ch=1:realCh
                    obj.maskGenerators{sl,ch}.buildBorder();
                end
            end

            obj.updateMasks();
        end
    end
end