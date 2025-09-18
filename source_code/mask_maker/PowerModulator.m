% model for managing analogue power modulation
% calibration map is Nx2 matrix. Column 1 is modulator voltage (pixel
% value) and column 2 is associated transmission

classdef PowerModulator < handle
    properties (GetAccess = private, SetAccess = private)    
        sliceHighPx
        sliceLowPx
        pxVals        
    end

    properties (GetAccess = public, SetAccess = private)
        UNIFORM = 1;
        EXP = 2;
        INVPOW = 3;

        mainChannel
        channels
        bidirectional
        slices
        sepChannels
        imDim
        
        % FOV power modulation
        fovTrans
        fovPreTrans
        fovTransNames
        fovPreTransNames

        % ROI power modulation
        roiRefPower
        roiTrans
        roiThresh
        roiThreshPxVal
        roiLowPxVal
        roiHighPxVal
        shapeMaskHandler
        threshMaskHandler
        
        % slice power modulation
        sliceRefPower
        sliceTrans
        slicePreTrans
        slicePreTransNames
        
        % variables for measuring slice power
        curveTypes
        z
        mpOrder
        sliceSep
        attenLen
        orientation
        numPx

        minTrans
        absMinTrans
        calMaps
        calNames

        % others
        roiUI
        sliceUI
        fovUI
        refImages
        imFactors
    end

    methods
        function obj = PowerModulator(params,ims,shapeMaskHandler,threshMaskHandler)
            % stack structure variables
            obj.sepChannels = params.splitChannels && (params.channels>1);
            obj.channels = params.channels;
            obj.mainChannel = params.mainChannel;
            obj.slices = params.slices;
            obj.imDim = [size(ims,1),size(ims,2)];
            obj.numPx = obj.imDim(1)*obj.imDim(2);
            obj.bidirectional = params.bidirectional;
            realCh = obj.channels*obj.sepChannels+1-obj.sepChannels;
            obj.imFactors = ones(obj.imDim(1),obj.imDim(2),obj.slices,(obj.channels*obj.sepChannels+1-obj.sepChannels));

            % volume power measurement variables
            obj.pxVals = zeros(obj.numPx,obj.slices,realCh);
            obj.mpOrder = 2;
            obj.sliceSep = 1;
            obj.z = 0:(obj.slices-1);
            obj.orientation = ones(obj.channels,1);
            obj.attenLen = ones(obj.channels,1);
            obj.sliceHighPx = round(1*obj.numPx);
            obj.sliceLowPx = round(0.9*obj.numPx);
            obj.curveTypes = obj.UNIFORM * ones(obj.channels,1);
            obj.sliceTrans = ones(obj.slices,realCh);
            obj.sliceRefPower = zeros(obj.slices,realCh);
            obj.slicePreTransNames = cell(realCh,1);
            obj.slicePreTrans = ones(obj.slices,realCh);

            % roi modulation variables
            obj.shapeMaskHandler = shapeMaskHandler;
            obj.shapeMaskHandler.powerMod = obj;
            obj.roiTrans = cell(obj.slices,realCh);
            obj.roiRefPower = cell(obj.slices,realCh);
            obj.roiThreshPxVal = cell(obj.slices,realCh);
            obj.roiLowPxVal = cell(obj.slices,realCh);
            obj.roiHighPxVal = cell(obj.slices,realCh);
            obj.roiThresh = 0.95;

            % fov modulation variables
            obj.fovTrans = ones(obj.imDim(1),obj.imDim(2),realCh);
            obj.fovPreTrans = ones(obj.imDim(1),obj.imDim(2),realCh);
            obj.fovTransNames = cell(realCh,1);
            obj.fovPreTransNames = cell(realCh,1);

            % variables to calculate final output
            obj.calMaps = cell(realCh,1);
            obj.calNames = cell(realCh,1);
            obj.minTrans = 0;
            obj.absMinTrans = 0;

            obj.refImages = zeros(obj.imDim(1),obj.imDim(2),obj.slices,obj.channels);
            obj.threshMaskHandler = threshMaskHandler;
            for ch=1:(obj.channels*obj.sepChannels+1-obj.sepChannels)
                for sl=1:obj.slices
                    obj.refImages(:,:,sl,ch) = ims(:,:,sl,(ch*(1+obj.bidirectional)));
                    obj.pxVals(:,sl,ch) = sort(reshape(obj.refImages(:,:,sl,ch),[],1));
                end
            end
            obj.calcSliceRefPower();
        end

        function setImages(obj,ims)
            if (obj.sepChannels)
                for ch=1:obj.channels
                    for sl=1:obj.slices
                        obj.refImages(:,:,sl,ch) = ims(:,:,sl,(ch*(1+obj.bidirectional)));
                        obj.pxVals(:,sl,ch) = sort(reshape(obj.refImages(:,:,sl,ch),[],1));
                    end
                end
            else
                for sl=1:obj.slices
                    obj.refImages(:,:,sl,1) = ims(:,:,sl,obj.mainChannel*(1+obj.bidirectional));
                    obj.pxVals(:,sl,1) = sort(reshape(obj.refImages(:,:,sl,1),[],1));
                end
            end

            obj.threshMaskHandler.setImages(obj.refImages.*repmat(permute(((obj.fovTrans./obj.fovPreTrans).^obj.mpOrder),[1,2,4,3]),1,1,obj.slices,1));
            obj.calcSliceRefPower();
            obj.calcAllSliceTrans();
            if (~isempty(obj.roiUI))
                obj.calcAllROIRefPower();
                obj.roiUI.chngSlCh();
            end
        end

        %% functions related to volume modulation
        function calcSliceRefPower(obj)
            for ch=1:(obj.channels*obj.sepChannels+1-obj.sepChannels)
                for sl=1:obj.slices
                    obj.sliceRefPower(sl,ch) = mean(obj.pxVals(obj.sliceLowPx:obj.sliceHighPx,sl,ch));
                end
                obj.sliceRefPower(:,ch) = obj.sliceRefPower(:,ch) / max(obj.sliceRefPower(:,ch));
            end
        end

        function calcSliceTrans(obj,ch)
            ch = ch*obj.sepChannels+1-obj.sepChannels;
            switch(obj.curveTypes(ch))
                case (obj.UNIFORM)
                    obj.sliceTrans(:,ch) = ones(obj.slices,1);
                case (obj.EXP)
                    obj.sliceTrans(:,ch) = exp(obj.orientation(ch)*(obj.z')/(obj.attenLen(ch)*obj.mpOrder));
                    obj.sliceTrans(:,ch) = obj.sliceTrans(:,ch)/max(obj.sliceTrans(:,ch));
                case (obj.INVPOW)
                    obj.sliceTrans(:,ch) = obj.slicePreTrans(:,ch)./(obj.sliceRefPower(:,ch).^(1/obj.mpOrder));
                    obj.sliceTrans(:,ch) = obj.sliceTrans(:,ch)/max(obj.sliceTrans(:,ch));
            end

            for sl=1:obj.slices
                obj.calcImFactor(sl,(ch*obj.sepChannels+1-obj.sepChannels));
            end
        end

        function calcAllSliceTrans(obj)
            if obj.sepChannels
                for ch=1:obj.channels
                    obj.calcSliceTrans(ch);
                end
            else
                obj.calcSliceTrans(1);
            end
        end

        function fitExp(obj,ch)
            X = [ones(obj.slices,1), obj.z'];
            Y = log(obj.sliceRefPower(:,ch)./(obj.slicePreTrans(:,ch).^obj.mpOrder));
            C = inv(X'*X)*X'*Y;
            if (C(2)~=0)
                obj.orientation(ch) = 1-2*(C(2)>0);
                obj.attenLen(ch) = 1/abs(C(2));
            end
            obj.calcSliceTrans(ch);
        end

        function power = getSlicePower(obj,ch)
            power = obj.sliceRefPower(:,ch);
        end

        function power = getSlicePrePower(obj,ch)
            power = obj.slicePreTrans(:,ch);
        end

        function t = getSliceTrans(obj,ch)
            t = obj.sliceTrans(:,ch);
        end

        function val = getSliceHighPx(obj)
            val = obj.sliceHighPx/obj.numPx;
        end

        function val = getSliceLowPx(obj)
            val = obj.sliceLowPx/obj.numPx;
        end

        function setSliceHighPx(obj,val)
            obj.sliceHighPx = round(obj.numPx*val);
            obj.calcSliceRefPower();
            obj.calcAllSliceTrans();
        end

        function setSliceLowPx(obj,val)
            obj.sliceLowPx = round(obj.numPx*val);
            if (obj.sliceLowPx == 0)
                obj.sliceLowPx = 1;
            end
            obj.calcSliceRefPower();
            obj.calcAllSliceTrans();
        end

        function setSliceSep(obj,val)
            obj.sliceSep = val;
            obj.z = linspace(0,(obj.sliceSep*(obj.slices-1)),obj.sliceSep);
            obj.calcAllSliceTrans();
        end

        function setOrientation(obj,ch,val)
            obj.orientation(ch) = val;
            obj.calcAllSliceTrans();
        end

        function setType(obj,ch,val)
            obj.curveTypes(ch) = val;
            obj.calcSliceTrans(ch);
        end

        function setSL(obj,ch,val)
            obj.attenLen(ch) = val;
            obj.calcSliceTrans(ch);
        end

        function success = setSlicePreTrans(obj,ch,curve,name)
            ch = ch*obj.sepChannels+1-obj.sepChannels;
            success = false;
            
            if (isempty(curve))
                obj.slicePreTransNames{ch} = [];
                obj.slicePreTrans(ch,:) = ones(1,obj.slices);
            else
                curve = curve(:,1);
                if ((length(curve)~=obj.slices) || (sum(isnan(curve))>0))
                    return;
                else
                    curve = abs(curve);
                    obj.slicePreTrans(:,ch) = curve/max(curve);
                    obj.slicePreTransNames{ch} = name;
                end
            end

            obj.calcSliceTrans(ch);
            success = true;
        end

        %% FOV Modulation functions
        function setFOVTrans(obj,ch,fovTrans,name)
            realCh = ch*obj.sepChannels+1-obj.sepChannels;
            if (isempty(fovTrans))
                obj.fovTrans(:,:,realCh) = ones(obj.imDim(1),obj.imDim(2));
                obj.fovTransNames{realCh} = [];
            else
                fovTrans = abs(fovTrans);
                fovTrans = imresize(fovTrans,obj.imDim);
                obj.fovTrans(:,:,realCh) = fovTrans / max(fovTrans,[],'all');
                obj.fovTransNames{realCh} = name;
            end

            for sl=1:obj.slices
                obj.calcImFactor(sl,realCh);
            end
            obj.threshMaskHandler.setImages(obj.refImages.*repmat(permute(((obj.fovTrans./obj.fovPreTrans).^obj.mpOrder),[1,2,4,3]),1,1,obj.slices,1));  
        end

        function setFOVPreTrans(obj,ch,fovTrans,name)
            realCh = ch*obj.sepChannels+1-obj.sepChannels;
            if (isempty(fovTrans))
                obj.fovPreTrans(:,:,realCh) = ones(obj.imDim(1),obj.imDim(2));
                obj.fovPreTransNames{realCh} = [];
            else
                fovTrans = abs(fovTrans);
                fovTrans = imresize(fovTrans,obj.imDim);
                obj.fovPreTrans(:,:,realCh) = fovTrans / max(fovTrans,[],'all');
                obj.fovPreTransNames{realCh} = name;
            end

            for sl=1:obj.slices
                obj.calcImFactor(sl,realCh);
            end

            obj.threshMaskHandler.setImages(obj.refImages.*repmat(permute(((obj.fovTrans./obj.fovPreTrans).^obj.mpOrder),[1,2,4,3]),1,1,obj.slices,1));
            if (~isempty(obj.roiUI))
                for sl=1:obj.slices
                    num = obj.shapeMaskHandler.numRois(sl,realCh);
                    if (num>0)
                        for ind = 1:num
                            obj.calcROIRefPower(sl,realCh,ind);
                        end
                    end
                end
                obj.roiUI.chngSlCh();
            end
        end

        %% functions to manage ROI modulation
        function addROI(obj,sl,ch)
            obj.roiTrans{sl,ch} = [obj.roiTrans{sl,ch},1];
            obj.roiRefPower{sl,ch} = [obj.roiRefPower{sl,ch},1];
            obj.roiThreshPxVal{sl,ch} = [obj.roiThreshPxVal{sl,ch}, 0];
            obj.roiLowPxVal{sl,ch} = [obj.roiLowPxVal{sl,ch}, 0];
            obj.roiHighPxVal{sl,ch} = [obj.roiHighPxVal{sl,ch}, 0];
            if (~isempty(obj.roiUI))
                obj.calcROIRefPower(sl,ch,obj.shapeMaskHandler.numRois(sl,ch));
                obj.roiUI.chngNumROIs();
            end
        end

        function deleteROI(obj,sl,ch,ind)
            temp = ones(1,length(obj.roiTrans{sl,ch}))>0;
            temp(ind) = false;
            obj.roiTrans{sl,ch} = obj.roiTrans{sl,ch}(temp);
            obj.roiRefPower{sl,ch} = obj.roiRefPower{sl,ch}(temp);
            obj.roiThreshPxVal{sl,ch} = obj.roiThreshPxVal{sl,ch}(temp);
            obj.roiLowPxVal{sl,ch} = obj.roiLowPxVal{sl,ch}(temp);
            obj.roiHighPxVal{sl,ch} = obj.roiHighPxVal{sl,ch}(temp);
            if (~isempty(obj.roiUI))
                obj.roiUI.chngNumROIs();
            end
        end

        function calcROIRefPower(obj,sl,ch,ind)
            px = obj.refImages(:,:,sl,ch) ./ (obj.fovPreTrans(:,:,ch) .^ obj.mpOrder);
            px = sort(px(obj.shapeMaskHandler.getRoiMask(sl,ch,ind)));
            threshPx = floor(obj.roiThresh*length(px));
            if (threshPx==0)
                threshPx = 1;
            end
            obj.roiThreshPxVal{sl,ch}(ind) = px(threshPx);
            obj.roiLowPxVal{sl,ch}(ind) = px(1);
            obj.roiHighPxVal{sl,ch}(ind) = px(end);
            obj.roiRefPower{sl,ch}(ind) = sum(px(threshPx:end))-obj.pxVals(1,sl,ch)*(length(px)-threshPx+1);
        end

        function calcAllROIRefPower(obj)
            for sl=1:obj.slices
                for ch=1:(obj.channels*obj.sepChannels+1-obj.sepChannels)
                    num = obj.shapeMaskHandler.numRois(sl,ch);
                    if (num>0)
                        for ind=1:num
                            obj.calcROIRefPower(sl,ch,ind);
                        end
                    end
                end
            end
        end

        function redrawROI(obj,sl,ch,ind)
            if (~isempty(obj.roiUI))
                obj.calcROIRefPower(sl,ch,ind);
                obj.roiUI.chngROI();
            end
        end

        function chngROI(obj)
            if (~isempty(obj.roiUI))
                obj.roiUI.chngROI();
            end
        end
        
        % modulate to get equal signal from rois
        function eqROISig(obj,sl,ch)
            obj.roiTrans{sl,ch} = 1./(obj.roiRefPower{sl,ch}.^(1/obj.mpOrder));
            obj.roiTrans{sl,ch} = obj.roiTrans{sl,ch}./max(obj.roiTrans{sl,ch});
            minVal = min(obj.roiTrans{sl,ch},[],'all');
            if (minVal < obj.minTrans)
                obj.roiTrans{sl,ch} = obj.roiTrans{sl,ch} * obj.minTrans / minVal;
                filter = obj.roiTrans{sl,ch} > 1;
                obj.roiTrans{sl,ch} = filter + (1-filter) .* obj.roiTrans{sl,ch};
            end
        end

        % full transmission to all rois
        function eqROITrans(obj,sl,ch)
            obj.roiTrans{sl,ch} = ones(size(obj.roiTrans{sl,ch}));
        end

        function setROIThresh(obj,val)
            obj.roiThresh = val;
            if (~isempty(obj.roiUI))
                obj.calcAllROIRefPower();
                obj.roiUI.chngROI();
            end
        end

        function setROITrans(obj,sl,ch,ind,val)
            if (val>1)
                obj.roiTrans{sl,ch}(ind) = 1;
            elseif (val*obj.sliceTrans(sl,ch)>obj.minTrans)
                obj.roiTrans{sl,ch}(ind) = val;
            else
                obj.roiTrans{sl,ch}(ind) = obj.minTrans;
            end

            if (~isempty(obj.roiUI))
                obj.roiUI.chngROI();
            end
        end

        %% Output transmission calibration functions
        function cal = getCal(obj,ch)
            cal = obj.calMaps{ch};
        end

        % voltage: column 1
        % transmission: column 2
        % voltage is normalized to range of [0,1]
        % transmission is normalize to max 1 (min may ~= 0)
        % transmission is always ascending down array
        function setCal(obj,ch,cal,calName)
            ch = ch*obj.sepChannels+1-obj.sepChannels;
            if (isempty(cal))
                obj.calNames{ch} = [];
                obj.calMaps{ch} = [];
                obj.absMinTrans = 0;
            else
                cal(:,1) = cal(:,1) - min(cal(:,1));
                cal(:,1) = cal(:,1)/max(cal(:,1));
                [temp,ind] = sort(cal(:,1));
                cal = [temp,cal(ind,2)];
                [~,minInd] = min(cal(:,2));
                [maxT,maxInd] = max(cal(:,2));
                cal = cal(minInd:maxInd,:);
                if (cal(end,2)<cal(1,2))
                    cal = flipud(cal);
                end
                cal(:,2) = cal(:,2)/maxT;
    
                obj.calMaps{ch} = cal;
                obj.calNames{ch} = calName;
                obj.absMinTrans = 0;
                for ch=1:obj.channels
                    if ~isempty(obj.calMaps{ch})
                        temp = min(obj.calMaps{ch}(:,2));
                        if (temp>obj.absMinTrans)
                            obj.absMinTrans = temp;
                        end
                    end
                end

                if (obj.minTrans < obj.absMinTrans)
                    obj.setMinTrans(obj.absMinTrans);
                end
            end
        end

        %% Functions for UI's
        function factor = getImFactor(obj,sl,ch)
            factor = obj.imFactors(:,:,sl,(ch*obj.sepChannels+1-obj.sepChannels));
        end

        function calcImFactor(obj,sl,ch)
            ch=ch*obj.sepChannels+1-obj.sepChannels;
            obj.imFactors(:,:,sl,ch) = (obj.sliceTrans(sl,ch)*obj.fovTrans(:,:,ch));
            minVal = min(obj.imFactors(:,:,sl,ch),[],'all');
            if (minVal < obj.minTrans)
                obj.imFactors(:,:,sl,ch) = obj.imFactors(:,:,sl,ch) * obj.minTrans/minVal;
            end
            filter = obj.imFactors(:,:,sl,ch) > 1;
            obj.imFactors(:,:,sl,ch) = filter+obj.imFactors(:,:,sl,ch).*(1-filter);
            obj.imFactors(:,:,sl,ch) = (obj.imFactors(:,:,sl,ch) ./ (obj.slicePreTrans(sl,ch)*obj.fovPreTrans(:,:,ch))).^obj.mpOrder;
        end

        function calcAllImFactors(obj)
            for sl=1:obj.slices
                for ch=1:(obj.channels*obj.sepChannels+1-obj.sepChannels)
                    obj.calcImFactor(sl,ch);
                end
            end
        end

        function setRoiUI(obj,helper_app)
            obj.roiUI = helper_app;
            if (~isempty(obj.roiUI))
                obj.calcAllROIRefPower();
                obj.roiUI.chngROI();
            end
        end

        function setSliceUI(obj,helper_app)
            obj.sliceUI = helper_app;
        end

        function setFovUI(obj,helper_app)
            obj.fovUI = helper_app;
        end

        %% file IO functions
        function save(obj,fname)
            powerModParams.sliceRefPower = obj.sliceRefPower;
            powerModParams.sliceTrans = obj.sliceTrans;
            powerModParams.slicePreTrans = obj.slicePreTrans;
            powerModParams.slicePreTransNames = obj.slicePreTransNames;
            powerModParams.z = obj.z;
            powerModParams.sliceHighPx = obj.sliceHighPx;
            powerModParams.sliceLowPx = obj.sliceLowPx;
            powerModParams.sliceSep = obj.sliceSep;
            powerModParams.attenLen = obj.attenLen;
            powerModParams.orientation = obj.orientation;
            powerModParams.curveTypes = obj.curveTypes;

            powerModParams.fovTrans = obj.fovTrans;
            powerModParams.fovPreTrans = obj.fovPreTrans;
            powerModParams.fovTransNames = obj.fovTransNames;
            powerModParams.fovPreTransNames = obj.fovPreTransNames;
            
            powerModParams.roiTrans = obj.roiTrans;
            powerModParams.roiThresh = obj.roiThresh;

            powerModParams.minTrans = obj.minTrans;
            powerModParams.mpOrder = obj.mpOrder;
            powerModParams.calMaps = obj.calMaps;
            powerModParams.calNames = obj.calNames;

            save(fname,"powerModParams",'-append');
        end

        function success = load(obj,data)
            success = false;
            if (~isfield(data,'powerModParams'))
                return;
            end

            realCh = obj.channels*obj.sepChannels+1-obj.sepChannels;
            params = data.powerModParams;
            if (isfield(params,'refPower'))
                obj.sliceRefPower = params.refPower';
                obj.sliceTrans = params.trans';
                obj.slicePreTrans = params.prePower';
                obj.sliceHighPx = params.highPx;
                obj.sliceLowPx = params.lowPx;
                obj.slicePreTransNames = params.prePowerNames;

                obj.fovTrans = ones(obj.imDim(1),obj.imDim(2),realCh);
                obj.fovTransNames = cell(realCh,1);
                obj.fovPreTrans = ones(obj.imDim(1),obj.imDim(2),realCh);
                obj.fovPreTransNames = cell(realCh,1);

                obj.roiThresh = 0.95;
                obj.roiTrans = cell(obj.slices,realCh);
                for sl=1:obj.slices
                    for ch=1:realCh
                        num = obj.shapeMaskHandler.numRois(sl,ch);
                        if (num>0)
                            obj.roiTrans{sl,ch} = ones(1,num);
                        end
                    end
                end
                obj.minTrans = 0;
            elseif (isfield(params,'sliceTrans'))
                obj.sliceRefPower = params.sliceRefPower;
                obj.sliceTrans = params.sliceTrans;
                obj.slicePreTrans = params.slicePreTrans;
                obj.slicePreTransNames = params.slicePreTransNames;
                obj.sliceHighPx = params.sliceHighPx;
                obj.sliceLowPx = params.sliceLowPx;
                
                obj.fovTrans = params.fovTrans;
                obj.fovTransNames = params.fovTransNames;
                obj.fovPreTrans = params.fovPreTrans;
                obj.fovPreTransNames = params.fovPreTransNames;

                obj.roiTrans = params.roiTrans;
                obj.roiThresh = params.roiThresh;

                obj.minTrans = params.minTrans;
            else
                return;
            end
            obj.calNames = params.calNames;
            obj.calMaps = params.calMaps;
            obj.z = params.z;
            obj.sliceSep = params.sliceSep;
            obj.attenLen = params.attenLen;
            obj.mpOrder = params.mpOrder;
            obj.orientation = params.orientation;
            obj.curveTypes = params.curveTypes;

            obj.roiRefPower = cell(obj.slices,realCh);
            obj.roiHighPxVal = cell(obj.slices,realCh);
            obj.roiLowPxVal = cell(obj.slices,realCh);
            obj.roiThreshPxVal = cell(obj.slices,realCh);
            for sl=1:obj.slices
                for ch=1:realCh
                    num=obj.shapeMaskHandler.numRois(sl,ch);
                    if (num>0)
                        obj.roiRefPower{sl,ch} = zeros(1,num);
                        obj.roiHighPxVal{sl,ch} = zeros(1,num);
                        obj.roiLowPxVal{sl,ch} = zeros(1,num);
                        obj.roiThreshPxVal{sl,ch} = zeros(1,num);
                    end
                end
            end

            obj.absMinTrans = 0;
            for ch=1:obj.channels
                if (~isempty(obj.calMaps{ch}))
                    temp = min(obj.calMaps{ch}(:,2),[],'all');
                    if (temp>obj.absMinTrans)
                        obj.absMinTrans = temp;
                    end
                end
            end

            if (obj.minTrans < obj.absMinTrans)
                obj.minTrans = obj.absMinTrans;
            end
            obj.calcAllImFactors();
            obj.threshMaskHandler.setImages(obj.refImages.*repmat(permute(((obj.fovTrans./obj.fovPreTrans).^obj.mpOrder),[1,2,4,3]),1,1,obj.slices,1));
            success = true;
        end

        function valid = checkValid(obj,powerMod)
            valid = false;
            if ((powerMod.channels ~= obj.channels) || (powerMod.slices ~= obj.slices) || (obj.sepChannels ~= powerMod.sepChannels))
                return;
            end

            obj.numPx = powerMod.numPx;
            obj.mainChannel = powerMod.mainChannel;
            valid = true;
        end

        %% Functions for building final AES mask
        function mask = buildAESMask(obj,sl,ch,modShapes)
            realCh = ch*obj.sepChannels+1-obj.sepChannels;
            mask = obj.threshMaskHandler.getFullMask(sl,realCh) .* obj.fovTrans(:,:,realCh) * obj.sliceTrans(sl,realCh);
            mask = mask .* (1 - obj.shapeMaskHandler.getChMask(sl,realCh));
            num = obj.shapeMaskHandler.numRois(sl,realCh);
            if (num>0)
                tempROIs = zeros(obj.imDim(1),obj.imDim(2),num);
                compMask = zeros(obj.imDim(1),obj.imDim(2));
                for ind=1:num
                    tempROIs(:,:,ind) = obj.shapeMaskHandler.getRoiMask(sl,realCh,ind) .* (1 - compMask);
                    compMask = compMask + tempROIs(:,:,ind);
                end

                if (modShapes)
                    powROIs = zeros(obj.imDim(1),obj.imDim(2),num);
                    for ind=1:num
                        powROIs(:,:,ind) = tempROIs(:,:,ind) .* obj.fovTrans(:,:,realCh);
                        powROIs(:,:,ind) = powROIs(:,:,ind) * obj.roiTrans{sl,realCh}(ind) * sum(tempROIs(:,:,ind),'all') / sum(powROIs(:,:,ind),'all');
                    end
                    powROIs = powROIs * obj.sliceTrans(sl,realCh) / max(powROIs,[],'all');
                    for ind=1:num
                        mask = mask + powROIs(:,:,ind);
                    end
                else
                    for ind=1:num
                        mask = mask + obj.sliceTrans(sl,realCh) * tempROIs(:,:,ind) .* obj.roiTrans{sl,realCh}(ind);
                    end
                end 
            end

            minVal = min(mask(mask>0),[],'all');
            if (minVal < obj.minTrans)
                mask = mask * obj.minTrans / minVal;
            end
            filter = (mask>1);
            mask = filter + (1-filter) .* mask;

            if (~isempty(obj.calMaps{ch}))
                vals = unique(mask);
                uniform =  (length(vals) < sqrt(obj.imDim(1)*obj.imDim(2)));
                if (uniform)
                    tMask = mask;
                    mask = zeros(size(mask));
                    for ii=1:length(vals)
                        if (vals(ii)==0)
                            continue;
                        end
                        
                        ind = 1;
                        while(obj.calMaps{ch}(ind,2) < vals(ii))
                            ind = ind + 1;
                        end

                        if (ind > 1)
                            v = obj.calMaps{ch}(ind,1) + (vals(ii)-obj.calMaps{ch}(ind,2)) * (obj.calMaps{ch}(ind,1)-obj.calMaps{ch}((ind-1),1)) / (obj.calMaps{ch}(ind,2)-obj.calMaps{ch}((ind-1),2));
                        else
                            v = obj.calMaps{ch}(1,1);
                        end
                        mask = mask + v*(tMask==vals(ii));
                    end
                else
                    for y=1:obj.imDim(1)
                        for x = 1:obj.imDim(2)
                            if (mask(y,x)==0)
                                continue;
                            end

                            ind = 1;
                            while (obj.calMaps{ch}(ind,2) < mask(y,x))
                                ind = ind + 1;
                            end
    
                            if (ind>1)
                                mask(y,x) = obj.calMaps{ch}(ind,1) + (mask(y,x)-obj.calMaps{ch}(ind,2)) * (obj.calMaps{ch}(ind,1)-obj.calMaps{ch}((ind-1),1)) / (obj.calMaps{ch}(ind,2)-obj.calMaps{ch}((ind-1),2));
                            else
                                mask(y,x) = obj.calMaps{ch}(1,1);
                            end
                        end
                    end
                end
            end
        end

        %% Other Functions
        function setMinTrans(obj,val)
            if (val > obj.absMinTrans)
                obj.minTrans = val;
            else
                obj.minTrans = obj.absMinTrans;
            end

            obj.calcAllSliceTrans();
            filter = obj.sliceTrans > obj.minTrans;
            obj.sliceTrans = obj.sliceTrans.*filter + (1-filter)*obj.minTrans;
            obj.calcAllImFactors();
            for ch=1:(obj.channels*obj.sepChannels+1-obj.sepChannels)
                for sl=1:obj.slices
                    minVal = min(obj.roiTrans{sl,ch});
                    if (minVal < obj.minTrans)
                        obj.roiTrans{sl,ch} = obj.roiTrans{sl,ch} * obj.minTrans / minVal;
                        filter = obj.roiTrans{sl,ch}>1;
                        obj.roiTrans{sl,ch} = filter + (1-filter) .* obj.roiTrans{sl,ch};
                    end
                end
            end

            if (~isempty(obj.roiUI))
                obj.roiUI.chngSlCh();
            end

            if (~isempty(obj.sliceUI))
                obj.sliceUI.updatePlots();
            end
        end

        function setMpOrder(obj,val)
            obj.mpOrder = val;
            obj.calcAllSliceTrans();

            if (~isempty(obj.roiUI))
                obj.roiUI.updateMpOrder();
            end

            if (~isempty(obj.sliceUI))
                obj.sliceUI.updateMpOrder();
            end

            if (~isempty(obj.fovUI))
                obj.sliceUI.updateMpOrder();
            end
        end   
    end
end