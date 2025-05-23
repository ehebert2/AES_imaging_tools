%%%%%%%%%%%%%%%%%
% Object for managing masks used for segmentation and motion correction
%%%%%%%%%%%%%%%%%

classdef MaskHandler < handle
    properties (Access = private)
        expMask
        expMaskPrev
        expMaskRef
        expMaskProp
        expMaskPropPrev

        roiSmpl
        roiSmplOutline
        roiSmplPrev
        roiSmplOutlinePrev
        roiSmplRef
        roiSmplProp
        roiSmplPropPrev
        roiSmplNames
        numRoiSmpl

        roiAes
        roiAesOutline
        roiAesPrev
        roiAesOutlinePrev
        roiAesRef
        roiAesProp
        roiAesPropPrev
        roiAesNames
        numRoiAes

        outlineKernel
    end

    properties (GetAccess = public)
        dim
        splitCh
        channels
        channel
        slices
        slice
        isVolume

        roiAesImage
        roiSmplImage
        roiSlctImage
        roiAesSlct
        roiSmplSlct

        showOutlineSmpl
        showOutlineAes
    end

    methods
        function obj = MaskHandler(params)
            obj.splitCh = params.splitChannels;
            obj.dim = params.dim;
            obj.channel = params.mainChannel;
            obj.channels = params.channels;
            obj.isVolume = params.volume;
            if (obj.isVolume)
                obj.slices = params.slices;
            else
                obj.slices = 1;
            end
            obj.slice = 1;
            obj.roiSmplNames = cell(obj.slices,obj.channels);
            obj.roiAesNames = cell(obj.slices,obj.channels);
            obj.roiSmpl = cell(obj.slices,obj.channels);
            obj.roiSmplRef = cell(obj.slices,obj.channels);
            obj.roiSmplOutline = cell(obj.slices,obj.channels);
            obj.roiAes = cell(obj.slices,obj.channels);
            obj.roiAesRef = cell(obj.slices,obj.channels);
            obj.roiAesOutline = cell(obj.slices,obj.channels);
            obj.expMask = zeros(obj.dim(1),obj.dim(2),obj.slices,obj.channels);
            obj.expMaskRef = zeros(obj.dim(1),obj.dim(2),obj.slices,obj.channels);
            obj.roiSmplProp = cell(obj.slices,obj.channels);
            obj.roiAesProp = cell(obj.slices,obj.channels);
            obj.expMaskProp = zeros(obj.slices,obj.channels,3);
            for sl=1:obj.slices
                for ch=1:obj.channels
                    obj.roiSmplNames{sl,ch} = {};
                    obj.roiAesNames{sl,ch} = {};
                    obj.roiSmpl{sl,ch} = ([]>0);
                    obj.roiSmplOutline{sl,ch} = ([]>0);
                    obj.roiSmplRef{sl,ch} = [];
                    obj.roiAes{sl,ch} = ([]>0);
                    obj.roiAesOutline{sl,ch} = ([]>0);
                    obj.roiAesRef{sl,ch} = [];
                    obj.roiSmplProp{sl,ch} = [];
                    obj.roiAesProp{sl,ch} = [];
                end
            end
            obj.numRoiSmpl = zeros(obj.slices,obj.channels);
            obj.numRoiAes = zeros(obj.slices,obj.channels);
            obj.showOutlineSmpl = false;
            obj.showOutlineAes = false;
            tempR = 2;
            tempX = (-tempR):tempR;
            tempX = repmat(tempX,(2*tempR+1),1);
            tempY = (-tempR):tempR;
            tempY = repmat(tempY',1,(2*tempR+1));
            obj.outlineKernel = ((tempX.^2+tempY.^2) <= (tempR^2));
            obj.outlineKernel = obj.outlineKernel * 1.01 / sum(obj.outlineKernel,'all');
        end

        %% functions for adding/removing/renaming masks
        function addRoiSmpl(obj, name, mask)
            obj.roiSmplRef{obj.slice,obj.channel} = cat(3,obj.roiSmplRef{obj.slice,obj.channel},mask);
            obj.roiSmpl{obj.slice,obj.channel} = cat(3,obj.roiSmpl{obj.slice,obj.channel},mask);
            obj.roiSmplOutline{obj.slice,obj.channel} = cat(3,obj.roiSmplOutline{obj.slice,obj.channel},obj.buildOutline(mask));
            obj.roiSmplNames{obj.slice,obj.channel} = [obj.roiSmplNames{obj.slice,obj.channel},name];
            obj.roiSmplProp{obj.slice,obj.channel} = [obj.roiSmplProp{obj.slice,obj.channel};[0,0,0]];
            obj.numRoiSmpl(obj.slice,obj.channel) = obj.numRoiSmpl(obj.slice,obj.channel) + 1;
            obj.roiSmplSlct = [obj.roiSmplSlct, 0];
        end

        function addRoiAes(obj, name, mask)
            obj.roiAesRef{obj.slice,obj.channel} = cat(3,obj.roiAesRef{obj.slice,obj.channel},mask);
            obj.roiAes{obj.slice,obj.channel} = cat(3,obj.roiAes{obj.slice,obj.channel},mask);
            obj.roiAesOutline{obj.slice,obj.channel} = cat(3,obj.roiAesOutline{obj.slice,obj.channel},obj.buildOutline(mask));
            obj.roiAesNames{obj.slice,obj.channel} = [obj.roiAesNames{obj.slice,obj.channel},name];
            obj.roiAesProp{obj.slice,obj.channel} = [obj.roiAesProp{obj.slice,obj.channel};[0,0,0]];
            obj.numRoiAes(obj.slice,obj.channel) = obj.numRoiAes(obj.slice,obj.channel) + 1;
            obj.roiAesSlct = [obj.roiAesSlct, 0];
        end

        function removeRoiSmpl(obj,toRemove)
            obj.roiSmplRef{obj.slice,obj.channel} = obj.roiSmplRef{obj.slice,obj.channel}(:,:,~toRemove);
            obj.roiSmpl{obj.slice,obj.channel} = obj.roiSmpl{obj.slice,obj.channel}(:,:,~toRemove);
            obj.roiSmplOutline{obj.slice,obj.channel} = obj.roiSmplOutline{obj.slice,obj.channel}(:,:,~toRemove);
            obj.roiSmplNames{obj.slice,obj.channel} = obj.roiSmplNames{obj.slice,obj.channel}(~toRemove);
            obj.roiSmplProp{obj.slice,obj.channel} = obj.roiSmplProp{obj.slice,obj.channel}(~toRemove,:);
            obj.numRoiSmpl(obj.slice,obj.channel) = obj.numRoiSmpl(obj.slice,obj.channel) - sum(toRemove);
            obj.roiSmplSlct = obj.roiSmplSlct(~toRemove);
        end

        function removeRoiAes(obj,toRemove)
            obj.roiAesRef{obj.slice,obj.channel} = obj.roiAesRef{obj.slice,obj.channel}(:,:,~toRemove);
            obj.roiAes{obj.slice,obj.channel} = obj.roiAes{obj.slice,obj.channel}(:,:,~toRemove);
            obj.roiAesOutline{obj.slice,obj.channel} = obj.roiAesOutline{obj.slice,obj.channel}(:,:,~toRemove);
            obj.roiAesNames{obj.slice,obj.channel} = obj.roiAesNames{obj.slice,obj.channel}(~toRemove);
            obj.roiAesProp{obj.slice,obj.channel} = obj.roiAesProp{obj.slice,obj.channel}(~toRemove,:);
            obj.numRoiAes(obj.slice,obj.channel) = obj.numRoiAes(obj.slice,obj.channel) - sum(toRemove);
            obj.roiAesSlct = obj.roiAesSlct(~toRemove);
        end

        function removeExpMask(obj)
            obj.expMask(:,:,obj.slice,obj.channel) = obj.expMask(:,:,obj.slice,obj.channel) * 0;
        end

        function renameRoiSmpl(obj,ind,name)
            obj.roiSmplNames{obj.slice,obj.channel}(ind) = name;
        end

        function renameRoiAes(obj,ind,name)
            obj.roiAesNames{obj.slice,obj.channel}(ind) = name;
        end

        function setExpMask(obj,mask)
            obj.expMask(:,:,obj.slice,obj.channel) = mask;
            obj.expMaskRef(:,:,obj.slice,obj.channel) = mask;
            obj.expMaskProp(obj.slice,obj.channel,:) = [0 0 0];
        end

        %% function for interfacing ui display in app
        function selectRoiAes(obj,indices)
            obj.roiAesSlct = 0 * obj.roiAesSlct;
            if (~isempty(indices))
                for ii=1:length(indices)
                    obj.roiAesSlct(indices{ii}) = true;
                end
            end
        end

        function selectRoiSmpl(obj,indices)
            obj.roiSmplSlct = 0 * obj.roiSmplSlct;
            if (~isempty(indices))
                for ii=1:length(indices)
                    obj.roiSmplSlct(indices{ii}) = true;
                end
            end
        end

        % called when user clicks on image to determin if mask is selected
        function selectCoord(obj,coord,cntl)
            if (obj.numRoiSmpl(obj.slice,obj.channel)>0)
                if (cntl)
                    obj.roiSmplSlct = (obj.roiSmplSlct + squeeze(obj.roiSmpl{obj.slice,obj.channel}(coord(2),coord(1),:))) > 0;
                else
                    obj.roiSmplSlct = squeeze(obj.roiSmpl{obj.slice,obj.channel}(coord(2),coord(1),:));
                end

                if (obj.numRoiAes(obj.slice,obj.channel)>0)
                    if (sum(obj.roiSmplSlct)==0)
                        if (cntl)
                            obj.roiAesSlct = (obj.roiAesSlct + squeeze(obj.roiAes{obj.slice,obj.channel}(coord(2),coord(1),:))) > 0;
                        else
                            obj.roiAesSlct = squeeze(obj.roiAes{obj.slice,obj.channel}(coord(2),coord(1),:));
                        end
                    else
                        obj.roiAesSlct = 0 * obj.roiAesSlct;
                    end
                end
            elseif(obj.numRoiAes(obj.slice,obj.channel)>0)
                if (cntl)
                    obj.roiAesSlct = (obj.roiAesSlct + squeeze(obj.roiAes{obj.slice,obj.channel}(coord(2),coord(1),:))) > 0;
                else
                    obj.roiAesSlct = squeeze(obj.roiAes{obj.slice,obj.channel}(coord(2),coord(1),:));
                end
            end
        end

        % called by main app to build image overlay highlighting masks
        function buildRoiImage(obj)
            if (isempty(obj.roiAes{obj.slice,obj.channel}))
                obj.roiAesImage = zeros(obj.dim);
                aesSlctImage = zeros(obj.dim);
            else
                obj.roiAesSlct = obj.roiAesSlct > 0;
                if (obj.showOutlineAes)
                    obj.roiAesImage = sum(obj.roiAesOutline{obj.slice,obj.channel},3);
                    aesSlctImage = sum(obj.roiAesOutline{obj.slice,obj.channel}(:,:,obj.roiAesSlct),3);
                else
                    obj.roiAesImage = sum(obj.roiAes{obj.slice,obj.channel},3);
                    aesSlctImage = sum(obj.roiAes{obj.slice,obj.channel}(:,:,obj.roiAesSlct),3);
                end
                obj.roiAesImage = obj.roiAesImage > 0;
            end

            if (isempty(obj.roiSmpl{obj.slice,obj.channel}))
                obj.roiSmplImage = zeros(obj.dim);
                smplSlctImage = zeros(obj.dim);
            else
                obj.roiSmplSlct = obj.roiSmplSlct > 0;
                if (obj.showOutlineSmpl)
                    obj.roiSmplImage = sum(obj.roiSmplOutline{obj.slice,obj.channel},3);
                    smplSlctImage = sum(obj.roiSmplOutline{obj.slice,obj.channel}(:,:,obj.roiSmplSlct),3);
                else
                    obj.roiSmplImage = sum(obj.roiSmpl{obj.slice,obj.channel},3);
                    smplSlctImage = sum(obj.roiSmpl{obj.slice,obj.channel}(:,:,obj.roiSmplSlct),3);
                end
                obj.roiSmplImage = obj.roiSmplImage > 0;
            end
            obj.roiSlctImage = (aesSlctImage + smplSlctImage) > 0;
        end

        function maskOutline = buildOutline(obj,mask)
            maskOutline = conv2(mask,obj.outlineKernel,'same');
            maskOutline = maskOutline >= 1;
            maskOutline = mask - maskOutline;
        end

        %% functions for adjusting mask position/sizes
        % hold current state incase you want to revert after making changes
        function storeState(obj)
            obj.roiSmplPropPrev = obj.roiSmplProp;
            obj.roiAesPropPrev = obj.roiAesProp;
            obj.expMaskPropPrev = obj.expMaskProp;
            obj.roiSmplPrev = obj.roiSmpl;
            obj.roiSmplOutlinePrev = obj.roiSmplOutline;
            obj.roiAesPrev = obj.roiAes;
            obj.roiAesOutlinePrev = obj.roiAesOutline;
            obj.expMaskPrev = obj.expMask;
        end

        % revert masks to previous stored state
        function revert(obj)
            obj.roiSmplProp = obj.roiSmplPropPrev;
            obj.roiAesProp = obj.roiAesPropPrev;
            obj.expMaskProp = obj.expMaskPropPrev;
            obj.roiSmpl = obj.roiSmplPrev;
            obj.roiSmplOutline = obj.roiSmplOutlinePrev;
            obj.roiAes = obj.roiAesPrev;
            obj.roiAesOutline = obj.roiAesOutlinePrev;
            obj.expMask = obj.expMaskPrev;
        end

        % called by main app to redraw masks
        function setMaskProp(obj, x, y, mgn, index, type)
            switch type
                case MaskType.EXPMASK
                    obj.expMaskProp(obj.slice,obj.channel,:) = [x,y,mgn];
                    obj.expMask(:,:,obj.slice,obj.channel) = obj.redrawMask(obj.expMaskRef(:,:,obj.slice,obj.channel),x,y,mgn);
                case MaskType.AESALL
                    obj.roiAesProp{obj.slice,obj.channel} = repmat([x,y,mgn],obj.numRoiAes(obj.slice,obj.channel),1);
                    for ii=1:obj.numRoiAes(obj.slice,obj.channel)
                        obj.roiAes{obj.slice,obj.channel}(:,:,ii) = obj.redrawMask(obj.roiAesRef{obj.slice,obj.channel}(:,:,ii),x,y,mgn);
                        obj.roiAesOutline{obj.slice,obj.channel}(:,:,ii) = obj.buildOutline(obj.roiAes{obj.slice,obj.channel}(:,:,ii));
                    end
                case MaskType.AES
                    obj.roiAesProp{obj.slice,obj.channel}(index,:) = [x,y,mgn];
                    obj.roiAes{obj.slice,obj.channel}(:,:,index) = obj.redrawMask(obj.roiAesRef{obj.slice,obj.channel}(:,:,index),x,y,mgn);
                    obj.roiAesOutline{obj.slice,obj.channel}(:,:,index) = obj.buildOutline(obj.roiAes{obj.slice,obj.channel}(:,:,index));
                case MaskType.SMPLALL
                    obj.roiSmplProp{obj.slice,obj.channel} = repmat([x,y,mgn],obj.numRoiSmpl(obj.slice,obj.channel),1);
                    for ii=1:obj.numRoiSmpl(obj.slice,obj.channel)
                        obj.roiSmpl{obj.slice,obj.channel}(:,:,ii) = obj.redrawMask(obj.roiSmplRef{obj.slice,obj.channel}(:,:,ii),x,y,mgn);
                        obj.roiSmplOutline{obj.slice,obj.channel}(:,:,ii) = obj.buildOutline(obj.roiSmpl{obj.slice,obj.channel}(:,:,ii));
                    end
                case MaskType.SMPL
                    obj.roiSmplProp{obj.slice,obj.channel}(index,:) = [x,y,mgn];
                    obj.roiSmpl{obj.slice,obj.channel}(:,:,index) = obj.redrawMask(obj.roiSmplRef{obj.slice,obj.channel}(:,:,index),x,y,mgn);
                    obj.roiSmplOutline{obj.slice,obj.channel}(:,:,index) = obj.buildOutline(obj.roiSmpl{obj.slice,obj.channel}(:,:,index));
            end
        end

        function flattenAesProp(obj)
            if (obj.numRoiAes(obj.slice,obj.channel)>1)
                prop = max(obj.roiAesProp{obj.slice,obj.channel});
            else
                prop = obj.roiAesProp{obj.slice,obj.channel};
            end
            obj.roiAesProp{obj.slice,obj.channel} = repmat(prop,obj.numRoiAes(obj.slice,obj.channel),1);
            obj.roiAesSlct = 0*obj.roiAesSlct+1;
            for ii=1:obj.numRoiAes(obj.slice,obj.channel)
                obj.roiAes{obj.slice,obj.channel}(:,:,ii) = obj.redrawMask(obj.roiAesRef{obj.slice,obj.channel}(:,:,ii),prop(1),prop(2),prop(3));
                obj.roiAesOutline{obj.slice,obj.channel}(:,:,ii) = obj.buildOutline(obj.roiAes{obj.slice,obj.channel}(:,:,ii));
            end
        end

        function flattenSmplProp(obj)
            if (obj.numRoiSmpl(obj.slice,obj.channel)>1)
                prop = max(obj.roiSmplProp{obj.slice,obj.channel});
            else
                prop = obj.roiSmplProp{obj.slice,obj.channel};
            end
            obj.roiSmplProp{obj.slice,obj.channel} = repmat(prop,obj.numRoiSmpl(obj.slice,obj.channel),1);
            obj.roiSmplSlct = 0*obj.roiSmplSlct+1;
            for ii=1:obj.numRoiSmpl(obj.slice,obj.channel)
                obj.roiSmpl{obj.slice,obj.channel}(:,:,ii) = obj.redrawMask(obj.roiSmplRef{obj.slice,obj.channel}(:,:,ii),prop(1),prop(2),prop(3));
                obj.roiSmplOutline{obj.slice,obj.channel}(:,:,ii) = obj.buildOutline(obj.roiSmpl{obj.slice,obj.channel}(:,:,ii));
            end
        end

        %% Getters and Setters
        function present = hasExpMask(obj)
            present = sum(obj.expMask(:,:,obj.slice,obj.channel),'all')>0;
        end

        function rois = getRoiAes(obj)
            rois = obj.roiAes{obj.slice,obj.channel};
        end

        function rois = getRoiSmpl(obj)
            rois = obj.roiSmpl{obj.slice,obj.channel};
        end

        function names = getAesNames(obj)
            names = obj.roiAesNames{obj.slice,obj.channel};
        end

        function names = getSmplNames(obj)
            names = obj.roiSmplNames{obj.slice,obj.channel};
        end

        function prop = getAesProp(obj,ind)
            prop = obj.roiAesProp{obj.slice,obj.channel}(ind,:);
        end

        function prop = getSmplProp(obj,ind)
            prop = obj.roiSmplProp{obj.slice,obj.channel}(ind,:);
        end

        function prop = getExpProp(obj)
            prop = obj.expMaskProp(obj.slice,obj.channel,:);
        end

        function mask = getExpMask(obj)
            mask = obj.expMask(:,:,obj.slice,obj.channel);
        end

        function num = getNumAes(obj)
            num = obj.numRoiAes(obj.slice,obj.channel);
        end

        function num = getNumSmpl(obj)
            num = obj.numRoiSmpl(obj.slice,obj.channel);
        end 

        function setOutlineAes(obj,enable)
            obj.showOutlineAes = enable;
            obj.buildRoiImage();
        end

        function setOutlineSmpl(obj,enable)
            obj.showOutlineSmpl = enable;
            obj.buildRoiImage();
        end
        
        %% Other Functions
        % Change color channel
        function setChannel(obj,slice,channel)
            if (obj.splitCh)
                obj.channel = channel;
            end
            obj.slice = slice;
            obj.roiAesSlct = zeros(obj.numRoiAes(obj.slice,obj.channel),1);
            obj.roiSmplSlct = zeros(obj.numRoiSmpl(obj.slice,obj.channel),1);  
        end

        % attach masks to params struct for processing files
        function params = attachParams(obj,params,reslicer)
            for sl = 1:obj.slices
                for ch = 1:obj.channels
                    obj.roiAes{sl,ch} = (obj.roiAes{sl,ch} > 0);
                    obj.roiSmpl{sl,ch} = (obj.roiSmpl{sl,ch} > 0);
                end
            end
            
            if (obj.splitCh)
                fullmask = zeros(size(obj.expMask));
                tempCh = obj.channel;
                tempSl = obj.slice;
                for sl=1:obj.slices
                    for ch = 1:obj.channels
                        obj.setChannel(sl,ch);
                        obj.buildRoiImage();
                        fullmask(:,:,sl,ch) = (obj.roiAesImage+obj.roiSmplImage+obj.roiSlctImage+obj.expMask(:,:,sl,ch))>0;
                    end
                end
                obj.setChannel(tempSl,tempCh);
            else
                fullmask = zeros(size(obj.expMask));
                tempSl = obj.slice;
                for sl=1:obj.slices
                    obj.setChannel(sl,1);
                    obj.buildRoiImage();
                    fullmask(:,:,sl,1) = (obj.roiAesImage+obj.roiSmplImage+obj.roiSlctImage+obj.expMask(:,:,sl,1))>0;
                    for ch=2:obj.channels
                        obj.roiAes{sl,ch} = obj.roiAes{sl,1};
                        obj.roiSmpl{sl,ch} = obj.roiSmpl{sl,1};
                        fullmask(:,:,sl,ch) = fullmask(:,:,sl,1);
                        obj.roiAesNames{sl,ch} = obj.roiAesNames{sl,1};
                        obj.roiSmplNames{sl,ch} = obj.roiSmplNames{sl,1};
                    end      
                end
                obj.setChannel(tempSl,1);
            end
            
            if (params.reslice && (~isempty(reslicer)))
                params.roiAes = reslicer.mapCells(obj.roiAes);
                params.roiSmpl = reslicer.mapCells(obj.roiSmpl);
                params.aesNames = reslicer.mapCells(obj.roiAesNames);
                params.smplNames = reslicer.mapCells(obj.roiSmplNames);
                params.expMask = (reslicer.mapImages(fullmask)>0);
                params.numRoiAes = reslicer.mapMatrix(obj.numRoiAes);
                params.numRoiSmpl = reslicer.mapMatrix(obj.numRoiSmpl);
            else
                params.roiAes = obj.roiAes;
                params.roiSmpl = obj.roiSmpl;
                params.aesNames = obj.roiAesNames;
                params.smplNames = obj.roiSmplNames;
                params.expMask = fullmask;
                params.numRoiAes = obj.numRoiAes;
                params.numRoiSmpl = obj.numRoiSmpl;
            end
            params.splitChannels = obj.splitCh;
            params.dim = obj.dim;
        end

        % parses masks from directory and adds them
        function loadMaskFolder(obj,inPath,type)
            tempChannels = 1;
            while (exist(fullfile(inPath,strcat('channel_',num2str(tempChannels))),'dir'))
                tempChannels=tempChannels+1;
            end
            tempChannels=tempChannels-1;
            if (tempChannels==0)
                tempDir = inPath;
            else
                tempDir = fullfile(inPath,'channel_1');
            end

            tempSlices = 1;
            while (exist(fullfile(tempDir,strcat('slice_',num2str(tempSlices))),'dir'))
                tempSlices=tempSlices+1;
            end
            tempSlices = tempSlices-1;
            
            if (tempChannels>0)
                if (tempSlices>0)
                    targetCh = min((obj.channels*obj.splitCh+(1-obj.splitCh)),tempChannels);
                    targetSl = min(obj.slices,tempSlices);
                    currentSlice = obj.slice;
                    currentChannel = obj.channel;
                    for sl=1:targetSl
                        for ch=1:targetCh
                            obj.setChannel(sl,ch);
                            tempPath = strcat(inPath,filesep,'channel_',num2str(ch),filesep,'slice_',num2str(sl));
                            obj.readMasks(tempPath,type);
                        end
                    end
                    obj.setChannel(currentSlice,currentChannel);
                else
                    targetCh = min((obj.channels*obj.splitCh+(1-obj.splitCh)),tempChannels);
                    currentChannel = obj.channel;
                    for ch=1:targetCh
                        obj.setChannel(obj.slice,ch);
                        obj.readMasks(strcat(inPath,filesep,'channel_',num2str(ch)),type);
                    end
                    obj.setChannel(obj.slice,currentChannel);
                end
            else
                if (tempSlices>0)
                    targetSl = min(obj.slices,tempSlices);
                    currentSlice = obj.slice;
                    for sl=1:targetSl
                        obj.setChannel(sl,obj.channel);
                        obj.readMasks(strcat(inPath,filesep,'slice_',num2str(sl)),type);
                    end
                    obj.setChannel(currentSlice,obj.channel);
                else
                    obj.readMasks(inPath,type);
                end
            end
        end
    end

    methods (Access = private)
        % helper for reading all masks in a folder
        function  readMasks(obj,path,type)
            temp = dir(fullfile(path,'*.tif'));
            names = cell(length(temp),1);
            masks = zeros(obj.dim(1),obj.dim(2),length(temp));
            keep = (zeros(length(temp),1)>0);
            for ii = 1:length(temp)
                temp2 = split(temp(ii).name,'.');
                names{ii} = temp2{1};
                tempTiff = double(imread(fullfile(path,strcat(names{ii},'.tif')),1));
                if (prod(size(tempTiff)==obj.dim))
                    masks(:,:,ii) = (tempTiff>0);
                    keep(ii) = true;
                end
            end

            if (sum(keep)>0)
                names = names(keep);
                masks = masks(:,:,keep);
                if (type==MaskType.AES)
                    for ii=1:length(names)
                        if(~any(strcmp(names{ii},obj.roiAesNames{obj.slice,obj.channel})))
                            obj.addRoiAes(names{ii},masks(:,:,ii));
                        end
                    end
                elseif (type==MaskType.SMPL)
                    for ii=1:length(names)
                        if (~any(strcmp(names{ii},obj.roiSmplNames{obj.slice,obj.channel})))
                            obj.addRoiSmpl(names{ii},masks(:,:,ii));
                        end
                    end
                end
            end
        end
    end

    methods (Static)
        % used to adjust mask position/dialation
        function adjMask = redrawMask(mask, x, y, mgn)
            adjMask = circshift(mask,[y,x]);
            if(x > 0)
                adjMask(:,1:x) = zeros(size(adjMask,1),x);
            elseif (x < 0)
                adjMask(:,(end+x+1):end) = zeros(size(adjMask,1),-x);
            end

            if(y > 0)
                adjMask(1:y,:) = zeros(y,size(adjMask,2));
            elseif (y < 0)
                adjMask((end+y+1):end,:) = zeros(-y,size(adjMask,2));
            end

            if (mgn>0)
                kernel = ones(mgn+1,mgn+1);
                adjMask = conv2(adjMask,kernel,'same');
                adjMask = adjMask > 0;
            elseif (mgn<0)
                w=1-mgn;
                kernel = ones(w,w)*1.01/w^2;
                adjMask = conv2(adjMask,kernel,'same');
                adjMask = adjMask > 1;
            end
        end

        % save masks to file 
        % (uses masks in params, not the ones in the object)
        function msg = saveMasks(params,outPath)
            [status, msg] = mkdir(outPath,'masks');
            if (~status)
                return;
            end
            maskPath = fullfile(outPath,'masks');

            [status, msg] = mkdir(maskPath,'AES');
            if (~status)
                return;
            end
            aesPath = fullfile(maskPath,'AES');
                
            [status, msg] = mkdir(maskPath,'sample');
            if (~status)
                return;
            end
            smplPath = fullfile(maskPath,'sample');

            tempSlices = size(params.roiAes,1);
            tempChannels = size(params.roiAes,2)*params.splitChannels + (1-params.splitChannels);
            for ch=1:tempChannels
                if (params.splitChannels)
                    chFolder = strcat('channel_',num2str(ch));
                    [status, msg] = mkdir(aesPath,chFolder);
                    if (~status)
                        return;
                    end
                    aesChPath = fullfile(aesPath,chFolder);

                    [status, msg] = mkdir(smplPath,chFolder);
                    if (~status)
                        return;
                    end
                    smplChPath = fullfile(smplPath,chFolder);
                else
                    aesChPath = aesPath;
                    smplChPath = smplPath;
                end

                for sl=1:tempSlices
                    if (params.volume)
                        slFolder = strcat('slice_',num2str(sl));
                        [status, msg] = mkdir(aesChPath,slFolder);
                        if (~status)
                            return;
                        end
                        aesSlPath = fullfile(aesChPath,slFolder);
    
                        [status, msg] = mkdir(smplChPath,slFolder);
                        if (~status)
                            return;
                        end
                        smplSlPath = fullfile(smplChPath,slFolder);
                    else
                        aesSlPath = aesChPath;
                        smplSlPath = smplChPath;
                    end

                    if (~isempty(params.roiAes{sl,ch}))
                        for ii=1:size(params.roiAes{sl,ch},3)
                            imwrite((squeeze(params.roiAes{sl,ch}(:,:,ii))+0),fullfile(aesSlPath,strcat(params.aesNames{sl,ch}{ii},'.tif')),'Compression','none');
                        end
                    end

                    if (~isempty(params.roiSmpl{sl,ch}))
                        for ii=1:size(params.roiSmpl{sl,ch},3)
                            imwrite((squeeze(params.roiSmpl{sl,ch}(:,:,ii))+0),fullfile(smplSlPath,strcat(params.smplNames{sl,ch}{ii},'.tif')),'Compression','none');
                        end
                    end
                end
            end
            msg = [];
        end
    end
end