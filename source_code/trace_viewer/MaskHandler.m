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
        roiSmplPrev
        roiSmplRef
        roiSmplProp
        roiSmplPropPrev
        roiSmplNames
        numRoiSmpl

        roiAes
        roiAesPrev
        roiAesRef
        roiAesProp
        roiAesPropPrev
        roiAesNames
        numRoiAes
    end

    properties (GetAccess = public)
        dim
        splitCh
        channels
        channel

        roiAesImage
        roiSmplImage
        roiSlctImage
        roiAesSlct
        roiSmplSlct

        smplHitTest
        aesHitTest
    end

    methods
        function obj = MaskHandler(splitChannels, channels, channel, dim)
            obj.splitCh = splitChannels;
            obj.dim = dim;
            obj.channel = channel;
            obj.channels = channels;
            obj.roiSmplNames = cell(obj.channels,1);
            obj.roiAesNames = cell(obj.channels,1);
            obj.roiSmpl = cell(obj.channels,1);
            obj.roiSmplRef = cell(obj.channels,1);
            obj.roiAes = cell(obj.channels,1);
            obj.roiAesRef = cell(obj.channels,1);
            obj.expMask = zeros(obj.dim(1),obj.dim(2),obj.channels);
            obj.expMaskRef = zeros(obj.dim(1),obj.dim(2),obj.channels);
            obj.roiSmplProp = cell(obj.channels,1);
            obj.roiAesProp = cell(obj.channels,1);
            obj.expMaskProp = zeros(obj.channels,3);
            for ch=1:obj.channels
                obj.roiSmplNames{ch} = {};
                obj.roiAesNames{ch} = {};
                obj.roiSmpl{ch} = ([]>0);
                obj.roiSmplRef{ch} = [];
                obj.roiAes{ch} = ([]>0);
                obj.roiAesRef{ch} = [];
                obj.roiSmplProp{ch} = [];
                obj.roiAesProp{ch} = [];
            end
            obj.numRoiSmpl = zeros(obj.channels,1);
            obj.numRoiAes = zeros(obj.channels,1);
            obj.smplHitTest = true;
            obj.aesHitTest = true;
        end

        %% functions for adding/removing/renaming masks
        function addRoiSmpl(obj, name, mask)
            obj.roiSmplRef{obj.channel} = cat(3,obj.roiSmplRef{obj.channel},mask);
            obj.roiSmpl{obj.channel} = cat(3,obj.roiSmpl{obj.channel},mask);
            obj.roiSmplNames{obj.channel} = [obj.roiSmplNames{obj.channel},name];
            obj.roiSmplProp{obj.channel} = [obj.roiSmplProp{obj.channel};[0,0,0]];
            obj.numRoiSmpl(obj.channel) = obj.numRoiSmpl(obj.channel) + 1;
            obj.roiSmplSlct = [obj.roiSmplSlct, 0];
        end

        function addRoiAes(obj, name, mask)
            obj.roiAesRef{obj.channel} = cat(3,obj.roiAesRef{obj.channel},mask);
            obj.roiAes{obj.channel} = cat(3,obj.roiAes{obj.channel},mask);
            obj.roiAesNames{obj.channel} = [obj.roiAesNames{obj.channel},name];
            obj.roiAesProp{obj.channel} = [obj.roiAesProp{obj.channel};[0,0,0]];
            obj.numRoiAes(obj.channel) = obj.numRoiAes(obj.channel) + 1;
            obj.roiAesSlct = [obj.roiAesSlct, 0];
        end

        function removeRoiSmpl(obj,toRemove)
            obj.roiSmplRef{obj.channel} = obj.roiSmplRef{obj.channel}(:,:,~toRemove);
            obj.roiSmpl{obj.channel} = obj.roiSmpl{obj.channel}(:,:,~toRemove);
            obj.roiSmplNames{obj.channel} = obj.roiSmplNames{obj.channel}(~toRemove);
            obj.roiSmplProp{obj.channel} = obj.roiSmplProp{obj.channel}(~toRemove,:);
            obj.numRoiSmpl(obj.channel) = obj.numRoiSmpl(obj.channel) - sum(toRemove);
            obj.roiSmplSlct = obj.roiSmplSlct(~toRemove);
        end

        function removeRoiAes(obj,toRemove)
            obj.roiAesRef{obj.channel} = obj.roiAesRef{obj.channel}(:,:,~toRemove);
            obj.roiAes{obj.channel} = obj.roiAes{obj.channel}(:,:,~toRemove);
            obj.roiAesNames{obj.channel} = obj.roiAesNames{obj.channel}(~toRemove);
            obj.roiAesProp{obj.channel} = obj.roiAesProp{obj.channel}(~toRemove,:);
            obj.numRoiAes(obj.channel) = obj.numRoiAes(obj.channel) - sum(toRemove);
            obj.roiAesSlct = obj.roiAesSlct(~toRemove);
        end

        function removeExpMask(obj)
            obj.expMask(:,:,obj.channel) = obj.expMask(:,:,obj.channel) * 0;
        end

        function renameRoiSmpl(obj,ind,name)
            obj.roiSmplNames{obj.channel}(ind) = name;
        end

        function renameRoiAes(obj,ind,name)
            obj.roiAesNames{obj.channel}(ind) = name;
        end

        function setExpMask(obj,mask)
            obj.expMask(:,:,obj.channel) = mask;
            obj.expMaskRef(:,:,obj.channel) = mask;
            obj.expMaskProp(obj.channel,:) = [0 0 0];
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
            if (obj.numRoiSmpl(obj.channel)>0 && obj.smplHitTest)
                if (cntl)
                    obj.roiSmplSlct = (obj.roiSmplSlct + squeeze(obj.roiSmpl{obj.channel}(coord(2),coord(1),:))) > 0;
                else
                    obj.roiSmplSlct = squeeze(obj.roiSmpl{obj.channel}(coord(2),coord(1),:));
                end

                if (obj.numRoiAes(obj.channel)>0 && obj.aesHitTest)
                    if (sum(obj.roiSmplSlct)==0)
                        if (cntl)
                            obj.roiAesSlct = (obj.roiAesSlct + squeeze(obj.roiAes{obj.channel}(coord(2),coord(1),:))) > 0;
                        else
                            obj.roiAesSlct = squeeze(obj.roiAes{obj.channel}(coord(2),coord(1),:));
                        end
                    else
                        obj.roiAesSlct = 0 * obj.roiAesSlct;
                    end
                end
            elseif(obj.numRoiAes(obj.channel)>0 && obj.aesHitTest)
                if (cntl)
                    obj.roiAesSlct = (obj.roiAesSlct + squeeze(obj.roiAes{obj.channel}(coord(2),coord(1),:))) > 0;
                else
                    obj.roiAesSlct = squeeze(obj.roiAes{obj.channel}(coord(2),coord(1),:));
                end
            end
        end

        % called by main app to build image overlay highlighting masks
        function buildRoiImage(obj)
            if (isempty(obj.roiAes{obj.channel}))
                obj.roiAesImage = zeros(obj.dim);
                aesSlctImage = zeros(obj.dim);
            else
                obj.roiAesImage = sum(obj.roiAes{obj.channel},3);
                obj.roiAesImage = obj.roiAesImage > 0;
                obj.roiAesSlct = obj.roiAesSlct > 0;
                aesSlctImage = sum(obj.roiAes{obj.channel}(:,:,obj.roiAesSlct),3);
            end

            if (isempty(obj.roiSmpl{obj.channel}))
                obj.roiSmplImage = zeros(obj.dim);
                smplSlctImage = zeros(obj.dim);
            else
                obj.roiSmplImage = sum(obj.roiSmpl{obj.channel},3);
                obj.roiSmplImage = obj.roiSmplImage > 0;
                obj.roiSmplSlct = obj.roiSmplSlct > 0;
                smplSlctImage = sum(obj.roiSmpl{obj.channel}(:,:,obj.roiSmplSlct),3);
            end
            obj.roiSlctImage = (aesSlctImage + smplSlctImage) > 0;
        end

        %% functions for adjusting mask position/sizes
        % hold current state incase you want to revert after making changes
        function storeState(obj)
            obj.roiSmplPropPrev = obj.roiSmplProp{obj.channel};
            obj.roiAesPropPrev = obj.roiAesProp{obj.channel};
            obj.expMaskPropPrev = obj.expMaskProp(obj.channel,:);
            obj.roiSmplPrev = obj.roiSmpl{obj.channel};
            obj.roiAesPrev = obj.roiAes{obj.channel};
            obj.expMaskPrev = obj.expMask(:,:,obj.channel);
        end

        % revert masks to previous stored state
        function revert(obj)
            obj.roiSmplProp{obj.channel} = obj.roiSmplPropPrev;
            obj.roiAesProp{obj.channel} = obj.roiAesPropPrev;
            obj.expMaskProp(obj.channel,:) = obj.expMaskPropPrev;
            obj.roiSmpl{obj.channel} = obj.roiSmplPrev;
            obj.roiAes{obj.channel} = obj.roiAesPrev;
            obj.expMask(:,:,obj.channel) = obj.expMaskPrev;
        end

        % called by main app to redraw masks
        function setMaskProp(obj, x, y, mgn, index, type)
            switch type
                case MaskType.EXPMASK
                    obj.expMaskProp(obj.channel,:) = [x,y,mgn];
                    obj.expMask(:,:,obj.channel) = obj.redrawMask(obj.expMaskRef(:,:,obj.channel),x,y,mgn);
                case MaskType.AESALL
                    obj.roiAesProp{obj.channel} = repmat([x,y,mgn],obj.numRoiAes(obj.channel),1);
                    for ii=1:obj.numRoiAes(obj.channel)
                        obj.roiAes{obj.channel}(:,:,ii) = obj.redrawMask(obj.roiAesRef{obj.channel}(:,:,ii),x,y,mgn);
                    end
                case MaskType.AES
                    obj.roiAesProp{obj.channel}(index,:) = [x,y,mgn];
                    obj.roiAes{obj.channel}(:,:,index) = obj.redrawMask(obj.roiAesRef{obj.channel}(:,:,index),x,y,mgn);
                case MaskType.SMPLALL
                    obj.roiSmplProp{obj.channel} = repmat([x,y,mgn],obj.numRoiSmpl(obj.channel),1);
                    for ii=1:obj.numRoiSmpl(obj.channel)
                        obj.roiSmpl{obj.channel}(:,:,ii) = obj.redrawMask(obj.roiSmplRef{obj.channel}(:,:,ii),x,y,mgn);
                    end
                case MaskType.SMPL
                    obj.roiSmplProp{obj.channel}(index,:) = [x,y,mgn];
                    obj.roiSmpl{obj.channel}(:,:,index) = obj.redrawMask(obj.roiSmplRef{obj.channel}(:,:,index),x,y,mgn);
            end
        end

        function flattenAesProp(obj)
            prop = max(obj.roiAesProp{obj.channel});
            obj.roiAesProp{obj.channel} = repmat(prop,obj.numRoiAes(obj.channel),1);
            obj.roiAesSlct = 0*obj.roiAesSlct+1;
            for ii=1:obj.numRoiAes(obj.channel)
                obj.roiAes{obj.channel}(:,:,ii) = obj.redrawMask(obj.roiAesRef{obj.channel}(:,:,ii),prop(1),prop(2),prop(3));
            end
        end

        function flattenSmplProp(obj)
            prop = max(obj.roiSmplProp{obj.channel});
            obj.roiSmplProp{obj.channel} = repmat(prop,obj.numRoiSmpl(obj.channel),1);
            obj.roiSmplSlct = 0*obj.roiSmplSlct+1;
            for ii=1:obj.numRoiSmpl(obj.channel)
                obj.roiSmpl{obj.channel}(:,:,ii) = obj.redrawMask(obj.roiSmplRef{obj.channel}(:,:,ii),prop(1),prop(2),prop(3));
            end
        end

        function enableRoiSmpl(obj,enable)
            obj.smplHitTest = enable;
        end

        function enableRoiAes(obj,enable)
            obj.aesHitTest = enable;
        end

        %% Getters and Setters
        function present = hasExpMask(obj)
            present = sum(obj.expMask(:,:,obj.channel),'all')>0;
        end

        function rois = getRoiAes(obj)
            rois = obj.roiAes{obj.channel};
        end

        function rois = getRoiSmpl(obj)
            rois = obj.roiSmpl{obj.channel};
        end

        function names = getAesNames(obj)
            names = obj.roiAesNames{obj.channel};
        end

        function name = getAesName(obj,index)
            name = obj.roiAesNames{obj.channel}(index);
        end

        function names = getSmplNames(obj)
            names = obj.roiSmplNames{obj.channel};
        end

        function name = getSmplName(obj,index)
            name = obj.roiSmplNames{obj.channel}(index);
        end

        function prop = getAesProp(obj,ind)
            prop = obj.roiAesProp{obj.channel}(ind,:);
        end

        function prop = getSmplProp(obj,ind)
            prop = obj.roiSmplProp{obj.channel}(ind,:);
        end

        function prop = getExpProp(obj)
            prop = obj.expMaskProp(obj.channel,:);
        end

        function mask = getExpMask(obj)
            mask = obj.expMask(:,:,obj.channel);
        end

        function num = getNumAes(obj)
            num = obj.numRoiAes(obj.channel);
        end

        function num = getNumSmpl(obj)
            num = obj.numRoiSmpl(obj.channel);
        end 
        
        %% Other Functions
        % Change color channel
        function setChannel(obj,channel)
            if (obj.splitCh)
                obj.channel = channel;
                obj.roiAesSlct = zeros(obj.numRoiAes(obj.channel),1);
                obj.roiSmplSlct = zeros(obj.numRoiSmpl(obj.channel),1);
            else
                obj.channel = 1;
            end
        end

        % attach masks to params struct for processing files
        function params = attachParams(obj,params)
            for ch = 1:obj.channel
                obj.roiAes{ch} = (obj.roiAes{ch} > 0);
                obj.roiSmpl{ch} = (obj.roiSmpl{ch} > 0);
            end
            
            if (obj.splitCh)
                fullmask = zeros(size(obj.expMask));
                tempCh = obj.channel;
                for ch = 1:obj.channels
                    obj.setChannel(ch);
                    obj.buildRoiImage();
                    fullmask(:,:,ch) = (obj.roiAesImage+obj.roiSmplImage+obj.roiSlctImage+obj.expMask(:,:,ch))>0;
                end
                obj.setChannel(tempCh);
            else
                fullmask = zeros(size(obj.expMask));
                obj.buildRoiImage();
                fullmask(:,:,1) = (obj.roiAesImage+obj.roiSmplImage+obj.roiSlctImage+obj.expMask(:,:,1))>0;
                for ch=2:obj.channels
                    obj.roiAes{ch} = obj.roiAes{1};
                    obj.roiSmpl{ch} = obj.roiSmpl{1};
                    fullmask(:,:,ch) = fullmask(:,:,1);
                    obj.roiAesNames{ch} = obj.roiAesNames{1};
                    obj.roiSmplNames{ch} = obj.roiSmplNames{1};
                end                
            end
            params.roiAes = obj.roiAes;
            params.roiSmpl = obj.roiSmpl;
            params.splitChannels = obj.splitCh;
            params.expMask = fullmask;
            params.aesNames = obj.roiAesNames;
            params.smplNames = obj.roiSmplNames;
            params.dim = obj.dim;
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
    end
end