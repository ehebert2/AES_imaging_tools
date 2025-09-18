classdef MaskListenerHandler < handle
    properties (SetAccess = private, GetAccess = private)
        rois
        roiMasks
        lineWidths
        ringWidths
        roiTypes
        xMap
        yMap
        
        channels
        slices
        splitChannel

        imageAxes
        callingApp
        imageDim
        compMasks
        roiMask
        currentROI
        

        writeDim1
        writeDim2
    end

    properties (Constant)
        NONE = -1;
        RECT = 1;
        CIRC = 2;
        ELPS = 3;
        RING = 4;
        LINE = 5;
        PLGN = 6;
    end

    properties (GetAccess = public)
        ROIParamLabels
        currentROIType
        currentROIInd
        numRois
        channel
        slice
    end

    properties (GetAccess = public, SetAccess = public)
        powerMod
    end

    methods
        function obj = MaskListenerHandler(callingApp,imageAxes,params)
            obj.callingApp = callingApp;
            obj.imageAxes = imageAxes;
            if (params.volume)
                obj.slices = params.slices;
            else
                obj.slices = 1;
            end
            obj.splitChannel = params.splitChannels;
            if (obj.splitChannel)
                obj.channels = params.channels;
            else
                obj.channels = 1;
            end
            obj.imageDim = [0,0];
            obj.imageDim(1) = params.dim(1);
            obj.imageDim(2) = params.dim(2);
            obj.roiMask = zeros(obj.imageDim(1),obj.imageDim(2));
            obj.compMasks = zeros(obj.imageDim(1),obj.imageDim(2),obj.slices,obj.channels) > 0;
            obj.numRois = zeros(obj.slices,obj.channels);
            obj.xMap = 1:obj.imageDim(2);
            obj.xMap = repmat(obj.xMap,obj.imageDim(1),1);
            obj.yMap = 1:obj.imageDim(1);
            obj.yMap = repmat(obj.yMap',1,obj.imageDim(2));
            obj.channel = 1;
            obj.slice = 1;
            
            obj.rois = cell(obj.slices,obj.channels);
            obj.roiTypes = cell(obj.slices,obj.channels);
            obj.lineWidths = cell(obj.slices,obj.channels);
            obj.ringWidths = cell(obj.slices,obj.channels);
            obj.roiMasks = cell(obj.slices,obj.channels);
            for sl=1:obj.slices
                for ch=1:obj.channels
                    obj.rois{sl,ch} = cell(0);
                    obj.roiTypes{sl,ch} = [];
                    obj.lineWidths{sl,ch} = [];
                    obj.ringWidths{sl,ch} = [];
                    obj.roiMasks{sl,ch} = []>0;
                end
            end

            obj.ROIParamLabels = cell(7,1);
            obj.ROIParamLabels{obj.RECT} = {'X min:', 'Y min:', 'Width:', 'Height:'};
            obj.ROIParamLabels{obj.CIRC} = {'X center:','Y center:','Radius:'};
            obj.ROIParamLabels{obj.ELPS} = {'X center:', 'Y center:', 'Semi axis 1:', 'Semi axis 2:', 'Angle:'};
            obj.ROIParamLabels{obj.RING} = {'X center:', 'Y center:', 'Semi axis 1:', 'Semi axis 2:', 'Angle:', 'Thickness'};
            obj.ROIParamLabels{obj.LINE} = {'Thickness:','Points'};
            obj.ROIParamLabels{obj.PLGN} = {'Points:'};
            obj.currentROIType = obj.NONE;
        end

        function delete(obj)
            for sl=1:obj.slices
                for ch=1:obj.channels
                    for ii=1:length(obj.rois{sl,ch})
                        delete(obj.rois{sl,ch}{ii});
                    end
                    obj.roiTypes{sl,ch} = [];
                    obj.lineWidths{sl,ch} = [];
                    obj.ringWidths{sl,ch} = [];
                end
            end
        end

        function addROI(obj,roiType)
            switch(roiType)
                case(obj.CIRC)
                    tempROI = drawcircle(obj.imageAxes);
                    if (isempty(tempROI.Center))
                        return;
                    end
                    tempROI.FaceAlpha = 0;
                    obj.roiMask = createMask(tempROI);
                case(obj.RECT)
                    tempROI = drawrectangle(obj.imageAxes);
                    if (isempty(tempROI.Position))
                        return;
                    end
                    tempROI.FaceAlpha = 0;
                    obj.roiMask = createMask(tempROI);
                case(obj.ELPS)
                    tempROI = drawellipse(obj.imageAxes);
                    if (isempty(tempROI.Center))
                        return;
                    end
                    tempROI.FaceAlpha = 0;
                    obj.roiMask = createMask(tempROI);
                case(obj.LINE)
                    tempROI = drawpolyline(obj.imageAxes);
                    if (isempty(tempROI.Position))
                        return;
                    end
                    addlistener(tempROI,'VertexAdded',@obj.roiMv);
                    addlistener(tempROI,'VertexDeleted',@obj.roiMv);
                    obj.drawLine(tempROI,1);
                case(obj.RING)
                    tempROI = drawellipse(obj.imageAxes);
                    if (isempty(tempROI.Center))
                        return;
                    end
                    tempROI.FaceAlpha = 0;
                    obj.drawRing(tempROI,1);
                case(obj.PLGN)
                    tempROI = drawpolygon(obj.imageAxes);
                    if (isempty(tempROI.Position))
                        return;
                    end
                    tempROI.FaceAlpha = 0;
                    addlistener(tempROI,'VertexAdded',@obj.roiMv);
                    addlistener(tempROI,'VertexDeleted',@obj.roiMv);
                    obj.roiMask = createMask(tempROI);
            end
            addlistener(tempROI,'ROIClicked',@obj.roiClck);
            addlistener(tempROI,'DeletingROI',@obj.roiDlt);
            addlistener(tempROI,'ROIMoved',@obj.roiMv);
            obj.rois{obj.slice,obj.channel} = [obj.rois{obj.slice,obj.channel}, {tempROI}];
            obj.roiTypes{obj.slice,obj.channel} = [obj.roiTypes{obj.slice,obj.channel} roiType];
            obj.lineWidths{obj.slice,obj.channel} = [obj.lineWidths{obj.slice,obj.channel},1];
            obj.ringWidths{obj.slice,obj.channel} = [obj.ringWidths{obj.slice,obj.channel}, 1];
            tempROI.LineWidth = 1;
            tempROI.DrawingArea = "unlimited";
            obj.currentROI = tempROI;
            obj.currentROIType = roiType;
            obj.numRois(obj.slice,obj.channel) = obj.numRois(obj.slice,obj.channel) + 1;
            obj.currentROIInd = obj.numRois(obj.slice,obj.channel);
            obj.roiMasks{obj.slice,obj.channel} = cat(3,obj.roiMasks{obj.slice,obj.channel},(obj.roiMask>0));
            
            obj.callingApp.setROILabels(obj.ROIParamLabels{obj.currentROIType},obj.currentROIType);
            obj.callingApp.enableROIFields(true);
            obj.setFields();
            obj.buildCompMask();
            obj.powerMod.addROI(obj.slice,obj.channel);
        end

        function deleteCurrentROI(obj)
            tempROI = obj.currentROI;
            obj.roiDlt(tempROI,[]);
            delete(tempROI);
        end

        function setFields(obj)
            switch (obj.currentROIType)
                case (obj.NONE)
                    vals = [];
                case (obj.CIRC)
                    vals = zeros(3,1);
                    coord = obj.currentROI.Center;
                    vals(1) = coord(1);
                    vals(2) = coord(2);
                    vals(3) = obj.currentROI.Radius;
                case (obj.RECT)
                    vals = obj.currentROI.Position;
                case (obj.ELPS)
                    vals = zeros(5,1);
                    coord = obj.currentROI.Center;
                    vals(1) = coord(1);
                    vals(2) = coord(2);
                    coord = obj.currentROI.SemiAxes;
                    vals(3) = coord(1);
                    vals(4) = coord(2);
                    vals(5) = obj.currentROI.RotationAngle;
                case (obj.LINE)
                    vals = zeros(2,1);
                    vals(1) = obj.lineWidths{obj.slice,obj.channel}(obj.currentROIInd);
                    coord = obj.currentROI.Position;
                    vals(2) = length(coord);
                case (obj.RING)
                    vals = zeros(6,1);
                    coord = obj.currentROI.Center;
                    vals(1) = coord(1);
                    vals(2) = coord(2);
                    
                    coord = obj.currentROI.SemiAxes;
                    vals(3) = coord(1);
                    vals(4) = coord(2);
                    vals(5) = obj.currentROI.RotationAngle;
                    vals(6) = obj.ringWidths{obj.slice,obj.channel}(obj.currentROIInd);
                case (obj.PLGN)
                    coord = obj.currentROI.Position;
                    vals = zeros(1,1);
                    vals(1) = length(coord);
            end
            obj.callingApp.setROIFields(vals);
        end

        function setChannel(obj,channel)
            if (obj.splitChannel && (channel~=obj.channel))
                obj.setROIsVisible(obj.slice,obj.channel,'off');
                obj.channel = channel;
                obj.setROIsVisible(obj.slice,obj.channel,'on');
                obj.currentROI = [];
                obj.deselectROI();
                obj.callingApp.enableROIFields(false);
            end
        end

        function setSlice(obj,slice)
            if (obj.slice~=slice)
                obj.setROIsVisible(obj.slice,obj.channel,'off');
                obj.slice = slice;
                obj.setROIsVisible(obj.slice,obj.channel,'on');
                obj.currentROI = [];
                obj.deselectROI();
                obj.callingApp.enableROIFields(false);
            end
        end

        function deselectROI(obj)
            obj.currentROI = [];
            obj.currentROIType = obj.NONE;
            obj.powerMod.chngROI();
        end

        function val = chngParam(obj,val,ind)
            switch (obj.currentROIType)
                case (obj.RECT)
                    coord = obj.currentROI.Position;
                    coord(ind) = val;
                    obj.currentROI.Position = coord;
                case (obj.CIRC)
                    if (ind<3)
                        coord = obj.currentROI.Center;
                        coord(ind) = val;
                        obj.currentROI.Center = coord;
                    else
                        val = abs(val);
                        obj.currentROI.Radius = val;
                    end
                case (obj.ELPS)
                    if (ind<3)
                        coord = obj.currentROI.Center;
                        coord(ind) = val;
                        obj.currentROI.Center = coord;
                    elseif (ind<5)
                        coord = obj.currentROI.SemiAxes;
                        coord(ind-2) = val;
                        obj.currentROI.SemiAxes = coord;
                    else
                        obj.currentROI.RotationAngle = val;
                    end
                case (obj.RING)
                    if (ind<3)
                        coord = obj.currentROI.Center;
                        coord(ind) = val;
                        obj.currentROI.Center = coord;
                    elseif (ind<5)
                        coord = obj.currentROI.SemiAxes;
                        coord(ind-2) = val;
                        obj.currentROI.SemiAxes = coord;
                    elseif (ind==5)
                        obj.currentROI.RotationAngle = val;
                    else
                        if (val<0)
                            val=0;
                        end
                        obj.ringWidths{obj.slice,obj.channel}(obj.currentROIInd) = val;
                    end
                case (obj.LINE)
                    if (ind==1)
                        val = abs(val);
                        val = 2*ceil((val+1)/2)-1;
                        obj.lineWidths{obj.slice,obj.channel}(obj.currentROIInd) = val;
                    else
                        if (val<1)
                            val=1;
                        end
                        val=round(val);
                        obj.currentROI.Position = obj.chngNumPnts(obj.currentROI.Position,val);
                    end
                case (obj.PLGN)
                    if (val<1)
                        val=1;
                    end
                    val=round(val);
                    obj.currentROI.Position = obj.chngNumPnts(obj.currentROI.Position,val);
            end
            obj.redrawROIMask(obj.currentROIInd);
            obj.buildCompMask();
        end

        % called when arrow buttons are used to move mask
        function shiftMask(obj,dspl)
            if (obj.currentROIType~=obj.NONE)
                switch obj.currentROIType
                    case obj.CIRC
                        center = obj.currentROI.Center;
                        rad = obj.currentROI.Radius/sqrt(2);
                        bounds = [-rad,rad,-rad,rad];
                        if (~obj.checkInbounds(dspl,center,bounds))
                            return;
                        end
                        center = center + dspl;
                        obj.currentROI.Center = center;
                        obj.callingApp.setROIFields(center);
                    case obj.RECT
                        pos = obj.currentROI.Position;
                        center = [pos(1),pos(2)];
                        bounds = [0,pos(3),0,pos(4)];
                        if (~obj.checkInbounds(dspl,center,bounds))
                            return;
                        end
                        pos(1) = pos(1) + dspl(1);
                        pos(2) = pos(2) + dspl(2);
                        obj.currentROI.Position = pos;
                        obj.callingApp.setROIFields([pos(1),pos(2)]);
                    case {obj.ELPS, obj.RING}
                        center = obj.currentROI.Center;
                        rad = min(obj.currentROI.SemiAxes)/sqrt(2);
                        bounds = [-rad,rad,-rad,rad];
                        if (~obj.checkInbounds(dspl,center,bounds))
                            return;
                        end
                        center = center + dspl;
                        obj.currentROI.Center = center;
                        obj.callingApp.setROIFields(center);
                    case {obj.LINE, obj.PLGN}
                        pos = obj.currentROI.Position;
                        center = mean(pos);
                        bounds = [0, 0, 0, 0];
                        if (~obj.checkInbounds(dspl,center,bounds))
                            return;
                        end
                        obj.currentROI.Position = pos + repmat(dspl,size(pos,1),1);
                end
                obj.redrawROIMask(obj.currentROIInd);
                obj.buildCompMask();
            end
        end

        function duplicateMask(obj,ch,sl1,sl2)
            if (obj.currentROIType~=obj.NONE)
                for sl=sl1:sl2
                    if ((ch==obj.channel) && (sl==obj.slice))
                        continue;
                    end

                    tempROI = copyobj(obj.currentROI,obj.imageAxes);
                    obj.rois{sl,ch} = [obj.rois{sl,ch}, {tempROI}];
                    obj.roiMasks{sl,ch} = cat(3,obj.roiMasks{sl,ch},obj.roiMasks{obj.slice,obj.channel}(:,:,obj.currentROIInd));
                    obj.roiTypes{sl,ch} = [obj.roiTypes{sl,ch},obj.currentROIType];
                    obj.lineWidths{sl,ch} = [obj.lineWidths{sl,ch},obj.lineWidths{obj.slice,obj.channel}(obj.currentROIInd)];
                    obj.ringWidths{sl,ch} = [obj.ringWidths{sl,ch}, obj.ringWidths{obj.slice,obj.channel}(obj.currentROIInd)];
                    addlistener(tempROI,'ROIClicked',@obj.roiClck);
                    addlistener(tempROI,'DeletingROI',@obj.roiDlt);
                    addlistener(tempROI,'ROIMoved',@obj.roiMv);
                    if ((obj.currentROIType==obj.LINE)||(obj.currentROIType==obj.PLGN))
                        addlistener(tempROI,'VertexAdded',@obj.roiMv);
                        addlistener(tempROI,'VertexDeleted',@obj.roiMv);
                    end
                    tempROI.Visible = 'off';
                    obj.compMasks(:,:,sl,ch) = (sum(obj.roiMasks{sl,ch},3)>0);
                end
            end
        end

        %% Main function for building mask
        function buildCompMask(obj)
            obj.compMasks(:,:,obj.slice,obj.channel) = (sum(obj.roiMasks{obj.slice,obj.channel},3)>0);
            obj.callingApp.updateOcc();
            obj.callingApp.updateDisplay();
        end

        %% functions for displaying mask/occupancy information
        function mask = getMask(obj)
            mask = obj.compMasks(:,:,obj.slice,obj.channel);
        end

        function mask = getSelectedMask(obj)
            if (obj.currentROIType == obj.NONE)
                mask = 0*obj.roiMask;
            else
                mask = obj.roiMasks{obj.slice,obj.channel}(:,:,obj.currentROIInd);
            end
        end

        function mask = getRoiMask(obj,sl,ch,ind)
            mask = obj.roiMasks{sl,ch}(:,:,ind);
        end

        function mask = getChMask(obj,sl,ch)
            mask = obj.compMasks(:,:,sl,ch);
        end

        %% File IO
        function success = loadROIs(obj,data)
            success = false;
            obj.writeDim1=0;
            obj.writeDim2=0;
            if (~isfield(data,'sepCh'))
                return;
            end

            sepCh = data.sepCh;
            for sl=1:obj.slices
                for ch=1:obj.channels
                    obj.rois{sl,ch} = cell(0);
                    obj.roiTypes{sl,ch} = [];
                    obj.lineWidths{sl,ch} = [];
                    obj.ringWidths{sl,ch} = [];
                end
            end

            if (isfield(data,'rois'))
                obj.copyROIs(data.rois,obj.NONE,1);
                if (~isfield(data,'types') || ~isfield(data,'lws') || ~isfield(data,'rws'))
                    return;
                end

                obj.roiTypes = obj.copyWidths(data.types);
                obj.lineWidths = obj.copyWidths(data.lws);
                obj.ringWidths = obj.copyWidths(data.rws);
            else
                if (isfield(data,'circ'))
                    obj.copyROIs(data.circ,obj.CIRC,1);
                end
    
                if (isfield(data,'rect'))
                    obj.copyROIs(data.rect,obj.RECT,1);
                end
    
                if (isfield(data,'elps'))
                    obj.copyROIs(data.elps,obj.ELPS,1);
                end
    
                if (isfield(data,'ring') && isfield(data,'rws'))
                    obj.copyROIs(data.ring,obj.RING,data.rws);
                end
    
                if (isfield(data,'line') && isfield(data,'lws'))
                    obj.copyROIs(data.line,obj.LINE,data.lws);
                end
    
                if (isfield(data,'plgn'))
                    obj.copyROIs(data.plgn,obj.PLGN,1);
                end
            end

            if ((~sepCh) && obj.splitChannel)
                for sl=obj.slices
                    for ch=2:obj.channels
                        obj.roiTypes{sl,ch} = obj.roiTypes{sl,1};
                        obj.lineWidths{sl,ch} = obj.lineWidths{sl,1};
                        obj.ringWidths{sl,ch} = obj.ringWidths{sl,1};
                        obj.rois{sl,ch} = cell(1,length(obj.rois{sl,1}));
                        for ii=1:length(obj.rois{sl,1})
                            obj.rois{sl,ch}{ii} = copyobj(obj.rois{sl,1}{ii},obj.imageAxes);
                        end
                    end
                end
            end

            tempCh = obj.channel;
            tempSl = obj.slice;
            for sl=1:obj.slices
                obj.slice = sl;
                for ch=1:obj.channels
                    obj.channel = ch;
                    if (isempty(obj.rois{sl,ch}))
                        obj.roiMasks{sl,ch} = [];
                        obj.compMasks(:,:,sl,ch) = 0*obj.compMasks(:,:,sl,ch);
                    else
                        obj.roiMasks{sl,ch} = (zeros(obj.imageDim(1),obj.imageDim(2),obj.numRois(sl,ch))>0);
                        for ii=1:length(obj.rois{sl,ch})
                            obj.rois{sl,ch}{ii}.Parent = obj.imageAxes;
                            addlistener(obj.rois{sl,ch}{ii},'ROIClicked',@obj.roiClck);
                            addlistener(obj.rois{sl,ch}{ii},'DeletingROI',@obj.roiDlt);
                            addlistener(obj.rois{sl,ch}{ii},'ROIMoved',@obj.roiMv);
                            if ((obj.roiTypes{sl,ch}(ii)==obj.LINE)||(obj.roiTypes{sl,ch}(ii)==obj.PLGN))
                                addlistener(obj.rois{sl,ch}{ii},'VertexAdded',@obj.roiMv);
                                addlistener(obj.rois{sl,ch}{ii},'VertexDeleted',@obj.roiMv);
                            end
                            obj.redrawROIMask(ii);
                        end
                        obj.compMasks(:,:,sl,ch) = (sum(obj.roiMasks{sl,ch},3)>0);

                        if ((ch~=tempCh)||(sl~=tempSl))
                            obj.setROIsVisible(sl,ch,'off');
                        end
                    end
                end
            end
            obj.channel = tempCh;
            obj.slice = tempSl;
            success = true;
        end

        function saveROIs(obj,fname)
            for sl=1:obj.slices
                for ch=1:obj.channels
                    if (obj.channel~=ch || obj.slice~=sl)
                        obj.setROIsVisible(sl,ch,'on');
                    end
                end
            end

            sepCh = obj.splitChannel;
            rois = obj.rois;
            types = obj.roiTypes;
            lws = obj.lineWidths;
            rws = obj.ringWidths;
            save(fname,"sepCh","rois","types","lws","rws");

            for sl=1:obj.slices
                for ch=1:obj.channels
                    if (obj.channel~=ch || obj.slice~=sl)
                        obj.setROIsVisible(sl,ch,'off');
                    end
                end
            end
        end

        function status = printROIMasks(obj,basename,path)
            for ch=1:obj.channels
                if (obj.channels>1)
                    folderTitle = strcat('channel_',num2str(ch));
                    [status, msg] = mkdir(path,folderTitle);
                    if (~status)
                        opts.WindowStyle = 'non-modal';
                        waitfor(errordlg(msg,'Folder Write Issue',opts));
                        return;
                    end
                    chPath = fullfile(path,folderTitle);
                else
                    chPath = path;
                end
                
                for sl=1:obj.slices
                    if (obj.slices>1)
                        folderTitle = strcat('slice_',num2str(sl));
                        [status, msg] = mkdir(chPath,folderTitle);
                        if (~status)
                            opts.WindowStyle = 'non-modal';
                            waitfor(errordlg(msg,'Folder Write Issue',opts));
                            return;
                        end
                        slPath = fullfile(chPath,folderTitle);
                    else
                        slPath = chPath;
                    end

                    for ii=1:obj.numRois(sl,ch)
                        ff = fullfile(slPath,strcat(basename,'_',num2str(ii),'.tif'));
                        imwrite((obj.roiMasks{sl,ch}(:,:,ii)+0),ff,'Compression','none');
                    end
                end
            end
            status = true;
        end

        %% Listeners
        % helper functions used for various listeners
        function chngROI(obj,ind)
            obj.currentROIInd = ind;
            obj.currentROIType = obj.roiTypes{obj.slice,obj.channel}(obj.currentROIInd);
            obj.currentROI = obj.rois{obj.slice,obj.channel}{ind};
            obj.callingApp.setROILabels(obj.ROIParamLabels{obj.currentROIType},obj.currentROIType);
            obj.callingApp.enableROIFields(true);
            obj.setFields();
            obj.powerMod.chngROI();
            obj.callingApp.updateDisplay();
        end

        function redrawROIMask(obj,ind)
            switch (obj.roiTypes{obj.slice,obj.channel}(ind))
                case{obj.RECT,obj.CIRC,obj.ELPS,obj.PLGN}
                    obj.roiMask = createMask(obj.rois{obj.slice,obj.channel}{ind});
                case(obj.RING)
                    obj.drawRing(obj.rois{obj.slice,obj.channel}{ind},obj.ringWidths{obj.slice,obj.channel}(ind));
                case(obj.LINE)
                    obj.drawLine(obj.rois{obj.slice,obj.channel}{ind},obj.lineWidths{obj.slice,obj.channel}(ind));
            end

            obj.roiMasks{obj.slice,obj.channel}(:,:,ind) = obj.roiMask > 0;
            obj.powerMod.redrawROI(obj.slice,obj.channel,ind);
        end

        function roiMv(obj,~,~)
            obj.redrawROIMask(obj.currentROIInd);
            obj.setFields();
            obj.buildCompMask();
        end

        function roiClck(obj,src,~)
            obj.callingApp.deselectThresh();
            for ii=1:length(obj.rois{obj.slice,obj.channel})
                if (eq(src,obj.rois{obj.slice,obj.channel}{ii}))
                    ind = ii;
                    break
                end
            end
            obj.chngROI(ind);
        end

        function roiDlt(obj,src,~)
            if (isempty(src))
                return;
            end
            
            keep = (ones(1,obj.numRois(obj.slice,obj.channel))>0);
            for ii=1:length(obj.rois{obj.slice,obj.channel})
                if (eq(src,obj.rois{obj.slice,obj.channel}{ii}))
                    obj.rois{obj.slice,obj.channel}{ii} = [];
                    keep(ii) = false;
                    break;
                end
            end
            obj.rois{obj.slice,obj.channel} = obj.rois{obj.slice,obj.channel}(keep);
            obj.roiTypes{obj.slice,obj.channel} = obj.roiTypes{obj.slice,obj.channel}(keep);
            obj.lineWidths{obj.slice,obj.channel} = obj.lineWidths{obj.slice,obj.channel}(keep);
            obj.ringWidths{obj.slice,obj.channel} = obj.ringWidths{obj.slice,obj.channel}(keep);
            obj.numRois(obj.slice,obj.channel) = length(obj.rois{obj.slice,obj.channel});
            obj.roiMasks{obj.slice,obj.channel} = obj.roiMasks{obj.slice,obj.channel}(:,:,keep);
            obj.currentROI = [];
            obj.currentROIType = obj.NONE;
            obj.callingApp.enableROIFields(false);
            obj.buildCompMask();
            obj.powerMod.deleteROI(obj.slice,obj.channel,ii);
        end

        %% other helper functions
        function setROIsVisible(obj,slice,channel,vis)
            for ii=1:length(obj.rois{slice,channel})
                obj.rois{slice,channel}{ii}.Visible = vis;
            end
        end

        function copyROIs(obj,arr,type,widths)
            if ((obj.writeDim1==0) || (obj.writeDim2==0))
                if (size(arr,1)>obj.slices)
                    obj.writeDim1 = obj.slices;
                else
                    obj.writeDim1 = size(arr,1);
                end

                if (size(arr,2)>obj.channels)
                    obj.writeDim2 = obj.channels;
                else
                    obj.writeDim2 = size(arr,2);
                end
            end
    
            for sl=1:obj.writeDim1
                for ch=1:obj.writeDim2
                    tempArray = cell(1,(length(obj.rois{sl,ch})+length(arr{sl,ch})));
                    if (~isempty(obj.rois{sl,ch}))
                        tempArray(1:length(obj.rois{sl,ch})) = obj.rois{sl,ch};
                    end

                    if (~isempty(arr{sl,ch}))
                        tempArray((length(obj.rois{sl,ch})+1):end) = arr{sl,ch};
                    end
                    obj.rois{sl,ch} = tempArray;
                    obj.numRois(sl,ch) = length(tempArray);
                    if (type~=obj.NONE)
                        obj.roiTypes{sl,ch} = [obj.roiTypes{sl,ch}, type*ones(1,length(arr{sl,ch}))];
                        if (type==obj.RING)
                            obj.ringWidths{sl,ch} = [obj.ringWidths{sl,ch}, widths{sl,ch}];
                            obj.lineWidths{sl,ch} = [obj.lineWidths{sl,ch}, ones(1,length(arr{sl,ch}))];
                        elseif (type==obj.LINE)
                            obj.lineWidths{sl,ch} = [obj.lineWidths{sl,ch}, widths{sl,ch}];
                            obj.ringWidths{sl,ch} = [obj.ringWidths{sl,ch}, ones(1,length(arr{sl,ch}))];
                        else
                            obj.lineWidths{sl,ch} = [obj.lineWidths{sl,ch}, ones(1,length(arr{sl,ch}))];
                            obj.ringWidths{sl,ch} = [obj.ringWidths{sl,ch}, ones(1,length(arr{sl,ch}))];
                        end
                    end
                end
            end
        end

        function newArray = copyWidths(obj,arr)
            newArray = cell(obj.slices,obj.channels);
            for sl=1:obj.slices
                for ch=1:obj.channels
                    newArray{sl,ch} = [];
                end
            end

            for sl=1:obj.writeDim1
                for ch=1:obj.writeDim2
                    newArray{sl,ch} = arr{sl,ch};
                end
            end
        end

        function newPnts = chngNumPnts(obj,pnts,target)
            diff = target-size(pnts,1);
            if (diff>0)
                newPnts = zeros(target,2);
                newPnts(1:size(pnts,1),:) = pnts;
                for ii=(size(pnts,1)+1):(size(pnts,1)+diff)
                    pt = newPnts((ii-1),:);
                    pt = pt+[15,15];
                    if (pt(1)>obj.imageDim(2))
                        pt(1)=obj.imageDim(2);
                    end

                    if (pt(2)>obj.imageDim(1))
                        pt(2)=obj.imageDim(1);
                    end
                    newPnts(ii,:)=pt;
                end
            else
                newPnts = pnts(1:target,:);
            end
        end

        function drawRing(obj,roi,width)
            ax = roi.SemiAxes;
            cent = roi.Center;
            x = (obj.xMap-cent(1));
            y = (obj.yMap-cent(2));
            theta = -pi/180*roi.RotationAngle;
            obj.roiMask = (1 < obj.ellipseMask(x,y,(ax(1)-width/2),(ax(2)-width/2),theta))...
                .* (1 > obj.ellipseMask(x,y,(ax(1)+width/2),(ax(2)+width/2),theta));
        end

        function drawLine(obj,roi,width)
            obj.roiMask = createMask(roi);
            obj.roiMask = imdilate(obj.roiMask,strel('disk',(width+1)/2));
        end

        function inbounds = checkInbounds(obj,dspl,center,box)
            inbounds = true;
            if ((center(1)+box(2))<1)
                inbounds = dspl(1) >= 0;
            elseif ((center(1)+box(1)+1)>obj.imageDim(2))
                inbounds = dspl(1) <= 0;
            end

            if ((center(2)+box(4))<1)
                inbounds = dspl(2) >= 0;
            elseif ((center(2)+box(3)+1)>obj.imageDim(1))
                inbounds = dspl(2) <= 0;
            end
        end
    end

    methods (Static)
        function shiftROI(roi,shiftVec)
            coord = roi.Center;
            coord = coord + shiftVec;
            roi.Center = coord;
        end

        function d = ellipseMask(x,y,a,b,theta)
            d = (a^2*sin(theta)^2+b^2*cos(theta)^2)*x.^2+2*(b^2-a^2)*sin(theta)*cos(theta)*x.*y+(a^2*cos(theta)^2+b^2*sin(theta)^2)*y.^2;
            d = d/(a^2*b^2);
        end
    end
end