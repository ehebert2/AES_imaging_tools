classdef MaskListenerHandler < handle
    properties (SetAccess = private, GetAccess = private)
        circROIs
        rectROIs
        elpsROIs
        lineROIs
        ringROIs
        plgnROIs
        lineWidths
        ringWidths
        xMap
        yMap
        
        channels
        channel
        slices
        slice
        splitChannel
        isVolume
        bidirectional

        imageAxes
        callingApp
        imageDim
        roiMasks
        roiMask

        currentROI
        currentROIInd

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
        numRois
    end

    properties (GetAccess = public, SetAccess = public)
        bidiShift
        htot
    end

    methods
        function obj = MaskListenerHandler(callingApp,imageAxes,params)
            obj.callingApp = callingApp;
            obj.imageAxes = imageAxes;
            obj.isVolume = params.volume;
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
            obj.bidirectional = params.bidirectional;
            obj.imageDim = [0,0];
            obj.imageDim(1) = params.dim(1);
            obj.imageDim(2) = params.dim(2);
            obj.roiMask = zeros(obj.imageDim(1),obj.imageDim(2));
            obj.roiMasks = zeros(obj.imageDim(1),obj.imageDim(2),obj.slices,obj.channels) > 0;
            obj.xMap = 1:obj.imageDim(2);
            obj.xMap = repmat(obj.xMap,obj.imageDim(1),1);
            obj.yMap = 1:obj.imageDim(1);
            obj.yMap = repmat(obj.yMap',1,obj.imageDim(2));
            obj.channel = 1;
            obj.slice = 1;
            obj.bidiShift = 0;
            obj.htot = 0;
            
            obj.circROIs = cell(obj.slices,obj.channels);
            obj.rectROIs = cell(obj.slices,obj.channels);
            obj.elpsROIs = cell(obj.slices,obj.channels);
            obj.lineROIs = cell(obj.slices,obj.channels);
            obj.ringROIs = cell(obj.slices,obj.channels);
            obj.plgnROIs = cell(obj.slices,obj.channels);
            obj.lineWidths = cell(obj.slices,obj.channels);
            obj.ringWidths = cell(obj.slices,obj.channels);
            for sl=1:obj.slices
                for ch=1:obj.channels
                    obj.rectROIs{sl,ch} = cell(0);
                    obj.circROIs{sl,ch} = cell(0);
                    obj.elpsROIs{sl,ch} = cell(0);
                    obj.ringROIs{sl,ch} = cell(0);
                    obj.lineROIs{sl,ch} = cell(0);
                    obj.plgnROIs{sl,ch} = cell(0);
                    obj.lineWidths{sl,ch} = [];
                    obj.ringWidths{sl,ch} = [];
                end
            end

            obj.ROIParamLabels = cell(7,1);
            obj.ROIParamLabels{obj.RECT} = {'X min:', 'Y min:', 'Width:', 'Height:'};
            obj.ROIParamLabels{obj.CIRC} = {'X center:','Y center:','Radius:'};
            obj.ROIParamLabels{obj.ELPS} = {'X center:', 'Y center:', 'Semi axis 1:', 'Semi axis 2:', 'Angle:'};
            obj.ROIParamLabels{obj.RING} = {'X center:', 'Y center:', 'Semi axis 1:', 'Semi axis 2:', 'Angle:', 'Thickness'};
            obj.ROIParamLabels{obj.LINE} = {'Thickness:','Points'};
            obj.ROIParamLabels{obj.PLGN} = {'Points:'};
        end

        function delete(obj)
            for sl=1:obj.slices
                for ch=1:obj.channels
                    for ii=1:length(obj.rectROIs{sl,ch})
                        delete(obj.rectROIs{sl,ch}{ii});
                    end

                    for ii=1:length(obj.circROIs{sl,ch})
                        delete(obj.circROIs{sl,ch}{ii});
                    end

                    for ii=1:length(obj.elpsROIs{sl,ch})
                        delete(obj.elpsROIs{sl,ch}{ii});
                    end

                    for ii=1:length(obj.ringROIs{sl,ch})
                        delete(obj.ringROIs{sl,ch}{ii});
                    end

                    for ii=1:length(obj.lineROIs{sl,ch})
                        delete(obj.lineROIs{sl,ch}{ii});
                    end

                    for ii=1:length(obj.plgnROIs{sl,ch})
                        delete(obj.plgnROIs{sl,ch}{ii});
                    end

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
                    tempROI.DrawingArea = "unlimited";
                    tempROI.FaceAlpha = 0;
                    obj.circROIs{obj.slice,obj.channel} = [obj.circROIs{obj.slice,obj.channel}, {tempROI}];
                    addlistener(tempROI,'ROIClicked',@obj.circClck);
                    addlistener(tempROI,'DeletingROI',@obj.circDlt);
                    addlistener(tempROI,'ROIMoved',@obj.circMv);
                case(obj.RECT)
                    tempROI = drawrectangle(obj.imageAxes);
                    if (isempty(tempROI.Position))
                        return;
                    end
                    tempROI.DrawingArea = "unlimited";
                    tempROI.FaceAlpha = 0;
                    obj.rectROIs{obj.slice,obj.channel} = [obj.rectROIs{obj.slice,obj.channel}, {tempROI}];
                    addlistener(tempROI,'ROIClicked',@obj.rectClck);
                    addlistener(tempROI,'DeletingROI',@obj.rectDlt);
                    addlistener(tempROI,'ROIMoved',@obj.rectMv);
                case(obj.ELPS)
                    tempROI = drawellipse(obj.imageAxes);
                    if (isempty(tempROI.Center))
                        return;
                    end
                    tempROI.DrawingArea = "unlimited";
                    tempROI.FaceAlpha = 0;
                    obj.elpsROIs{obj.slice,obj.channel} = [obj.elpsROIs{obj.slice,obj.channel}, {tempROI}];
                    addlistener(tempROI,'ROIClicked',@obj.elpsClck);
                    addlistener(tempROI,'DeletingROI',@obj.elpsDlt);
                    addlistener(tempROI,'ROIMoved',@obj.elpsMv);
                case(obj.LINE)
                    tempROI = drawpolyline(obj.imageAxes);
                    if (isempty(tempROI.Position))
                        return;
                    end
                    tempROI.DrawingArea = "unlimited";
                    obj.lineROIs{obj.slice,obj.channel} = [obj.lineROIs{obj.slice,obj.channel}, {tempROI}];
                    obj.lineWidths{obj.slice,obj.channel} = [obj.lineWidths{obj.slice,obj.channel},1];
                    addlistener(tempROI,'ROIClicked',@obj.lineClck);
                    addlistener(tempROI,'DeletingROI',@obj.lineDlt);
                    addlistener(tempROI,'ROIMoved',@obj.lineMv);
                    addlistener(tempROI,'VertexAdded',@obj.lineMv);
                    addlistener(tempROI,'VertexDeleted',@obj.lineMv);
                    obj.currentROIInd = length(obj.lineROIs{obj.slice,obj.channel});
                case(obj.RING)
                    tempROI = drawellipse(obj.imageAxes);
                    if (isempty(tempROI.Center))
                        return;
                    end
                    tempROI.DrawingArea = "unlimited";
                    tempROI.FaceAlpha = 0;
                    obj.ringROIs{obj.slice,obj.channel} = [obj.ringROIs{obj.slice,obj.channel}, {tempROI}];
                    obj.ringWidths{obj.slice,obj.channel} = [obj.ringWidths{obj.slice,obj.channel}, 1];
                    addlistener(tempROI,'ROIClicked',@obj.ringClck);
                    addlistener(tempROI,'DeletingROI',@obj.ringDlt);
                    addlistener(tempROI,'ROIMoved',@obj.ringMv);
                    obj.currentROIInd = length(obj.ringROIs{obj.slice,obj.channel});
                case(obj.PLGN)
                    tempROI = drawpolygon(obj.imageAxes);
                    if (isempty(tempROI.Position))
                        return;
                    end
                    tempROI.DrawingArea = "unlimited";
                    tempROI.FaceAlpha = 0;
                    obj.plgnROIs{obj.slice,obj.channel} = [obj.plgnROIs{obj.slice,obj.channel}, {tempROI}];
                    addlistener(tempROI,'ROIClicked',@obj.plgnClck);
                    addlistener(tempROI,'DeletingROI',@obj.plgnDlt);
                    addlistener(tempROI,'ROIMoved',@obj.plgnMv);
                    addlistener(tempROI,'VertexAdded',@obj.plgnMv);
                    addlistener(tempROI,'VertexDeleted',@obj.plgnMv);
            end
            tempROI.LineWidth = 1;
            obj.currentROI = tempROI;
            obj.currentROIType = roiType;
            obj.callingApp.setROILabels(obj.ROIParamLabels{obj.currentROIType});
            obj.callingApp.enableROIFields(true);
            obj.setFields();
            obj.buildCompMask(true);
        end

        function deleteCurrentROI(obj)
            tempROI = obj.currentROI;
            switch (obj.currentROIType)
                case (obj.RECT)
                    obj.rectDlt(tempROI,[]);
                case (obj.CIRC)
                    obj.circDlt(tempROI,[]);
                case (obj.ELPS)
                    obj.elpsDlt(tempROI,[]);
                case (obj.RING)
                    obj.ringDlt(tempROI,[]);
                case (obj.LINE)
                    obj.lineDlt(tempROI,[]);
                case (obj.PLGN)
                    obj.plgnDlt(tempROI,[]);
            end
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
                obj.chngActiveROIs(obj.slice,channel);
            end
        end

        function setSlice(obj,slice)
            if (obj.slice~=slice)
                obj.chngActiveROIs(slice,obj.channel);
            end
        end

        function deselectROI(obj)
            obj.currentROI = [];
            obj.currentROIType = obj.NONE;
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
            obj.buildCompMask(true);
        end

        function hit = imgClk(obj,x,y)
            hit=obj.roiMask(y,x);
            if (~hit)
                obj.currentROIType = obj.NONE;
                obj.currentROI = [];
                obj.callingApp.enableROIFields(false);
            end
        end

        function shiftMask(obj,dspl)
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
                    obj.buildCompMask(true);
                    obj.callingApp.updateDisplay();
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
                    obj.buildCompMask(true);
                    obj.callingApp.updateDisplay();
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
                    obj.buildCompMask(true);
                    obj.callingApp.updateDisplay();
                    obj.callingApp.setROIFields(center);
                case {obj.LINE, obj.PLGN}
                    pos = obj.currentROI.Position;
                    center = mean(pos);
                    bounds = [0, 0, 0, 0];
                    if (~obj.checkInbounds(dspl,center,bounds))
                        return;
                    end
                    obj.currentROI.Position = pos + repmat(dspl,size(pos,1),1);
                    obj.buildCompMask(true);
                    obj.callingApp.updateDisplay();
            end
        end

        function duplicateMask(obj,ch,sl1,sl2)
            if (obj.currentROIType~=obj.NONE)
                for sl=sl1:sl2
                    if ((ch==obj.channel) && (sl==obj.slice))
                        continue;
                    end

                    tempROI = copyobj(obj.currentROI,obj.imageAxes);
                    switch(obj.currentROIType)
                        case(obj.CIRC)
                            obj.circROIs{sl,ch} = [obj.circROIs{sl,ch}, {tempROI}];
                            addlistener(tempROI,'ROIClicked',@obj.circClck);
                            addlistener(tempROI,'DeletingROI',@obj.circDlt);
                            addlistener(tempROI,'ROIMoved',@obj.circMv);
                            obj.shiftROI(tempROI,2*obj.imageDim);
                        case(obj.RECT)
                            obj.rectROIs{sl,ch} = [obj.rectROIs{sl,ch}, {tempROI}];
                            addlistener(tempROI,'ROIClicked',@obj.rectClck);
                            addlistener(tempROI,'DeletingROI',@obj.rectDlt);
                            addlistener(tempROI,'ROIMoved',@obj.rectMv);
                            coord = tempROI.Position;
                            coord(1) = coord(1) + 2*obj.imageDim(1);
                            coord(2) = coord(2) + 2*obj.imageDim(2);
                            tempROI.Position = coord;
                        case(obj.ELPS)
                            obj.elpsROIs{sl,ch} = [obj.elpsROIs{sl,ch}, {tempROI}];
                            addlistener(tempROI,'ROIClicked',@obj.elpsClck);
                            addlistener(tempROI,'DeletingROI',@obj.elpsDlt);
                            addlistener(tempROI,'ROIMoved',@obj.elpsMv);
                            obj.shiftROI(tempROI,2*obj.imageDim);
                        case(obj.LINE)
                            obj.lineROIs{sl,ch} = [obj.lineROIs{sl,ch}, {tempROI}];
                            obj.lineWidths{sl,ch} = [obj.lineWidths{sl,ch},obj.lineWidths{obj.slice,obj.channel}(obj.currentROIInd)];
                            addlistener(tempROI,'ROIClicked',@obj.lineClck);
                            addlistener(tempROI,'DeletingROI',@obj.lineDlt);
                            addlistener(tempROI,'ROIMoved',@obj.lineMv);
                            addlistener(tempROI,'VertexAdded',@obj.lineMv);
                            addlistener(tempROI,'VertexDeleted',@obj.lineMv);
                            coord = tempROI.Position;
                            coord = coord + repmat((2*obj.imageDim),size(coord,1),1);
                            tempROI.Position = coord;
                        case(obj.RING)
                            obj.ringROIs{sl,ch} = [obj.ringROIs{sl,ch}, {tempROI}];
                            obj.ringWidths{sl,ch} = [obj.ringWidths{sl,ch}, obj.ringWidths{obj.slice,obj.channel}(obj.currentROIInd)];
                            addlistener(tempROI,'ROIClicked',@obj.ringClck);
                            addlistener(tempROI,'DeletingROI',@obj.ringDlt);
                            addlistener(tempROI,'ROIMoved',@obj.ringMv);
                            obj.shiftROI(tempROI,2*obj.imageDim);
                        case(obj.PLGN)
                            obj.plgnROIs{sl,ch} = [obj.plgnROIs{sl,ch}, {tempROI}];
                            addlistener(tempROI,'ROIClicked',@obj.plgnClck);
                            addlistener(tempROI,'DeletingROI',@obj.plgnDlt);
                            addlistener(tempROI,'ROIMoved',@obj.plgnMv);
                            addlistener(tempROI,'VertexAdded',@obj.plgnMv);
                            addlistener(tempROI,'VertexDeleted',@obj.plgnMv);
                            coord = tempROI.Position;
                            coord = coord + repmat((2*obj.imageDim),size(coord,1),1);
                            tempROI.Position = coord;
                    end
                end

                tempCh = obj.channel;
                tempSl = obj.slice;
                obj.channel = ch;
                for sl=sl1:sl2
                    if (tempCh~=ch || tempSl~=sl)
                        obj.shiftROIs(sl,ch,-2*obj.imageDim);
                        obj.slice = sl;
                        obj.buildCompMask(false);
                        obj.shiftROIs(sl,ch,2*obj.imageDim);
                    end
                end
                obj.channel = tempCh;
                obj.slice = tempSl;
            end
        end

        %% Main function for building mask
        function buildCompMask(obj,update)
            obj.roiMask = 0*obj.roiMask;
            obj.drawMasks(obj.rectROIs);
            obj.drawMasks(obj.circROIs);
            obj.drawMasks(obj.elpsROIs);
            obj.drawMasks(obj.plgnROIs);
            for ii=1:length(obj.lineROIs{obj.slice,obj.channel})
                tempMask = createMask(obj.lineROIs{obj.slice,obj.channel}{ii});
                tempMask = imdilate(tempMask,strel('disk',ceil(obj.lineWidths{obj.slice,obj.channel}(ii)/2)));
                obj.roiMask = obj.roiMask + tempMask;
            end

            for ii=1:length(obj.ringROIs{obj.slice,obj.channel})
                ax = obj.ringROIs{obj.slice,obj.channel}{ii}.SemiAxes;
                cent = obj.ringROIs{obj.slice,obj.channel}{ii}.Center;
                x = (obj.xMap-cent(1));
                y = (obj.yMap-cent(2));
                theta = -pi/180*obj.ringROIs{obj.slice,obj.channel}{ii}.RotationAngle;
                tempMask = (1 < obj.ellipseMask(x,y,(ax(1)-obj.ringWidths{obj.slice,obj.channel}(ii)/2),(ax(2)-obj.ringWidths{obj.slice,obj.channel}(ii)/2),theta))...
                    .* (1 > obj.ellipseMask(x,y,(ax(1)+obj.ringWidths{obj.slice,obj.channel}(ii)/2),(ax(2)+obj.ringWidths{obj.slice,obj.channel}(ii)/2),theta));
                obj.roiMask = obj.roiMask + tempMask;
            end
            obj.roiMasks(:,:,obj.slice,obj.channel) = obj.roiMask > 0;
            if (update)
                obj.callingApp.updateOcc();
            end
        end

        %% functions for displaying mask/occupancy information
        function mask = getMask(obj)
            mask = obj.roiMasks(:,:,obj.slice,obj.channel);
        end

        function mask = getChMask(obj,sl,ch)
            mask = obj.roiMasks(:,:,sl,ch);
        end

        %% File IO
        function success = loadROIs(obj,data)
            success = false;
            obj.writeDim1=0;
            obj.writeDim2=0;
            if (isfield(data,'sepCh'))
                sepCh = data.sepCh;
                dup = (~sepCh) && obj.splitChannel;
                if (isfield(data,'circ'))
                    obj.circROIs = obj.copyROIs(data.circ,dup);
                    obj.attachListeners(obj.circROIs,@obj.circMv,@obj.circClck,@obj.circDlt);
                end

                if (isfield(data,'rect'))
                    obj.rectROIs = obj.copyROIs(data.rect,dup);
                    obj.attachListeners(obj.rectROIs,@obj.rectMv,@obj.rectClck,@obj.rectDlt);
                end

                if (isfield(data,'elps'))
                    obj.elpsROIs = obj.copyROIs(data.elps,dup);
                    obj.attachListeners(obj.elpsROIs,@obj.elpsMv,@obj.elpsClck,@obj.elpsDlt);
                end

                if (isfield(data,'ring'))
                    obj.ringROIs = obj.copyROIs(data.ring,dup);
                    obj.attachListeners(obj.ringROIs,@obj.ringMv,@obj.ringClck,@obj.ringDlt);
                    if (isfield(data,'rws'))
                        obj.ringWidths = obj.copyWidths(data.rws,dup);
                    else
                        return;
                    end
                end

                if (isfield(data,'line'))
                    obj.lineROIs = obj.copyROIs(data.line,dup);
                    for sl=1:obj.slices
                        for ch=1:obj.channels
                            for ii=1:length(obj.lineROIs{sl,ch})
                                obj.lineROIs{sl,ch}{ii}.Parent = obj.imageAxes;
                                addlistener(obj.lineROIs{sl,ch}{ii},'ROIClicked',@obj.lineClck);
                                addlistener(obj.lineROIs{sl,ch}{ii},'DeletingROI',@obj.lineDlt);
                                addlistener(obj.lineROIs{sl,ch}{ii},'ROIMoved',@obj.lineMv);
                                addlistener(obj.lineROIs{sl,ch}{ii},'VertexAdded',@obj.lineMv);
                                addlistener(obj.lineROIs{sl,ch}{ii},'VertexDeleted',@obj.lineMv);
                            end
                        end
                    end
                    if (isfield(data,'lws'))
                        obj.lineWidths = obj.copyWidths(data.lws,dup);
                    else
                        return;
                    end
                end

                if (isfield(data,'plgn'))
                    obj.plgnROIs = obj.copyROIs(data.plgn,dup);
                    for sl=1:obj.slices
                        for ch=1:obj.channels
                            for ii=1:length(obj.plgnROIs{sl,ch})
                                obj.plgnROIs{sl,ch}{ii}.Parent = obj.imageAxes;
                                addlistener(obj.plgnROIs{sl,ch}{ii},'ROIClicked',@obj.plgnClck);
                                addlistener(obj.plgnROIs{sl,ch}{ii},'DeletingROI',@obj.plgnDlt);
                                addlistener(obj.plgnROIs{sl,ch}{ii},'ROIMoved',@obj.plgnMv);
                                addlistener(obj.plgnROIs{sl,ch}{ii},'VertexAdded',@obj.plgnMv);
                                addlistener(obj.plgnROIs{sl,ch}{ii},'VertexDeleted',@obj.plgnMv);
                            end
                        end
                    end
                end

                tempCh = obj.channel;
                tempSl = obj.slice;
                for sl=1:obj.slices
                    for ch=1:obj.channels
                        obj.channel = ch;
                        obj.slice = sl;
                        obj.buildCompMask(false);
                        if (tempCh~=ch || tempSl~=sl)
                            obj.shiftROIs(sl,ch,2*obj.imageDim);
                        end
                    end
                end
                obj.channel = tempCh;
                obj.slice = tempSl;
            else
                return;
            end
            success = true;
        end

        function saveROIs(obj,fname)
            for sl=1:obj.slices
                for ch=1:obj.channels
                    if (obj.channel~=ch || obj.slice~=sl)
                        obj.shiftROIs(sl,ch,-2*obj.imageDim);
                    end
                end
            end

            sepCh = obj.splitChannel;
            circ = obj.circROIs;
            rect = obj.rectROIs;
            elps = obj.elpsROIs;
            line = obj.lineROIs;
            ring = obj.ringROIs;
            plgn = obj.plgnROIs;
            lws = obj.lineWidths;
            rws = obj.ringWidths;
            save(fname,"sepCh","circ","rect","elps","line","ring","plgn","lws","rws");

            for sl=1:obj.slices
                for ch=1:obj.channels
                    if (obj.channel~=ch || obj.slice~=sl)
                        obj.shiftROIs(sl,ch,2*obj.imageDim);
                    end
                end
            end
        end

        function status = printROIMasks(obj,basename,path)
            obj.numRois = zeros(obj.slices,obj.channels);
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
                
                ind = 0;
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

                    if ((sl~=obj.slice) || (ch~=obj.channel))
                        obj.shiftROIs(sl,ch,-2*obj.imageDim);
                    end

                    ind = obj.printMaskArr(obj.rectROIs{sl,ch},basename,slPath,ind);
                    ind = obj.printMaskArr(obj.circROIs{sl,ch},basename,slPath,ind);
                    ind = obj.printMaskArr(obj.elpsROIs{sl,ch},basename,slPath,ind);
                    ind = obj.printMaskArr(obj.plgnROIs{sl,ch},basename,slPath,ind);

                    for ii=1:length(obj.ringROIs{sl,ch})
                        ind = ind+1;
                        ax = obj.ringROIs{sl,ch}{ii}.SemiAxes;
                        cent = obj.ringROIs{sl,ch}{ii}.Center;
                        x = (obj.xMap-cent(1));
                        y = (obj.yMap-cent(2));
                        theta = -pi/180*obj.ringROIs{sl,ch}{ii}.RotationAngle;
                        tempMask = (1 < obj.ellipseMask(x,y,(ax(1)-obj.ringWidths{sl,ch}(ii)/2),(ax(2)-obj.ringWidths{sl,ch}(ii)/2),theta))...
                            .* (1 > obj.ellipseMask(x,y,(ax(1)+obj.ringWidths{sl,ch}(ii)/2),(ax(2)+obj.ringWidths{sl,ch}(ii)/2),theta));
                        ff = fullfile(slPath,strcat(basename,'_',num2str(ind),'.tif'));
                        imwrite((tempMask+0),ff, 'Compression','none');
                    end

                    for ii=1:length(obj.lineROIs{sl,ch})
                        ind = ind+1;
                        tempMask = createMask(obj.lineROIs{obj.slice,obj.channel}{ii});
                        tempMask = imdilate(tempMask,strel('disk',ceil(obj.lineWidths{obj.slice,obj.channel}(ii)/2))) > 0;
                        ff = fullfile(slPath,strcat(basename,'_',num2str(ind),'.tif'));
                        imwrite((tempMask+0),ff, 'Compression','none');
                    end

                    if ((sl~=obj.slice) || (ch~=obj.channel))
                        obj.shiftROIs(sl,ch,2*obj.imageDim);
                    end
                    obj.numRois(sl,ch) = ind;
                end
            end
            status = true;
        end

        %% Listeners
        % helper functions used for various listeners
        function chngROI(obj,roi,type)
            if (isempty(obj.currentROI)||~eq(obj.currentROI,roi))
                if (type==obj.LINE)
                    for ii=1:length(obj.lineROIs{obj.slice,obj.channel})
                        if (eq(roi,obj.lineROIs{obj.slice,obj.channel}{ii}))
                            obj.currentROIInd = ii;
                            break;
                        end
                    end
                elseif (type==obj.RING)
                    for ii=1:length(obj.ringROIs{obj.slice,obj.channel})
                        if (eq(roi,obj.ringROIs{obj.slice,obj.channel}{ii}))
                            obj.currentROIInd = ii;
                            break;
                        end
                    end
                end
                obj.currentROIType = type;
                obj.currentROI = roi;
                obj.callingApp.setROILabels(obj.ROIParamLabels{obj.currentROIType});
                obj.callingApp.enableROIFields(true);
            end
            obj.setFields();
            obj.callingApp.updateDisplay();
        end

        % Circle ROI Listeners
        function circMv(obj,src,~)
            obj.buildCompMask();
            obj.chngROI(src,obj.CIRC);
        end

        function circClck(obj,src,~)
            obj.callingApp.deselectThresh();
            obj.chngROI(src,obj.CIRC);
        end
        
        function circDlt(obj,src,~)
            for ii=1:length(obj.circROIs{obj.slice,obj.channel})
                if (eq(src,obj.circROIs{obj.slice,obj.channel}{ii}))
                    obj.circROIs{obj.slice,obj.channel}{ii} = [];
                    break;
                end
            end
            obj.currentROI = [];
            obj.circROIs{obj.slice,obj.channel} = obj.circROIs{obj.slice,obj.channel}(~cellfun(@isempty,obj.circROIs{obj.slice,obj.channel}));
            obj.currentROIType = obj.NONE;
            obj.callingApp.enableROIFields(false);
            obj.buildCompMask(true);
            obj.callingApp.updateDisplay();
        end

        % Rectangle ROI Listeners
        function rectMv(obj,src,~)
            obj.buildCompMask(true);
            obj.chngROI(src,obj.RECT);
        end
        
        function rectClck(obj,src,~)
            obj.callingApp.deselectThresh();
            obj.chngROI(src,obj.RECT);
        end

        function rectDlt(obj,src,~)
            for ii=1:length(obj.rectROIs{obj.slice,obj.channel})
                if(eq(src,obj.rectROIs{obj.slice,obj.channel}{ii}))
                    obj.rectROIs{obj.slice,obj.channel}{ii} = [];
                    break;
                end
            end
            obj.currentROI = [];
            obj.rectROIs{obj.slice,obj.channel} = obj.rectROIs{obj.slice,obj.channel}(~cellfun(@isempty,obj.rectROIs{obj.slice,obj.channel}));
            obj.currentROIType = obj.NONE;
            obj.callingApp.enableROIFields(false);
            obj.buildCompMask(true);
            obj.callingApp.updateDisplay();
        end

        % Ellipse ROI Listeners
        function elpsMv(obj,src,~)
            obj.buildCompMask(true);
            obj.chngROI(src,obj.ELPS);
        end

        function elpsClck(obj,src,~)
            obj.callingApp.deselectThresh();
            obj.chngROI(src,obj.ELPS);
        end

        function elpsDlt(obj,src,~)
            for ii=1:length(obj.elpsROIs{obj.slice,obj.channel})
                if (eq(src,obj.elpsROIs{obj.slice,obj.channel}{ii}))
                    obj.elpsROIs{obj.slice,obj.channel}{ii} = [];
                    break;
                end
            end
            obj.currentROI = [];
            obj.elpsROIs{obj.slice,obj.channel} = obj.elpsROIs{obj.slice,obj.channel}(~cellfun(@isempty,obj.elpsROIs{obj.slice,obj.channel}));
            obj.currentROIType = obj.NONE;
            obj.callingApp.enableROIFields(false);
            obj.buildCompMask(true);
            obj.callingApp.updateDisplay();
        end

        % Line ROI Listeners
        function lineMv(obj,src,~)
            obj.buildCompMask(true);
            obj.chngROI(src,obj.LINE);
        end

        function lineClck(obj,src,~)
            obj.callingApp.deselectThresh();
            obj.chngROI(src,obj.LINE);
        end

        function lineDlt(obj,src,~)
            for ii=1:length(obj.lineROIs{obj.slice,obj.channel})
                if (eq(src,obj.lineROIs{obj.slice,obj.channel}{ii}))
                    obj.lineROIs{obj.slice,obj.channel}{ii} = [];
                    break;
                end
            end
            obj.currentROI = [];
            keep = ~cellfun(@isempty,obj.lineROIs{obj.slice,obj.channel});
            obj.lineROIs{obj.slice,obj.channel} = obj.lineROIs{obj.slice,obj.channel}(keep);
            obj.lineWidths{obj.slice,obj.channel} = obj.lineWidths{obj.slice,obj.channel}(keep);
            obj.currentROIType = obj.NONE;
            obj.callingApp.enableROIFields(false);
            obj.buildCompMask(true);
            obj.callingApp.updateDisplay();
        end

        % Ring ROI Listeners
        function ringMv(obj,src,~)
            obj.buildCompMask(true);
            obj.chngROI(src,obj.RING);
        end

        function ringClck(obj,src,~)
            obj.callingApp.deselectThresh();
            obj.chngROI(src,obj.RING);
        end

        function ringDlt(obj,src,~)
            for ii=1:length(obj.ringROIs{obj.slice,obj.channel})
                if (eq(src,obj.ringROIs{obj.slice,obj.channel}{ii}))
                    obj.ringROIs{obj.slice,obj.channel}{ii} = [];
                    break;
                end
            end
            obj.currentROI = [];
            keep = ~cellfun(@isempty,obj.ringROIs{obj.slice,obj.channel});
            obj.ringROIs{obj.slice,obj.channel} = obj.ringROIs{obj.slice,obj.channel}(keep);
            obj.ringWidths{obj.slice,obj.channel} = obj.ringWidths{obj.slice,obj.channel}(keep);
            obj.currentROIType = obj.NONE;
            obj.callingApp.enableROIFields(false);
            obj.buildCompMask(true);
            obj.callingApp.updateDisplay();
        end

        % Polygon ROI Listeners
        function plgnMv(obj,src,~)
            obj.buildCompMask(true);
            obj.chngROI(src,obj.PLGN);
        end

        function plgnClck(obj,src,~)
            obj.callingApp.deselectThresh();
            obj.chngROI(src,obj.PLGN);
        end

        function plgnDlt(obj,src,~)
            for ii=1:length(obj.plgnROIs{obj.slice,obj.channel})
                if (eq(src,obj.plgnROIs{obj.slice,obj.channel}{ii}))
                    obj.plgnROIs{obj.slice,obj.channel}{ii} = [];
                    break;
                end
            end
            obj.currentROI = [];
            obj.plgnROIs{obj.slice,obj.channel} = obj.plgnROIs{obj.slice,obj.channel}(~cellfun(@isempty,obj.plgnROIs{obj.slice,obj.channel}));
            obj.currentROIType = obj.NONE;
            obj.callingApp.enableROIFields(false);
            obj.buildCompMask(true);
            obj.callingApp.updateDisplay();
        end

        %% other helper functions
        function shiftROIs(obj,slice,channel,shiftVec)
            for ii=1:length(obj.rectROIs{slice,channel})
                coord = obj.rectROIs{slice,channel}{ii}.Position;
                coord(1) = coord(1) + shiftVec(1);
                coord(2) = coord(2) + shiftVec(2);
                obj.rectROIs{slice,channel}{ii}.Position = coord;
            end

            for ii=1:length(obj.circROIs{slice,channel})
                obj.shiftROI(obj.circROIs{slice,channel}{ii},shiftVec);
            end

            for ii=1:length(obj.elpsROIs{slice,channel})
                obj.shiftROI(obj.elpsROIs{slice,channel}{ii},shiftVec);
            end

            for ii=1:length(obj.ringROIs{slice,channel})
                obj.shiftROI(obj.ringROIs{slice,channel}{ii},shiftVec);
            end

            for ii=1:length(obj.lineROIs{slice,channel})
                coord = obj.lineROIs{slice,channel}{ii}.Position;
                coord = coord + repmat(shiftVec,size(coord,1),1);
                obj.lineROIs{slice,channel}{ii}.Position = coord;
            end

            for ii=1:length(obj.plgnROIs{slice,channel})
                coord = obj.plgnROIs{slice,channel}{ii}.Position;
                coord = coord + repmat(shiftVec,size(coord,1),1);
                obj.plgnROIs{slice,channel}{ii}.Position = coord;
            end
        end

        function chngActiveROIs(obj,slice,channel)
            obj.shiftROIs(obj.slice,obj.channel,2*obj.imageDim);
            obj.slice=slice;
            obj.channel=channel;
            obj.shiftROIs(obj.slice,obj.channel,-2*obj.imageDim);
            obj.currentROI = [];
            obj.currentROIType = obj.NONE;
            obj.callingApp.enableROIFields(false);
        end

        function newArray = copyROIs(obj,arr,dup)
            newArray = cell(obj.slices,obj.channels);
            for sl=1:obj.slices
                for ch=1:obj.channels
                    newArray{sl,ch} = cell(0);
                end
            end

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
    
            if (dup)
                for sl=1:obj.writeDim1
                    for ch=1:obj.channels
                        newArray{sl,ch} = arr{sl,1};
                    end
                end
            else
                for sl=1:obj.writeDim1
                    for ch=1:obj.writeDim2
                        newArray{sl,ch} = arr{sl,ch};
                    end
                end
            end
        end

        function newArray = copyWidths(obj,arr,dup)
            newArray = cell(obj.slices,obj.channels);
            for sl=1:obj.slices
                for ch=1:obj.channels
                    newArray{sl,ch} = [];
                end
            end

            if (dup)
                for sl=1:obj.writeDim1
                    for ch=1:obj.channels
                        newArray{sl,ch} = arr{sl,1};
                    end
                end
            else
                for sl=1:obj.writeDim1
                    for ch=1:obj.writeDim2
                        newArray{sl,ch} = arr{sl,ch};
                    end
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

        function attachListeners(obj,rois,mvFun,clckFun,delFun)
            for sl=1:obj.slices
                for ch=1:obj.channels
                    for ii=1:length(rois{sl,ch})
                        rois{sl,ch}{ii}.Parent = obj.imageAxes;
                        addlistener(rois{sl,ch}{ii},'ROIClicked',clckFun);
                        addlistener(rois{sl,ch}{ii},'DeletingROI',delFun);
                        addlistener(rois{sl,ch}{ii},'ROIMoved',mvFun);
                    end
                end
            end
        end

        function drawMasks(obj,rois)
            for ii=1:length(rois{obj.slice,obj.channel})
                obj.roiMask = obj.roiMask + createMask(rois{obj.slice,obj.channel}{ii});
            end
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

        function ind = printMaskArr(rois,basename,path,ind)
            for ii=1:length(rois)
                ind = ind + 1;
                tempMask = createMask(rois{ii});
                ff = fullfile(path,strcat(basename,'_',num2str(ind),'.tif'));
                imwrite((tempMask+0),ff, 'Compression','none');
            end
        end
    end
end