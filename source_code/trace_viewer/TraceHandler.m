%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Object for handling and plotting trace data after parsing with parsing object
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef TraceHandler < handle
    properties (GetAccess = public)
        % constant data
        aesTracesRef
        smplTracesRef
        bgTracesRef
        overlapTraces
        mtnTraces
        numFiles
        aesNumPx
        smplNumPx
        traceZeroed
        
        % booleans telling if traces are available
        hasAes
        hasSmpl
        hasBg
        hasMtn
        hasOverlap
        
        % filter variables
        filter
        filtering
        filterFwhm

        % processed data
        aesTraces
        smplTraces
        aesTotTraces
        smplTotTraces
        bgTraces
        
        % additional ui variables
        overlapEdges
        yLimAes
        yLimSmpl
        numAesTraces
        numSmplTraces

        % ui variables
        fps
        x
        frames
        sbtrBG
        shOver
        photonCal
        usingFnc
        usrFnc

        % information about file and video structure
        status
        splitChannel
        channels
        channel
        fIndex

        % names
        aesNames
        smplNames
    end

    % Trace types
    properties (Constant)
        AES = 1;
        SAMPLE = 2;
        BG = 3;
        MOTION = 4;
        NONE = 5;
    end

    methods
        function obj = TraceHandler(outputParser)
            obj.frames = outputParser.frames;
            obj.numFiles = length(obj.frames);
            obj.aesTracesRef = outputParser.aesTraces;
            obj.smplTracesRef = outputParser.smplTraces;
            obj.bgTracesRef = outputParser.bgTraces;
            obj.overlapTraces = outputParser.overlapTraces;
            obj.mtnTraces = outputParser.mtnTraces;
            obj.splitChannel = outputParser.splitChannels;
            obj.channels = outputParser.channels;
            obj.traceZeroed = outputParser.traceZeroed;
            obj.channel = 1;
            obj.fps = 1;
            obj.fIndex = 1;
            obj.photonCal = 0;
            obj.usingFnc = false;
            obj.usrFnc = [];
            obj.x = (1:obj.frames(obj.fIndex))'/obj.fps;
            obj.filter = [];
            obj.sbtrBG = ~isempty(obj.bgTracesRef);
            obj.aesNumPx = outputParser.aesNumPx;
            obj.smplNumPx = outputParser.smplNumPx;
            obj.hasAes = ~isempty(obj.aesTracesRef);
            obj.hasSmpl = ~isempty(obj.smplTracesRef);
            obj.hasBg = ~isempty(obj.bgTracesRef);
            obj.hasMtn = ~isempty(obj.mtnTraces);
            obj.hasOverlap = ~isempty(obj.overlapTraces);

            obj.yLimAes = cell(obj.channels,1);
            obj.yLimSmpl = cell(obj.channels,1);
            obj.numAesTraces = zeros(obj.channels,1);
            obj.numSmplTraces = zeros(obj.channels,1);
            obj.aesTraces = cell(obj.channels,1);
            obj.smplTraces = cell(obj.channels,1);
            obj.bgTraces = cell(obj.channels,1);
            for ch=1:obj.channels
                if (obj.hasAes)
                    obj.numAesTraces(ch) = length(obj.aesTracesRef{1}{ch});
                    obj.yLimAes{ch} = zeros(obj.numAesTraces(ch),2);
                    obj.aesTraces{ch} = cell(obj.numAesTraces(ch),1);
                end
                
                if (obj.hasSmpl)
                    obj.numSmplTraces(ch) = length(obj.smplTracesRef{1}{ch});
                    obj.yLimSmpl{ch} = zeros(obj.numSmplTraces(ch),2);
                    obj.smplTraces{ch} = cell(obj.numSmplTraces(ch),1);
                end
            end
            obj.aesTotTraces = zeros(obj.frames(obj.fIndex),obj.channels);
            obj.smplTotTraces = zeros(obj.frames(obj.fIndex),obj.channels);
            
            % getting bounds of rectangles to draw denoting out of bounds
            % periods of traces
            obj.shOver = ~isempty(obj.overlapTraces);
            if (obj.shOver)
                obj.overlapEdges = cell(obj.numFiles,1);
                for fl=1:obj.numFiles
                    obj.overlapEdges{fl} = cell(obj.channels,1);
                    if (obj.splitChannel)
                        for ch=1:obj.channels
                            if (obj.numSmplTraces(ch)>0)
                                obj.overlapEdges{fl}{ch} = cell(obj.numSmplTraces(ch),1);
                                for ii=1:obj.numSmplTraces(ch)
                                    obj.overlapEdges{fl}{ch}{ii} = obj.findEdges(obj.overlapTraces{fl}{ch}(:,ii));
                                end
                            end
                        end
                    else
                        obj.overlapEdges{fl} = cell(obj.channels,1);
                        obj.overlapEdges{fl}{1} = cell(obj.numSmplTraces(1),1);
                        for ii=1:obj.numSmplTraces(1)
                            obj.overlapEdges{fl}{1}{ii} = obj.findEdges(obj.overlapTraces{fl}{1}(:,ii));
                        end
                        
                        for ch=2:obj.channels
                            obj.overlapEdges{fl}{ch} = obj.overlapEdges{fl}{1};
                        end
                    end
                end
            end

            obj.aesNames = outputParser.aesNames;
            obj.smplNames = outputParser.smplNames;
            obj.buildTraces();
        end

        %% functions for getting current traces (after processing)
        function trace = getAesTrace(obj,index)
            trace = obj.aesTraces{obj.channel}{index};
        end

        function trace = getSmplTrace(obj,index)
            trace = obj.smplTraces{obj.channel}{index};
        end

        function trace = getBgTrace(obj)
            trace = obj.bgTraces{ch};
        end

        function trace = getOverlapTrace(obj,index)
            trace = obj.overlapTraces{obj.fIndex}{ch}(:,index);
        end

        function trace = getMtnTrace(obj)
            trace = obj.mtnTraces{obj.fIndex};
        end

        function frames=getFrames(obj)
            frames=obj.frames(obj.fIndex);
        end

        %% setter functions
        function setChannel(obj,channel)
            obj.channel = channel;
        end

        function setFile(obj,fIndex)
            obj.fIndex = fIndex;
            obj.aesTotTraces = zeros(obj.frames(obj.fIndex),obj.channels);
            obj.smplTotTraces = zeros(obj.frames(obj.fIndex),obj.channels);
            obj.x = (1:obj.frames(obj.fIndex))'/obj.fps;
            obj.setFilter(obj.filterFwhm);
            obj.buildTraces();
        end

        %% functions changing how data is plotted
        function setFPS(obj,fps)
            if (fps ~= obj.fps)
                obj.fps = fps;
                obj.x = (1:obj.frames(obj.fIndex))'/fps;
            end
        end

        function setFilter(obj,fwhm)
            if (fwhm ~= 0)
                obj.filterFwhm = fwhm;
                obj.filtering = true;
                l = length(obj.x);
                a = fwhm^2/log(2);
                f = (linspace(-1,1,l)*obj.fps/2)';
                obj.filter = exp(-f.^2/a);
                obj.filter = ifftshift(obj.filter);
            else
                obj.filtering = false;
            end
        end

        function subtractBG(obj,enable)
            obj.sbtrBG = (enable && ~isempty(obj.bgTracesRef));
        end

        function showOverlap(obj,enable)
            obj.shOver = (enable && ~isempty(obj.overlapTraces));
        end

        function setPhotonCal(obj,val)
            obj.photonCal = val;
        end

        % validates and changes the user function (fnc is function handle)
        function status=setUsrFnc(obj,enable,fnc)
            status = false;
            if (obj.hasSmpl)
                if (enable)
                    if (~isempty(fnc))
                        try
                            ch=1;
                            while(obj.numSmplTraces(ch)==0)
                                ch=ch+1;
                            end
                            fnc(obj.x,obj.smplTraces{ch}{1});
                        catch
                            disp('Function execution error');
                            return;
                        end
                        obj.usrFnc = fnc;
                    elseif (isempty(obj.usrFnc))
                        return;
                    end
                    obj.usingFnc=true;
                else
                    obj.usingFnc = false;
                end
                status = true;
            end
        end

        %% process traces
        function buildTraces(obj)
            if (obj.sbtrBG)
                for ch=1:obj.channels
                    if (obj.hasAes)
                        for ii=1:obj.numAesTraces(ch)
                            obj.aesTraces{ch}{ii} = obj.aesTracesRef{obj.fIndex}{ch}{ii}(:,1) - obj.bgTracesRef{obj.fIndex}{ch}(:,1);
                        end
                    end

                    if (obj.hasSmpl)
                        for ii=1:obj.numSmplTraces(ch)
                            obj.smplTraces{ch}{ii} = obj.smplTracesRef{obj.fIndex}{ch}{ii}(:,1) - obj.bgTracesRef{obj.fIndex}{ch}(:,1);
                        end
                    end
                end
            else
                for ch=1:obj.channels
                    if (obj.hasAes)
                        for ii=1:obj.numAesTraces(ch)
                            obj.aesTraces{ch}{ii} = obj.aesTracesRef{obj.fIndex}{ch}{ii}(:,1);
                        end
                    end

                    if (obj.hasSmpl)
                        for ii=1:obj.numSmplTraces(ch)
                            obj.smplTraces{ch}{ii} = obj.smplTracesRef{obj.fIndex}{ch}{ii}(:,1);
                        end
                    end
                end
            end

            if (obj.photonCal ~= 0)
                for ch=1:obj.channels
                    if (obj.hasAes)
                        for ii=1:obj.numAesTraces(ch)
                            obj.aesTraces{ch}{ii} = obj.aesTraces{ch}{ii} * obj.photonCal * obj.aesNumPx{ch}(ii);
                        end
                    end

                    if (obj.hasSmpl)
                        for ii=1:obj.numSmplTraces(ch)
                            obj.smplTraces{ch}{ii} = obj.smplTraces{ch}{ii} * obj.photonCal * obj.smplNumPx{ch}(ii);
                        end
                    end

                    if (obj.hasBg)
                        obj.bgTraces{ch} = obj.bgTracesRef{obj.fIndex}{ch}(:,1) * obj.photonCal;
                    end
                end
            elseif (obj.hasBg)
                for ch=1:obj.channels
                    obj.bgTraces{ch} = obj.bgTracesRef{obj.fIndex}{ch}(:,1);
                end
            end

            if (obj.filtering)
                for ch=1:obj.channels
                    if (obj.hasAes)
                        for ii=1:obj.numAesTraces(ch)
                            obj.aesTraces{ch}{ii} = abs(ifft(fft(obj.aesTraces{ch}{ii}).*obj.filter));
                        end
                    end

                    if (obj.hasSmpl)
                        for ii=1:obj.numSmplTraces(ch)
                            obj.smplTraces{ch}{ii} = abs(ifft(fft(obj.smplTraces{ch}{ii}).*obj.filter));
                        end
                    end

                    if (obj.hasBg)
                        obj.bgTraces{ch} = abs(ifft(fft(obj.bgTraces{ch}).*obj.filter));
                    end
                end
            end

            if (obj.usingFnc)
                for ch=1:obj.channels
                    for ii=1:obj.numSmplTraces(ch)
                        obj.smplTraces{ch}{ii} = obj.usrFnc(obj.x,obj.smplTraces{ch}{ii});
                    end
                end
            end

            obj.aesTotTraces = obj.aesTotTraces * 0;
            obj.smplTotTraces = obj.smplTotTraces*0;
            for ch=1:obj.channels
                if (obj.hasAes)
                    for ii=1:obj.numAesTraces(ch)
                        obj.yLimAes{ch}(ii,2) = max(obj.aesTraces{ch}{ii});
                        if (obj.sbtrBG)
                            obj.yLimAes{ch}(ii,1) = 0;
                        else
                            obj.yLimAes{ch}(ii,1) = min(obj.aesTraces{ch}{ii});
                        end
                        obj.aesTotTraces(:,ch) = obj.aesTotTraces(:,ch)+obj.aesTraces{ch}{ii};
                    end
                end

                if (obj.hasSmpl)
                    for ii=1:obj.numSmplTraces(ch)
                        obj.yLimSmpl{ch}(ii,2) = max(obj.smplTraces{ch}{ii});
                        if (obj.sbtrBG && ~obj.usingFnc)
                            obj.yLimSmpl{ch}(ii,1) = 0;
                        else
                            obj.yLimSmpl{ch}(ii,1) = min(obj.smplTraces{ch}{ii});
                        end
                        obj.smplTotTraces(:,ch) = obj.smplTotTraces(:,ch)+obj.smplTraces{ch}{ii};
                    end
                end
            end
        end

        %% plot traces (status==false means no trace was plotted)
        function status = plotTrace(obj, ax, channel, tracetype, index)
            if (((tracetype==obj.AES)&&~obj.hasAes)||((tracetype==obj.SAMPLE)&&~obj.hasSmpl)||((tracetype==obj.BG)&&~obj.hasBg)||((tracetype==obj.MOTION)&&~obj.hasMtn))
                tracetype=obj.NONE;
            end

            switch (tracetype)
                case (obj.AES)
                    if (index==0)
                        plot(ax,obj.x,obj.aesTotTraces(:,channel));
                        ymin = sum(obj.yLimAes{channel}(:,1));
                        ymax = sum(obj.yLimAes{channel}(:,2));
                        title(ax, 'Combined AES Region');
                    else
                        plot(ax,obj.x,obj.aesTraces{channel}{index});
                        ymin = obj.yLimAes{channel}(index,1);
                        ymax = obj.yLimAes{channel}(index,2);
                        title(ax, obj.aesNames{channel}{index});
                    end
                    if (obj.photonCal==0)
                        ylabel(ax,'Mean Pixel Value');
                    else
                        ylabel(ax,'Photons')
                    end
                    legend(ax,'off');
                case (obj.SAMPLE)
                    if (index==0)
                        plot(ax,obj.x,obj.smplTotTraces(:,channel));
                        ymin = sum(obj.yLimSmpl{channel}(:,1));
                        ymax = sum(obj.yLimSmpl{channel}(:,2));
                        title(ax, 'Combined Sample Region');
                    else
                        plot(ax,obj.x,obj.smplTraces{channel}{index});
                        ymin = obj.yLimSmpl{channel}(index,1);
                        ymax = obj.yLimSmpl{channel}(index,2);
    
                        if (obj.shOver && ~isempty(obj.overlapEdges{obj.fIndex}{channel}{index}))
                            for ii = 1:size(obj.overlapEdges{obj.fIndex}{channel}{index},1)
                                rectangle(ax,'Position',[obj.x(obj.overlapEdges{obj.fIndex}{channel}{index}(ii,1)),ymin,obj.x(obj.overlapEdges{obj.fIndex}{channel}{index}(ii,2)),(ymax-ymin)],'FaceColor',[1,0,0,0.5],'LineStyle','none');
                            end
                        end
                        title(ax, obj.smplNames{channel}{index});
                    end
                    if (obj.photonCal==0)
                        ylabel(ax,'Mean Pixel Value');
                    else
                        ylabel(ax,'Photons')
                    end
                    legend(ax,'off');
                case (obj.BG)
                    plot(ax,obj.x,obj.bgTraces{channel});
                    ymin = min(obj.bgTraces{channel});
                    ymax = max(obj.bgTraces{channel});
                    title(ax, 'Background');
                    if (obj.photonCal==0)
                        ylabel(ax,'Mean Pixel Value');
                    else
                        ylabel(ax, 'Mean Equivalent Photons');
                    end
                    legend(ax,'off');
                case (obj.MOTION)
                    plot(ax,obj.x,obj.mtnTraces{obj.fIndex}(:,1),obj.x,obj.mtnTraces{obj.fIndex}(:,2));
                    title(ax, 'Motion');
                    xlim(ax,[obj.x(1),obj.x(end)]);
                    ymin = min(obj.mtnTraces{obj.fIndex},[],'all');
                    ymax = max(obj.mtnTraces{obj.fIndex},[],'all');
                    legend(ax, 'X', 'Y');
                    ylabel(ax,'Displacement (px)')
                case (obj.NONE)
                    plot(ax,[0,0], [0,0]);
                    title(ax, 'No Trace');
                    ymin=-1;
                    ymax=1;
            end
            xlim(ax,[obj.x(1),obj.x(end)]);
            ylim(ax,[ymin,ymax]);
            xlabel(ax,'Time (s)');
            status = (tracetype~=obj.NONE);
        end
    end

    methods (Static)
        % helper function to find boundaries of rectangles for plotting
        function edges = findEdges(inBounds)
            edgeStart = zeros(length(inBounds),1);
            edgeEnd = zeros(length(inBounds),1);
            edgeEndInd = 0;
            if (inBounds(1)==0)
                edgeStart(1) = 1;
                edgeStartInd = 1;
            else
                edgeStartInd = 0;
            end

            for ii = 2:length(inBounds)
                if (inBounds(ii)>inBounds(ii-1))
                    edgeEndInd = edgeEndInd+1;
                    edgeEnd(edgeEndInd) = (ii-1);
                elseif (inBounds(ii)<inBounds(ii-1))
                    edgeStartInd = edgeStartInd+1;
                    edgeStart(edgeStartInd) = ii;
                end
            end

            if (inBounds(end)==0)
                edgeEndInd = edgeEndInd+1;
                edgeEnd(edgeEndInd) = length(inBounds);
            end

            if (edgeStartInd>0)
                temp = edgeEnd(1:edgeEndInd)-edgeStart(1:edgeStartInd);
                temp = temp + (temp==0);
                edges = [edgeStart(1:edgeStartInd),temp];
            else
                edges = [];
            end
        end
    end
end