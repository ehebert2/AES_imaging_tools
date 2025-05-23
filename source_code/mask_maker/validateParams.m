%%%%%%%%%%%%%%%%%%%
% used to make sure all fields expected in motion struct are present and
% make a degree of sense
%%%%%%%%%%%%%%%%%%%

function params = validateParams(params)
    if (~isfield(params,'channels'))
        params.channels = 1;
    end
    
    if (~isfield(params,'mainChannel'))
        params.mainChannel = 1;
    else
        params.mainChannel = min(params.mainChannel,params.channels);
    end

    if (~isfield(params,'splitChannels'))
        params.splitChannels = false;
    elseif (params.channels == 1)
        params.splitChannels = false;
    end

    if (~isfield(params,'slices'))
        params.slices = 1;
    elseif (params.slices < 1)
        params.slices = 1;
    end

    if (~isfield(params,'reslice'))
        params.reslice = false;
    elseif (params.reslice)
        if (~isfield(params,'resliceMap'))
            params.reslice = false;
        else
            params.slicesOut = max(params.resliceMap(:,:,1),[],'all');
            params.channelsOut = max(params.resliceMap(:,:,2),[],'all');
            if ((params.slicesOut*params.channelsOut)~=(params.slices*params.channels))
                params.reslice = false;
            end
        end
    end

    if (~params.reslice)
        params.slicesOut = params.slices;
        params.channelsOut = params.channels;
        params.resliceMap = [];
    end
    slices2 = params.slicesOut;
    channels2 = params.channelsOut;

    if (~isfield(params,'volume'))
        params.volume = params.slices>1;
    end

    if (~isfield(params,'expMask'))
        params.expMask = [];
    end

    if (~isfield(params,'compress'))
        params.compress = false;
    end

    if (~isfield(params,'roiAes'))
        params.roiAes = cell(slices2,channels2);
        if (isempty(params.expMask))
            params.compress = false;
        end
        params.aesNames = cell(slices2,channels2);
    elseif (~isfield(params,'aesNames'))
        params.aesNames = cell(slices2,channels2);
    end

    if (~isfield(params,'roiSmpl'))
        params.roiSmpl = cell(slices2,channels2);
        params.smplNames = cell(slices2,channels2);
    elseif (~isfield(params,'smplNames'))
        params.smplNames = cell(slices2,channels2);
    end

    if (~params.splitChannels)
        if (size(params.roiAes,2)~=channels2)
            temp = params.roiAes(:,1);
            params.roiAes = cell(slices2,channels2);
            for ch=1:channels2
                params.roiAes(:,ch) = temp;
            end
        end

        if (size(params.roiSmpl,2)~=channels2)
            temp = params.roiSmpl(:,1);
            params.roiSmpl = cell(slices2,channels2);
            for ch=1:channels2
                params.roiSmpl(:,ch) = temp;
            end
        end
    end

    if ~isfield(params,'numRoiAes')
        params.numRoiAes = zeros(slices2,channels2);
        for sl=1:slices2
            for ch=1:channels2
                if (~isempty(params.roiAes{sl,ch}))
                    params.numRoiAes(sl,ch) = size(params.roiAes{sl,ch},3);
                end
            end
        end
    end

    if ~isfield(params,'numRoiSmpl')
        params.numRoiSmpl = zeros(slices2,channels2);
        for sl=1:slices2
            for ch=1:channels2
                if (~isempty(params.roiSmpl{sl,ch}))
                    params.numRoiSmpl(sl,ch) = size(params.roiSmpl{sl,ch},3);
                end
            end
        end
    end

    if (~isfield(params,'saveMasks'))
        params.saveMasks = false;
    end

    if (isempty(params.roiAes{1,1}))
        if (isempty(params.roiSmpl{1,1}))
            if (isempty(params.expMask))
                if (~isfield(params,'dim'))
                    params.dim = [];
                end
                params.saveMasks = false;
            else
                params.dim = [size(params.expMask,1),size(params.expMask,2)];
            end
        else
            params.dim = [size(params.roiSmpl{1,1},1),size(params.roiSmpl{1,1},2)];
        end
    else
        params.dim = [size(params.roiAes{1,1},1),size(params.roiAes{1,1},2)];
    end

    if (~isfield(params,'saveVid'))
        params.saveVid = false;
    end

    if (~isfield(params,'zeroVid'))
        params.zeroVid = false;
    end

    if (isfield(params,'aesTrace'))
        if (isempty(params.roiAes))
            params.aesTrace = false;
        end
    else
        params.aesTrace = false;
    end

    if (isfield(params,'smplTrace'))
        if (isempty(params.roiSmpl))
            params.smplTrace = false;
        end
    else
        params.smplTrace = false;
    end

    if (~isfield(params,'zeroTrace'))
        params.zeroTrace = false;
    end

    if (~isfield(params,'mtn'))
        params.mtn = false;
    end

    if (~params.mtn)
        params.mtnIntlWndw = 1;
        params.mtnWndw = 3;
        params.sptlWndw = 0;
        params.passes = 1;
        params.intlPasses = 1;
        params.fft = false;
        params.mtnOverlap = false;
    else
        if (~isfield(params,'mtnIntlWndw'))
            params.mtnIntlWndw = 1;
        end
    
        if (~isfield(params,'mtnWndw'))
            params.mtnWndw = 3;
        end
    
        if (~isfield(params,'sptlWndw'))
            params.sptlWndw = 0;
        end
    
        if (~isfield(params,'passes'))
            params.passes = 1;
        end
    
        if(~isfield(params,'intlPasses'))
            params.intlPasses = 1;
        end
    
        if (~isfield(params,'fft'))
            params.fft = false;
        end

        if (~isfield(params,'mtnOverlap'))
            params.mtnOverlap = false;
        end
    end    

    if (~isfield(params,'bgTrace'))
        params.bgTrace = false;
    end

    if (~isfield(params,'yBidirectional'))
        params.yBidirectional = false;
        params.bidiShift = 0;
    else
        if (~isfield(params,'bidiShift'))
            params.bidiShift = 0;
        end
    end

    if (~isfield(params,'zBidirectional'))
        params.zBidirectional = false;
    elseif (~params.volume)
        params.zBidirectional = false;
    end

    params.bidirectional = params.zBidirectional || params.yBidirectional;

    if (~isfield(params,'multiVid'))
        params.multiVid = false;
    end

    if (~isfield(params,'proj'))
        params.proj = false;
    end
end