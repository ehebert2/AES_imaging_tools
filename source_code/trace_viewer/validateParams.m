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

    if (~isfield(params,'expMask'))
        params.expMask = [];
    end

    if (~isfield(params,'compress'))
        params.compress = false;
    end

    if (~isfield(params,'roiAes'))
        params.roiAes = cell(params.channels);
        if (isempty(params.expMask))
            params.compress = false;
        end
        params.aesNames = cell(params.channels);
    elseif (~isfield(params,'aesNames'))
        params.aesNames = cell(params.channels);
    end

    if (~isfield(params,'roiSmpl'))
        params.roiSmpl = cell(params.channels);
        params.smplNames = cell(params.channels);
    elseif (~isfield(params,'smplNames'))
        params.smplNames = cell(params.channels);
    end

    if (~params.splitChannels)
        if (length(params.roiAes)~=params.channels)
            temp = params.roiAes{1};
            params.roiAes = cell(params.channels,1);
            for ch=1:params.channels
                params.roiAes{ch} = temp;
            end
        end

        if (length(params.roiSmpl)~=params.channels)
            temp = params.roiSmpl{1};
            params.roiSmpl = cell(params.channels,1);
            for ch=1:params.channels
                params.roiSmpl{ch} = temp;
            end
        end
    end

    if ~isfield(params,'numRoiAes')
        params.numRoiAes = zeros(params.channels,1);
        for ch=1:params.channels
            if (~isempty(params.roiAes{ch}))
                params.numRoiAes(ch) = size(params.roiAes{ch},3);
            end
        end
    end

    if ~isfield(params,'numRoiSmpl')
        params.numRoiSmpl = zeros(params.channels,1);
        for ch=1:params.channels
            if (~isempty(params.roiSmpl{ch}))
                params.numRoiSmpl(ch) = size(params.roiSmpl{ch},3);
            end
        end
    end

    if (~isfield(params,'saveMasks'))
        params.saveMasks = false;
    end

    if (isempty(params.roiAes{1}))
        if (isempty(params.roiSmpl{1}))
            if (isempty(params.expMask))
                if (~isfield(params,'dim'))
                    params.dim = [];
                end
                params.saveMasks = false;
            else
                params.dim = [size(params.expMask,1),size(params.expMask,2)];
            end
        else
            params.dim = [size(params.roiSmpl{1},1),size(params.roiSmpl{1},2)];
        end
    else
        params.dim = [size(params.roiAes{1},1),size(params.roiAes{1},2)];
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
        if (params.splitChannels)
            params.mtnOverlap = (zeros(params.channels,1)>0);
        else
            params.mtnOverlap = false;
        end
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

        if (params.splitChannels)
            if (length(params.mtnOverlap)==1)
                if (params.mtnOverlap)
                    params.mtnOverlap = (zeros(params.channels,1)>0);
                    for ch=1:params.channels
                        if ((params.numRoiAes(ch)>0)&&(params.numRoiSmpl(ch)>0))
                            params.mtnOverlap(ch)=true;
                        end
                    end
                else
                    params.mtnOverlap = (zeros(params.channels,1)>0);
                end
            end
        else
            if (length(params.mtnOverlap)>1)
                params.mtnOverlap = params.mtnOverlap(1);
            end

            if ((params.numRoiAes(1)==0)||(params.numRoiSmpl(1)==0))
                params.mtnOverlap = false;
            end
        end
    end    

    if (~isfield(params,'bgTrace'))
        params.bgTrace = false;
    end

    if (~isfield(params,'bidirectional'))
        params.bidirectional = false;
        params.bidiShift = 0;
    else
        if (~isfield(params,'bidiShift'))
            params.bidiShift = 0;
        end
    end

    if (~isfield(params,'multiVid'))
        params.multiVid = false;
    end

    if (~isfield(params,'proj'))
        params.proj = false;
    end
end