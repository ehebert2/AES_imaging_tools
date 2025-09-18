%%%%%%%%%%%%%%%%%%%%%
% For enoding processing settings into json format
%%%%%%%%%%%%%%%%%%%%%
function jString = paramsToJson(vidNames,params)
    % save video names
    for ii=1:length(vidNames)
        temp=split(vidNames{ii},'.');
        vidNames{ii} = temp{1};
    end

    if (params.multiVid)
        ledger.Videos = vidNames;
    else
        ledger.Videos = vidNames(1);
    end

    % information on stack structure
    reslice.reslice = params.reslice;
    if (params.reslice)
        Stack_Properties.channels = params.channelsOut;
        Stack_Properties.slices = params.slicesOut;
        reslice.input_slices = params.slices;
        reslice.input_channels = params.channels;
        reslice.output_slices = params.slicesOut;
        reslice.output_channels = params.channelsOut;
        map = params.resliceMap;
        keys = cell(1,params.slices*params.channels);
        values = cell(1,params.slices*params.channels);
        for sl=1:params.slices
            for ch=1:params.channels
                index = (sl-1)*params.channels+ch;
                keys{index} = strcat('slice_',num2str(sl),',channel_',num2str(ch));
                values{index} = strcat('slice_',num2str(map(sl,ch,1)),',channel_',num2str(map(sl,ch,2)));
            end
        end
        reslice.map = containers.Map(keys,values);
    else
        Stack_Properties.channels = params.channels;
        Stack_Properties.slices = params.slices;
    end
    Stack_Properties.dimensions.x = params.dim(2);
    Stack_Properties.dimensions.y = params.dim(1);
    Stack_Properties.main_channel = params.mainChannel;
    Stack_Properties.seperate_channels = params.splitChannels;
    bidirectional.y = params.yBidirectional;
    bidirectional.z = params.zBidirectional;
    bidirectional.shift = params.bidiShift;
    Stack_Properties.bidirectional = bidirectional;
    Stack_Properties.reslice = reslice;

    % information about ROIs
    ROIs.masks_saved = params.saveMasks;
    num_vals = size(params.aesNames,1)*size(params.aesNames,2);
    aes_keys = cell(1,num_vals);
    aes_values = cell(1,num_vals);
    smpl_keys = cell(1,num_vals);
    smpl_values = cell(1,num_vals);
    for sl=1:size(params.aesNames,1)
        for ch=1:size(params.aesNames,2)
            index = (sl-1)*size(params.aesNames,2)+ch;
            aes_keys{index} = strcat('slice_',num2str(sl),',channel_',num2str(ch));
            smpl_keys{index} = aes_keys{index};
            if (~isempty(params.aesNames{sl,ch}))
                temp = params.aesNames{sl,ch}{1};
                for ii=2:length(params.aesNames{sl,ch})
                    temp = strcat(temp,',',params.aesNames{sl,ch}{ii});
                end
                aes_values{index} = temp;
            end

            if (~isempty(params.smplNames{sl,ch}))
                temp = params.smplNames{sl,ch}{1};
                for ii=2:length(params.smplNames{sl,ch})
                    temp = strcat(temp,',',params.smplNames{sl,ch}{ii});
                end
                smpl_values{index} = temp;
            end
        end
    end
    ROIs.aes_names = containers.Map(aes_keys,aes_values);
    ROIs.smpl_names = containers.Map(smpl_keys,smpl_values);

    % information about saved video
    Video_Properties.saved_video = params.saveVid;
    Video_Properties.saved_projection = params.proj;
    Video_Properties.compressed = params.compress;
    Video_Properties.filled_bg = params.fillBG;
    if (params.zeroVid)
        Video_Properties.background_subtracted = 'electrical background';
    elseif (params.zeroVidFl)
        Video_Properties.background_subtracted = 'fluorescent background';
    else
        Video_Properties.background_subtracted = 'none';
    end
    
    % information about extracted traces
    Trace_Extraction.sample = params.smplTrace;
    Trace_Extraction.aes = params.aesTrace;
    if (params.zeroTrace)
        Trace_Extraction.background_subtracted = 'electrical background';
    elseif (params.zeroTraceFl)
        Trace_Extraction.background_subtracted = 'fluorescent background';
    else
        Trace_Extraction.background_subtracted = 'none';
    end
    Trace_Extraction.electrical_bg = params.bgTrace;
    Trace_Extraction.fluorescent_bg = params.flTrace;
    Trace_Extraction.displacement = params.mtnTrace;
    Trace_Extraction.in_bounds = params.mtnOverlap;
    if (params.checkAcc)
        Trace_Extraction.max_acc = params.maxAcc;
    else
        Trace_Extraction.max_acc = 'none';
    end

    % information about motion correction
    Motion_Correction.motion_corrected = params.mtn;
    if (params.mtn)
        Motion_Correction.fft_xcor = params.fft;
        Motion_Correction.reference_frames = params.mtnIntlWndw;
        Motion_Correction.temporal_average = params.mtnWndw;
        Motion_Correction.spatial_average = params.sptlWndw;
        Motion_Correction.reference_passes = params.intlPasses;
        Motion_Correction.full_passes = params.passes;
    end

    % attach to main struct and encode json
    ledger.Stack_Properties = Stack_Properties;
    ledger.ROIs = ROIs;
    ledger.Video_Properties = Video_Properties;
    ledger.Trace_Extraction = Trace_Extraction;
    ledger.Motion_Correction = Motion_Correction;
    jString = jsonencode(ledger,PrettyPrint=true);    
end