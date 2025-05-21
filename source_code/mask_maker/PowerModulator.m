classdef PowerModulator < handle
    properties (GetAccess = private, SetAccess = private)    
        highPx
        lowPx
        pxVals        
        calMaps
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
        curveTypes
        
        refPower
        trans
        prePower
        v
        z

        mpOrder
        sliceSep
        attenLen
        orientation
        numPx

        calNames
        prePowerNames
    end

    methods
        function obj = PowerModulator(params,ims)
            obj.sepChannels = params.splitChannels && (params.channels>1);
            obj.channels = params.channels;
            obj.mainChannel = params.mainChannel;
            obj.slices = params.slices;
            obj.mpOrder = 2;
            obj.sliceSep = 1;
            obj.z = 0:(obj.slices-1);
            obj.orientation = ones(obj.channels,1);
            obj.attenLen = ones(obj.channels,1);
            obj.bidirectional = params.bidirectional;

            obj.numPx = size(ims,1)*size(ims,2);
            obj.pxVals = zeros(obj.numPx,(obj.channels*obj.sepChannels+1-obj.sepChannels),obj.slices);
            obj.highPx = round(1*obj.numPx);
            obj.lowPx = round(0.9*obj.numPx);
            obj.curveTypes = obj.UNIFORM * ones(obj.channels,1);
            obj.trans = ones(obj.channels,obj.slices);
            obj.v = ones(obj.channels,obj.slices);
            obj.refPower = zeros((obj.channels*obj.sepChannels+1-obj.sepChannels),obj.slices);
            obj.calNames = cell((obj.channels*obj.sepChannels+1-obj.sepChannels),1);
            obj.prePowerNames = cell((obj.channels*obj.sepChannels+1-obj.sepChannels),1);
            obj.prePower = ones((obj.channels*obj.sepChannels+1-obj.sepChannels),obj.slices);
            obj.calMaps = cell((obj.channels*obj.sepChannels+1-obj.sepChannels),1);
            obj.setIms(ims);
        end

        function setIms(obj,ims)
            if (obj.sepChannels)
                for ch=1:obj.channels
                    for sl=1:obj.slices
                        obj.pxVals(:,ch,sl) = sort(reshape(ims(:,:,sl,(ch*(1+obj.bidirectional)-obj.bidirectional)),[],1));
                    end
                end
            else
                for sl=1:obj.slices
                    obj.pxVals(:,1,sl) = sort(reshape(ims(:,:,sl,(obj.mainChannel*(1+obj.bidirectional)-obj.bidirectional)),[],1));
                end
            end

            obj.calcRefPower();
            obj.calcAllTrans();
            obj.calcAllV();
        end

        function calcRefPower(obj)
            if (obj.sepChannels)
                for ch=1:obj.channels
                    for sl = 1:obj.slices
                        obj.refPower(ch,sl) = mean(obj.pxVals(obj.lowPx:obj.highPx,ch,sl));
                    end
                    obj.refPower(ch,:) = obj.refPower(ch,:) / max(obj.refPower(ch,:));
                end
            else
                for sl = 1:obj.slices
                    obj.refPower(1,sl) = mean(obj.pxVals(obj.lowPx:obj.highPx,1,sl));
                end
                obj.refPower(1,:) = obj.refPower(1,:) / max(obj.refPower(1,:));
            end
        end

        function calcTrans(obj,ch)
            ch = ch*obj.sepChannels+1-obj.sepChannels;
            if (obj.curveTypes(ch)==obj.UNIFORM)
                obj.trans(ch,:) = ones(1,obj.slices);
            elseif (obj.curveTypes(ch)==obj.EXP)
                obj.trans(ch,:) = exp(obj.orientation(ch)*obj.z/(obj.attenLen(ch)*obj.mpOrder));
                obj.trans(ch,:) = obj.trans(ch,:)/max(obj.trans(ch,:));
            elseif (obj.curveTypes(ch)==obj.INVPOW)
                obj.trans(ch,:) = obj.prePower(ch,:)./(obj.refPower(ch,:).^(1/obj.mpOrder));
                obj.trans(ch,:) = obj.trans(ch,:)/max(obj.trans(ch,:));
            end

            if (~obj.sepChannels)
                obj.trans = repmat(obj.trans(1,:),obj.channels,1);
            end
        end

        function calcAllTrans(obj)
            if obj.sepChannels
                for ch=1:obj.channels
                    obj.calcTrans(ch);
                end
            else
                obj.calcTrans(1);
            end
        end

        function fitExp(obj,ch)
            X = [ones(obj.slices,1), obj.z'];
            Y = log(obj.refPower(ch,:)./(obj.prePower(ch,:).^2))';
            C = inv(X'*X)*X'*Y;
            if (C(2)~=0)
                obj.orientation(ch) = 1-2*(C(2)>0);
                obj.attenLen(ch) = 1/abs(C(2));
            end
            obj.calcTrans(ch);
            obj.calcV(ch);
        end

        function calcV(obj,ch)
            if (isempty(obj.calNames{ch}))
                obj.v(ch,:) = obj.trans(ch,:);
            else
                for sl=1:obj.slices
                    if (obj.trans(ch,sl) <= obj.calMaps{ch}(1,2))
                        obj.v(ch,sl) = obj.calMaps{ch}(1,1);
                    else
                        ind = 2;
                        while((obj.calMaps{ch}(ind,2)<obj.trans(ch,sl))&&(ind<size(obj.calMaps{ch},1)))
                            ind=ind+1;
                        end

                        obj.v(ch,sl) = obj.calMaps{ch}(ind,1) + (obj.trans(ch,sl)-obj.calMaps{ch}(ind,2)) * (obj.calMaps{ch}(ind,1)-obj.calMaps{ch}(ind-1,1)) / (obj.calMaps{ch}(ind,2)-obj.calMaps{ch}(ind-1,2));
                    end
                end
            end
        end

        function calcAllV(obj)
            if (obj.sepChannels)
                for ch=1:obj.channels
                    obj.calcV(ch);
                end
            else
                obj.calcV(1);
                obj.v = repmat(obj.v(1,:),obj.channels,1);
            end
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

        function save(obj,fname)
            powerModParams.refPower = obj.refPower;
            powerModParams.trans = obj.trans;
            powerModParams.prePower = obj.prePower;
            powerModParams.calMaps = obj.calMaps;
            powerModParams.v = obj.v;
            powerModParams.z = obj.z;

            powerModParams.mpOrder = obj.mpOrder;
            powerModParams.sliceSep = obj.sliceSep;
            powerModParams.attenLen = obj.attenLen;
            powerModParams.orientation = obj.orientation;
            powerModParams.curveTypes = obj.curveTypes;
            powerModParams.highPx = obj.highPx;
            powerModParams.lowPx = obj.lowPx;

            powerModParams.prePowerNames = obj.prePowerNames;
            powerModParams.calNames = obj.calNames;

            save(fname,"powerModParams",'-append');
        end

        function success = load(obj,data)
            success = false;
            if (~isfield(data,'powerModParams'))
                return;
            end

            params = data.powerModParams;
            obj.refPower = params.refPower;
            obj.trans = params.trans;
            obj.prePower = params.prePower;
            obj.calMaps = params.calMaps;
            obj.v = params.v;
            obj.z = params.z;

            obj.mpOrder = params.mpOrder;
            obj.sliceSep = params.sliceSep;
            obj.attenLen = params.attenLen;
            obj.orientation = params.orientation;
            obj.curveTypes = params.curveTypes;
            obj.highPx = params.highPx;
            obj.lowPx = params.lowPx;

            obj.prePowerNames = params.prePowerNames;
            obj.calNames = params.calNames;

            success = true;
        end

        function factor = getImFactor(obj,ch,sl)
            factor = (obj.trans(ch,sl)/obj.prePower(ch,sl))^obj.mpOrder;
        end

        function power = getPower(obj,ch)
            power = obj.refPower(ch,:);
        end

        function power = getPrePower(obj,ch)
            power = obj.prePower(ch,:);
        end

        function t = getTrans(obj,ch)
            t = obj.trans(ch,:);
        end

        function volt = getV(obj,ch)
            volt = obj.v(ch,:);
        end

        function cal = getCal(obj,ch)
            cal = obj.calMaps{ch};
        end

        function val = getHighPx(obj)
            val = obj.highPx/obj.numPx;
        end

        function val = getLowPx(obj)
            val = obj.lowPx/obj.numPx;
        end

        function setHighPx(obj,val)
            obj.highPx = round(obj.numPx*val);
            obj.calcRefPower();
            obj.calcAllTrans();
            obj.calcAllV();
        end

        function setLowPx(obj,val)
            obj.lowPx = round(obj.numPx*val);
            if (obj.lowPx == 0)
                obj.lowPx = 1;
            end
            obj.calcRefPower();
            obj.calcAllTrans();
            obj.calcAllV();
        end

        function setSliceSep(obj,val)
            obj.sliceSep = val;
            obj.z = linspace(0,(obj.sliceSep*(obj.slices-1)),obj.sliceSep);
            obj.calcAllTrans();
            obj.calcAllV();
        end

        function setMpOrder(obj,val)
            obj.mpOrder = val;
            obj.calcAllTrans();
            obj.calcAllV();
        end

        function setOrientation(obj,ch,val)
            obj.orientation(ch) = val;
            obj.calcAllTrans();
            obj.calcAllV();
        end

        function setType(obj,ch,val)
            obj.curveTypes(ch) = val;
            obj.calcTrans(ch);
            obj.calcAllV();
        end

        function setSL(obj,ch,val)
            obj.attenLen(ch) = val;
            obj.calcTrans(ch);
            obj.calcAllV();
        end

        function setCal(obj,ch,cal,calName)
            ch = ch*obj.sepChannels+1-obj.sepChannels;
            if (isempty(cal))
                obj.calNames{ch} = [];
                obj.calMaps{ch} = [];
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
            end

            obj.calcV(ch);
            if (~obj.sepChannels)
                obj.v = repmat(obj.v(ch,:),obj.channels,1);
            end
        end

        function success = setPrePower(obj,ch,pow,name)
            ch = ch*obj.sepChannels+1-obj.sepChannels;
            success = false;
            
            if (isempty(pow))
                obj.prePowerNames{ch} = [];
                obj.prePower(ch,:) = ones(1,obj.slices);
            else
                pow = pow(:,1);
                if ((length(pow)~=obj.slices) || (sum(isnan(pow))>0))
                    return;
                else
                    pow = abs(pow);
                    obj.prePower(ch,:) = pow'/max(pow);
                    obj.prePowerNames{ch} = name;
                end
            end

            obj.calcTrans(ch);
            obj.calcV(ch);
            if (~obj.sepChannels)
                obj.v = repmat(obj.v(ch,:),obj.channels,1);
            end

            success = true;
        end
    end
end