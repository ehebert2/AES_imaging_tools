classdef MaskGenerator < handle
    properties (GetAccess = private, SetAccess = private)
        ignoreSize = 3;
        neighbors = [-1, 0; 1, 0; 0, -1; 0, 1];
        validNeighborMap
        numNeighbors
        frontier
        searchMask
        blobs
        blob
        x
        y
        numPx
        kernel
        roiInd
        filterInd
    end

    properties (GetAccess = public, SetAccess = private)
        image
        sizes
        border
        filterMask
        excludedFilterMask
        fullMask
        roiSelected
        selectedMask
        maxImVal
        imBinX
        imBinY
        otsuCurve
        otsuThresh
        
        % threshold types
        GLOBAL = 0;
        OTSU = 1;
        ADAPTIVE = 2;
    end

    properties (GetAccess = public, SetAccess = public)
        thresholdType
        minFeature
        threshold
        baseMask
        excludedMask
        maxHole
        sensitivity
    end

    methods
        function obj = MaskGenerator(image)
            obj.maxImVal = max(image,[],'all');
            obj.thresholdType = obj.GLOBAL;
            obj.threshold = obj.maxImVal;
            obj.sensitivity = 0.5;
            obj.minFeature = 0;
            obj.setBorder(0);
            obj.numNeighbors = size(obj.neighbors,1);

            obj.image = image;
            obj.baseMask = 0*image;
            obj.fullMask = 0*image;
            obj.excludedMask = 0*image;
            obj.filterMask = 0*image;
            obj.excludedFilterMask = 0*image;
            obj.selectedMask = 0*image;
            obj.numPx = numel(image);
            obj.x=repmat(1:size(image,2),size(image,1),1);
            obj.y=repmat((1:size(image,1))',1,size(image,2));
            obj.frontier = zeros(obj.numPx,2);
            obj.blob = image*0;
            obj.roiSelected = false;
            obj.maxHole = 0;
            obj.filterInd = 1;

            obj.validNeighborMap = zeros(size(image,1),size(image,2),obj.numNeighbors)>0;
            for jj=1:obj.numNeighbors
                neighborCoordY = obj.y+obj.neighbors(jj,1);
                neighborCoordX = obj.x+obj.neighbors(jj,2);
                obj.validNeighborMap(:,:,jj) = ((neighborCoordY > 0) .* (neighborCoordX > 0) .* (neighborCoordY <= size(image,1)) .* (neighborCoordX <= size(image,2))) > 0;
            end

            obj.calculateOtsu();
        end

        function buildBaseMask(obj)
            obj.thresholdMask();
            obj.blobDet();
            obj.quickFilter();
        end

        function thresholdMask(obj)
            if (obj.thresholdType == obj.ADAPTIVE)
                minVal = min(obj.image,[],'all');
                obj.baseMask = imbinarize(((obj.image-minVal)/(max(obj.image,[],'all')-minVal)),'adaptive','ForegroundPolarity','bright','Sensitivity',obj.sensitivity);
            elseif (obj.thresholdType == obj.OTSU)
                obj.baseMask = obj.image > obj.otsuThresh;
            else
                obj.baseMask = obj.image > obj.threshold;
            end
            obj.excludedMask = 0*obj.excludedMask;

            if ((obj.maxHole > 1) && ((obj.threshold < obj.maxImVal) || (obj.thresholdType ~= obj.GLOBAL)))
                obj.removeHoles();
            end
        end

        function blobDet(obj)
            obj.blobs = [];
            obj.searchMask = (obj.baseMask+obj.excludedMask)>0;
            obj.frontier = 0*obj.frontier;

            seedInd = 1;
            while((seedInd < obj.numPx) && ~obj.searchMask(seedInd))
                seedInd = seedInd+1;
                if (seedInd == obj.numPx)
                    obj.sizes = [];
                    return;
                end
            end

            blockSize = 250;
            blobInd = 0;
            obj.blobs = zeros(size(obj.blob,1),size(obj.blob,2),blockSize)>0;
            obj.sizes = zeros(blockSize,1);

            while(seedInd<obj.numPx)
                frontierInd = 1;
                obj.frontier(1,:) = [obj.y(seedInd),obj.x(seedInd)];
                obj.searchMask(obj.frontier(1,1),obj.frontier(1,2)) = false;
                obj.blob = (0*obj.blob)>0;
                obj.blob(obj.frontier(1,1),obj.frontier(1,2)) = true;
                while (frontierInd > 0)
                    coord = obj.frontier(frontierInd,:);
                    frontierInd = frontierInd-1;
                    for ii=1:obj.numNeighbors
                        if (obj.validNeighborMap(coord(1),coord(2),ii))
                            tempCoord = [(coord(1)+obj.neighbors(ii,1)),(coord(2)+obj.neighbors(ii,2))];
                            if (obj.searchMask(tempCoord(1),tempCoord(2)))
                                frontierInd = frontierInd+1;
                                obj.frontier(frontierInd,:) = tempCoord;
                                obj.searchMask(tempCoord(1),tempCoord(2)) = false;
                                obj.blob(tempCoord(1),tempCoord(2)) = true;
                            end
                        end
                    end
                end

                while (~obj.searchMask(seedInd) && (seedInd<obj.numPx))
                    seedInd = seedInd+1;
                end
                blobSize = sum(obj.blob,'all');
                if (blobSize>obj.ignoreSize)
                    if ((blobInd)==length(obj.sizes))
                        newLen = length(obj.sizes)+blockSize;
                        tempSizes = zeros(newLen,1);
                        tempBlobs = zeros(size(obj.blob,1),size(obj.blob,2),newLen)>0;
                        tempSizes(1:blobInd) = obj.sizes;
                        tempBlobs(:,:,1:blobInd) = obj.blobs;
                        obj.sizes = tempSizes;
                        obj.blobs = tempBlobs;
                    end
                    blobInd = blobInd + 1;
                    obj.sizes(blobInd) = blobSize;
                    obj.blobs(:,:,blobInd) = obj.blob;
                end
            end
            obj.sizes = sqrt(obj.sizes(1:blobInd));
            obj.blobs = obj.blobs(:,:,1:blobInd);

            [obj.sizes,ind] = sort(obj.sizes);
            obj.blobs = obj.blobs(:,:,ind);
            obj.sizes = sqrt(obj.sizes);
        end

        function removeHoles(obj)
            obj.searchMask = ~(obj.baseMask>0);
            obj.frontier = 0*obj.frontier;

            seedInd = 1;
            while((seedInd < obj.numPx) && ~obj.searchMask(seedInd))
                seedInd = seedInd+1;
                if (seedInd == obj.numPx)
                    return;
                end
            end

            while(seedInd<obj.numPx)
                frontierInd = 1;
                obj.frontier(1,:) = [obj.y(seedInd),obj.x(seedInd)];
                obj.searchMask(obj.frontier(1,1),obj.frontier(1,2)) = false;
                obj.blob = (0*obj.blob)>0;
                obj.blob(obj.frontier(1,1),obj.frontier(1,2)) = true;
                while (frontierInd > 0)
                    coord = obj.frontier(frontierInd,:);
                    frontierInd = frontierInd-1;
                    for ii=1:obj.numNeighbors
                        if (obj.validNeighborMap(coord(1),coord(2),ii))
                            tempCoord = [(coord(1)+obj.neighbors(ii,1)),(coord(2)+obj.neighbors(ii,2))];
                            if (obj.searchMask(tempCoord(1),tempCoord(2)))
                                frontierInd = frontierInd+1;
                                obj.frontier(frontierInd,:) = tempCoord;
                                obj.searchMask(tempCoord(1),tempCoord(2)) = false;
                                obj.blob(tempCoord(1),tempCoord(2)) = true;
                            end
                        end
                    end
                end

                while (~obj.searchMask(seedInd) && (seedInd<obj.numPx))
                    seedInd = seedInd+1;
                end

                if (sum(obj.blob,'all')<obj.maxHole)
                    obj.baseMask = obj.baseMask + obj.blob;
                end
            end
        end

        function quickFilter(obj)
            if (isempty(obj.sizes))
                obj.excludedMask = zeros(size(obj.image));
                obj.baseMask = zeros(size(obj.image));
                return;
            end

            obj.filterInd=1;
            while ((obj.filterInd <= length(obj.sizes)) && (obj.sizes(obj.filterInd) < obj.minFeature))
                obj.filterInd=obj.filterInd+1;
            end

            if (obj.filterInd==1)
                obj.excludedFilterMask = zeros(size(obj.image))>0;
            else
                obj.excludedFilterMask = (sum(obj.blobs(:,:,1:(obj.filterInd-1)),3)>0);
            end

            if (obj.filterInd > length(obj.sizes))
                obj.filterMask = zeros(size(obj.image))>0;
            else
                obj.filterMask = (sum(obj.blobs(:,:,obj.filterInd:end),3)>0);
            end

            obj.excludedMask = obj.excludedFilterMask;
            obj.baseMask = obj.filterMask;
        end

        function buildBorder(obj)
            obj.fullMask = conv2(obj.baseMask,obj.kernel,'same');
            obj.fullMask = obj.fullMask > 0;
            obj.excludedMask = 0*obj.excludedMask;
        end

        function setBorder(obj,border)
            obj.border = border;
            if (border==1)
                obj.kernel = ones(3,3);
            else
                tempX = repmat(-border:border,(2*border+1),1);
                tempY = tempX';
                obj.kernel = ((tempX.^2+tempY.^2) <= border^2);
            end
        end

        function setImage(obj,image)
            obj.image = image;
            obj.maxImVal = max(image,[],'all');
            obj.calculateOtsu();
            obj.buildBaseMask();
        end

        function clckDetect(obj, x, y, hold)  
            if (isempty(obj.sizes) || (obj.filterInd==length(obj.sizes)))
                obj.roiSelected = false;
                obj.selectedMask = obj.selectedMask*0;
                return;
            end

            if (hold && obj.roiSelected)
                temp = obj.blobs(y,x,:);
                temp(1:(obj.filterInd-1))=0;
                obj.roiInd = (obj.roiInd+temp)>0;
            else
                obj.roiInd = obj.blobs(y,x,:);
                obj.roiInd(1:(obj.filterInd-1))=0;
            end

            obj.roiSelected = sum(obj.roiInd)>0;
            if (obj.roiSelected)
                obj.selectedMask = sum(obj.blobs(:,:,obj.roiInd),3)>0;
            else
                obj.selectedMask = obj.selectedMask*0;
            end
        end

        function deleteROI(obj)
            if (obj.roiSelected)
                obj.blobs = obj.blobs(:,:,~obj.roiInd);
                obj.sizes = obj.sizes(~obj.roiInd);
                obj.quickFilter();
                obj.buildBorder();
                obj.deselectROI();
            end
        end

        function deselectROI(obj)
            obj.roiSelected = false;
            obj.selectedMask = obj.selectedMask*0;
        end

        function printBlobs(obj,basename,path,ind)
            if (isempty(obj.sizes))
                return;
            end

            tempInd = 1;
            while((tempInd <= length(obj.sizes)) && (obj.sizes(tempInd)<obj.minFeature))
                tempInd = tempInd + 1;
            end

            if (tempInd > length(obj.sizes))
                return;
            end

            for ii=tempInd:length(obj.sizes)
                ff = fullfile(path,strcat(basename,'_',num2str(ind+ii-tempInd+1),'.tif'));
                imwrite((obj.blobs(:,:,ii)+0),ff, 'Compression','none');
            end
        end

        function calculateOtsu(obj)
            [obj.imBinY,obj.imBinX] = hist(obj.image(:),50);
            obj.otsuCurve = obj.imBinX(1:(end-1))*0;
            muPre = obj.imBinX.*obj.imBinY;
            for ii=1:(length(obj.imBinX)-1)
                w0 = sum(obj.imBinY(1:ii));
                w1 = sum(obj.imBinY((ii+1):end));
                mu0 = sum(muPre(1:ii))/w0;
                mu1 = sum(muPre((ii+1):end))/w1;
                obj.otsuCurve(ii) = w0*w1*(mu0-mu1)^2;
            end
            [~,ind] = max(obj.otsuCurve);
            obj.otsuThresh = obj.imBinX(ind);
        end
    end
end