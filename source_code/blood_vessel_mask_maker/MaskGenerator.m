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
    end

    properties (GetAccess = public, SetAccess = private)
        image
        sizes
        baseMask
        fullMask
        excludedMask
        border
    end

    properties (GetAccess = public, SetAccess = public)
        minFeature
        threshold
    end

    methods
        function obj = MaskGenerator(image)
            obj.threshold = max(image,[],'all');
            obj.minFeature = 0;
            obj.setBorder(0);
            obj.numNeighbors = size(obj.neighbors,1);

            obj.image = image;
            obj.baseMask = 0*image;
            obj.fullMask = 0*image;
            obj.excludedMask = 0*image;
            obj.numPx = numel(image);
            obj.x=repmat(1:size(image,2),size(image,1),1);
            obj.y=repmat((1:size(image,1))',1,size(image,2));
            obj.frontier = zeros(obj.numPx,2);
            obj.blob = image*0;

            obj.validNeighborMap = boolean(zeros(size(image,1),size(image,2),obj.numNeighbors));
            for jj=1:obj.numNeighbors
                neighborCoordY = obj.y+obj.neighbors(jj,1);
                neighborCoordX = obj.x+obj.neighbors(jj,2);
                obj.validNeighborMap(:,:,jj) = ((neighborCoordY > 0) .* (neighborCoordX > 0) .* (neighborCoordY <= size(image,1)) .* (neighborCoordX <= size(image,2))) > 0;
            end
        end

        function buildBaseMask(obj)
            obj.thresholdMask();
            obj.blobDet();
            obj.quickFilter();
        end

        function thresholdMask(obj)
            if (obj.threshold < 0)
                sensitivity = -(obj.threshold+1);
                obj.baseMask = imbinarize(obj.image,'adaptive','ForegroundPolarity','bright','Sensitivity',sensitivity);
            else
                obj.baseMask = obj.image > obj.threshold;
            end
            obj.excludedMask = 0*obj.excludedMask;
        end

        function blobDet(obj)
            obj.blobs = boolean([]);
            obj.searchMask = obj.baseMask+obj.excludedMask;
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
            obj.blobs = boolean(zeros(size(obj.blob,1),size(obj.blob,2),blockSize));
            obj.sizes = zeros(blockSize,1);

            while(seedInd<obj.numPx)
                frontierInd = 1;
                obj.frontier(1,:) = [obj.y(seedInd),obj.x(seedInd)];
                obj.searchMask(obj.frontier(1,1),obj.frontier(1,2)) = false;
                obj.blob = boolean(0*obj.blob);
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
                        tempBlobs = boolean(zeros(size(obj.blob,1),size(obj.blob,2),newLen));
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

        function quickFilter(obj)
            if (isempty(obj.sizes))
                obj.excludedMask = zeros(size(obj.image));
                obj.baseMask = zeros(size(obj.image));
                return;
            end

            ind=1;
            while ((ind <= length(obj.sizes)) && (obj.sizes(ind) < obj.minFeature))
                ind=ind+1;
            end

            if (ind==1)
                obj.excludedMask = boolean(zeros(size(obj.image)));
            else
                obj.excludedMask = (sum(obj.blobs(:,:,1:(ind-1)),3)>0);
            end

            if (ind > length(obj.sizes))
                obj.baseMask = boolean(zeros(size(obj.image)));
            else
                obj.baseMask = (sum(obj.blobs(:,:,ind:end),3)>0);
            end
        end

        function buildBorder(obj)
            obj.fullMask = conv2(obj.baseMask,obj.kernel,'same');
            obj.fullMask = obj.fullMask > 0;
        end

        function setBorder(obj,border)
            obj.border = border;
            obj.kernel = ones(2*obj.border+1,2*obj.border+1);
        end
    end
end