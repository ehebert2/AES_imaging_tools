%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matlab class for dealing with binary files used to store traces (matrix)
% can create an object to read through frame by frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef AESFile < handle
    properties (Constant)
        DOUBLE = -1;
    end

    properties (GetAccess=private)
        filename
        fin
    end

    properties (SetAccess=private)
        frames
        points
        current
        dataType
        bits
    end

    methods
        function obj = AESFile(filename)
            obj.filename = filename;
            obj.fin = fopen(filename,'r');
            a = fread(obj.fin,[1,2],'uint','b');
            obj.frames = a(1);
            obj.points = a(2);
            obj.current = 0;
            a = fread(obj.fin,[1,2],'int','b');
            obj.bits = a(1);
            switch (a(1))
                case 16
                    obj.dataType = 'int16';
                case 32
                    obj.dataType = 'int';
                case 1
                    obj.dataType = 'bit1';
                    a(2) = false;
                case AESFile.DOUBLE
                    obj.dataType = 'double';
                    a(2) = true;
                    obj.bits = 64;
            end

            if (a(2)==0)
                obj.dataType = strcat('u',obj.dataType);
            end
        end

        function frame = getNextFrame(obj)
            if (obj.current < obj.frames)
                frame = fread(obj.fin,[1,obj.points],obj.dataType,'b');
                obj.current = obj.current + 1;
            else
                frame = [];
            end
        end

        function frame = getFrame(obj,ind)
            if ((ind < 1) || (ind > obj.frames))
                frame = [];
            elseif (ind == (obj.current + 1))
                frame = fread(obj.fin,[1,obj.points],obj.dataType,'b');
                obj.current = obj.current + 1;
            else
                if (obj.bits == 1)
                    fseek(obj.fin,(16+floor((ind-1)/8)*obj.points),-1);
                    if (mod(ind-1,8)~=0)
                        fread(obj.fin,[1,mod(ind-1,8)],obj.dataType,'b');
                    end
                else
                    fseek(obj.fin,(16+(ind-1)*(obj.bits/8)*obj.points),-1);
                end
                frame = fread(obj.fin,[1,obj.points],obj.dataType,'b');
                obj.current = ind;
            end
        end

        function close(obj)
            fclose(obj.fin);
        end
    end
    
    methods (Static)
        % reads full matrix into memory without needing to instantiate
        % object
        function data = readFullFile(filename)
            fin = fopen(filename,'r');
            a = fread(fin,[1,2],'uint','b');
            frames = a(1);
            points = a(2);

            a = fread(fin,[1,2],'int','b');
            switch (a(1))
                case 16
                    type = 'int16';
                case 32
                    type = 'int';
                case 1
                    type = 'bit1';
                    a(2) = false;
                case AESFile.DOUBLE
                    type = 'double';
                    a(2) = true;
            end

            if (a(2)==0)
                type = strcat('u',type);
            end

            data = zeros(frames,points);
            for ii = 1:frames
                data(ii,:) = fread(fin,[1,points],type,'b');
            end

            fclose(fin);
        end

        % sets up file writer and writes initial parameters to file
        % you need to close the writer manually and also keep data bit
        % depth and signed the same throughout
        function writer = getWriter(filepath,cols,frames,bits,signed)
            writer = fopen(filepath,'w');
            fwrite(writer,[frames,cols],'uint','b');
            switch bits
                case 1
                    signed = false;
                case AESFile.DOUBLE
                    signed = true;
            end

            fwrite(writer,[bits,signed],'int','b');
        end

        % writes full matrix to file
        function success = writeToFile(filepath,A,bits,signed)
            success = false;
            fout = fopen(filepath, 'w');
            if (length(size(A))~=2)
                return;
            end

            fwrite(fout,size(A),'uint','b');
            
            switch bits
                case 16
                    type = 'int16';
                case 32
                    type = 'int';
                case 1
                    type = 'bit1';
                    signed = false;
                case AESFile.DOUBLE
                    type = 'double';
                    signed = false;
            end
            if (~signed)
                type = strcat('u',type);
            end

            fwrite(fout,[bits,signed],'int','b');
            for ii=1:size(A,1)
                fwrite(fout,A(ii,:),type,'b');
            end
            fclose(fout);
            success = true;
        end
    end
end