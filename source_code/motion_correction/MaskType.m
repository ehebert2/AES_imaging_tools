%%%%%%%%%%%%%%
% container for mask type ids used when adjusting ROIs
%%%%%%%%%%%%%%

classdef MaskType
    properties (Constant)
        EXPMASK = 1;
        AESALL = 2;
        AES = 3;
        SMPLALL = 4;
        SMPL = 5;
        NULL = 0;
        NOMASK = 6;
    end
end