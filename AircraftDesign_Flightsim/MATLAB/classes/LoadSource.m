classdef (Abstract) LoadSource < handle
    methods (Abstract)
        fm = get_force_moment(obj, x, u)
    end
end