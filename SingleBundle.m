classdef SingleBundle
   properties
      R {mustBeNumeric} % axon radius
      L {mustBeNumeric} % size of the cell
   end
   methods
      function r = roundOff(obj)
         r = round([obj.Value],2);
      end
      function r = multiplyBy(obj,n)
         r = [obj.Value] * n;
      end
   end
end 