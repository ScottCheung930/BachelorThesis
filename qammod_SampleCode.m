%  qammod Quadrature amplitude modulation
%  
%     Y = qammod(X,M) outputs the complex envelope of the modulation of the
%     message signal X using quadrature amplitude modulation. M is the
%     alphabet size and must be an integer power of two. The message signal X
%     must consist of integers between 0 and M-1. X can be a scalar, a
%     vector, a matrix or an array with 3 dimensions.
%  
%     Y = qammod(X,M,SYMBOL_ORDER) specifies how the function maps an integer
%     or group of log2(M) input bits to the corresponding symbol. If
%     SYMBOL_ORDER is set to 'gray', then the function uses a Gray-coded
%     ordering. If SYMBOL_ORDER is set to 'bin', then the function uses a
%     natural binary-coded ordering. If SYMBOL_ORDER is an integer-valued
%     vector with M elements, the function uses the ordering specified by
%     this vector. This vector must have unique elements in the range [0,
%     M-1]. The first element of this vector corresponds to the top-leftmost
%     point of the constellation, with subsequent elements running down
%     column-wise, from left to right. The last element corresponds to the
%     bottom-rightmost point. The default value is 'gray'.
%  
%     Y = qammod(X,M,...,Name,Value) specifies additional name-value pair
%     arguments described below:
%  
%     'InputType'          One of the strings: 'integer', or 'bit'. 'integer'
%                          indicates that the message signal is integer
%                          valued between 0 and M-1. 'bit' indicates that the
%                          message signal is binary (0 or 1). In this case,
%                          the number of rows (dimension 1) must be an
%                          integer multiple of log2(M). A group of log2(M)
%                          bits are mapped onto a symbol, with the first bit
%                          representing the MSB and the last bit representing
%                          the LSB. The default value is 'integer'.
%  
%     'UnitAveragePower'   A logical scalar value. If true, the QAM
%                          constellation is scaled to average power of 1. If
%                          false, the QAM constellation with minimum distance
%                          of 2 between constellation points is used. The
%                          default value is false.
%  
%     'OutputDataType'     Output fixed-point type as a signed, unscaled
%                          numerictype object in MATLAB simulation, and as a
%                          signed, scaled numerictype object during C code or
%                          MEX generation. When this argument is not
%                          specified, if the input datatype is double or
%                          built-in integers, the output datatype is double;
%                          if the input datatype is single, the output
%                          datatype is single. When the input is fixed-point,
%                          this parameter must be specified.
%  
%     'PlotConstellation'  A logical scalar value. If true, the QAM
%                          constellation is plotted. The default value is
%                          false. The input X is processed and the modulated
%                          signal is returned in output Y.
 
%%   Example 1:
     % 32-QAM modulation. Default: Integer input, Gray coding, minimum
     % distance of 2 between constellation points
     x = (0:31)';
     y = qammod(x, 32);
  
%%   Example 2:
     % 16-QAM modulation, with LTE specific symbol mapping and constellation
     % scaled to average power of 1. Default: Integer input
     x = randi([0, 15], 20, 4, 2);
     lteSymMap = [11 10 14 15 9 8 12 13 1 0 4 5 3 2 6 7];
     y = qammod(x, 16, lteSymMap, 'UnitAveragePower', true);
  
%%   Example 3: 
     % 64-QAM modulation with binary mapping, bit input and signed
     % fixed-point output data type with 16 bits of word length and 10 
     % bits of fraction length. Default: minimum distance of 2 between
     % constellation points
     x = randi([0, 1], 10*log2(64), 3);
     y = qammod(x, 64, 'bin', 'InputType', 'bit', 'OutputDataType', numerictype(1,16,10));
 
%%   Example 4:
     % Visualize the constellation for 16-QAM modulation, with gray
     % mapping, bit input and constellation scaled to average power of 1.
     x = randi([0, 1], log2(16), 1);
     y = qammod(x, 16, 'InputType', 'bit', 'UnitAveragePower', true, 'PlotConstellation', true);
  