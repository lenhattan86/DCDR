tic
data = rand(5e6,1);  % pre-processing (5M elements, ~40MB)
javaaddpath 'C:\Yair\Code\'  % path to MyJavaThread.class
%start(MyJavaThread('F:\test.data',data));  % start running in parallel
data = fft(data);  % post-processing (Java I/O runs in parallel)
toc