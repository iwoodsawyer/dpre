%MAKE_FACTOR Compilation of DPRE mex-files
MATLAB_PATH = matlabroot;
COMPILE_OPTIONS = '';
v = ver('matlab');
matver = sscanf(v.Version, '%d.%d.%d')';
COMPILE_OPTIONS = [ COMPILE_OPTIONS ' -DMATLAB_VERSION=0x' sprintf('%02d%02d', matver(1), matver(2)) ];
MATLAB_VERSION = matver(1) + matver(2)/100;

if MATLAB_VERSION < 7.14
    error('Requires Matlab 2012a (7.14) or higher with Lapack version 3.3.0 or higher.')
end

if strcmpi('GLNX86', computer) || strcmpi('GLNXA64', computer) ...
        || strcmpi('MACI', computer) || strcmpi('MAC', computer) ...
        || strcmpi('MACI64', computer)
    % GNU/Linux (x86-32 or x86-64) or MacOS (Intel or PPC)
    if MATLAB_VERSION < 7.05
        BLAS_PATH = '';
    else
        BLAS_PATH = ' -lmwblas';
    end
    LAPACK_PATH = ' -lmwlapack';
    SLICOT_PATH = ' -lmwslicot';
    if strcmpi('GLNX86', computer) || strcmpi('GLNXA64', computer)
        COMPILE_OPTIONS = [COMPILE_OPTIONS,' -DSkip_f2c_Undefs',' -DNON_UNIX_STDIO'];
    else
        COMPILE_OPTIONS = [COMPILE_OPTIONS,' -DSkip_f2c_Undefs'];
    end
elseif strcmpi('PCWIN', computer) || strcmpi('PCWIN64', computer)
    % Windows (x86-32 or x86-64)
    if strcmpi('PCWIN', computer)
        error('Unsupported platform')
    else
        cc = mex.getCompilerConfigurations('C','Selected');
        MANUFACTURER = cc.Manufacturer;
        switch lower(MANUFACTURER)
            case {'gnu'}
                BLAS_PATH = fullfile(MATLAB_PATH, 'extern', 'lib', 'win64', 'mingw64', 'libmwblas.lib');
                LAPACK_PATH = fullfile(MATLAB_PATH, 'extern', 'lib', 'win64', 'mingw64', 'libmwlapack.lib');
                SLICOT_PATH = fullfile(pwd, 'mingw64', 'libmwslicot.lib');
            case {'microsoft'}
                BLAS_PATH = fullfile(MATLAB_PATH, 'extern', 'lib', 'win64', 'microsoft', 'libmwblas.lib');
                LAPACK_PATH = fullfile(MATLAB_PATH, 'extern', 'lib', 'win64', 'microsoft', 'libmwlapack.lib');
                SLICOT_PATH = fullfile(pwd, 'microsoft', 'libmwslicot.lib');
            otherwise
                disp('Try "mex -setup", because BLAS/LAPACK library is not available!')
        end
    end
    if MATLAB_VERSION < 7.05
        BLAS_PATH = ''; % On <= 7.4, BLAS in included in LAPACK
    else
        BLAS_PATH = [' "' BLAS_PATH '"'];
    end
    LAPACK_PATH = [' "' LAPACK_PATH '"'];
    SLICOT_PATH = [' "' SLICOT_PATH '"'];
    COMPILE_OPTIONS = [COMPILE_OPTIONS,' -DMSDOS',' -DUSE_CLOCK',' -DNO_ONEXIT'];
else
    error('Unsupported platform')
end

% Large array dims for 64 bits platforms appeared in Matlab 7.3
if (strcmpi('GLNXA64', computer) || strcmpi('PCWIN64', computer) ...
        || strcmpi('MACI64', computer))
    if ~(MATLAB_VERSION < 9.04)
        COMPILE_OPTIONS = [ COMPILE_OPTIONS ' -R2018a' ];
    elseif ~(MATLAB_VERSION < 7.03)
        COMPILE_OPTIONS = [ COMPILE_OPTIONS ' -largeArrayDims' ];
    end
end

% Comment next line to suppress optimization
COMPILE_OPTIONS = [ ' -O' COMPILE_OPTIONS ];

% Comment next line to suppress compilation debugging info
%COMPILE_OPTIONS = [ ' -v' COMPILE_OPTIONS ];

disp('Compiling pdare...')
eval(['mex ', COMPILE_OPTIONS, ' dpre.c', BLAS_PATH, LAPACK_PATH, SLICOT_PATH]);
