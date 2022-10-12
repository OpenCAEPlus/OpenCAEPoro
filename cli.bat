@ECHO OFF
SET system=windows
SET compiler=intel
SET build=Debug
SET target=all

SET CMAKE="C:\\PROGRAM FILES (X86)\\MICROSOFT VISUAL STUDIO\\2019\\COMMUNITY\\COMMON7\\IDE\\COMMONEXTENSIONS\\MICROSOFT\\CMAKE\\CMake\\bin\\cmake.exe"
SET CTEST="C:\\PROGRAM FILES (X86)\\MICROSOFT VISUAL STUDIO\\2019\\COMMUNITY\\COMMON7\\IDE\\COMMONEXTENSIONS\\MICROSOFT\\CMAKE\\CMake\\bin\\ctest.exe"



SET argC=0
FOR %%x IN (%*) DO SET /A argC+=1

IF %argC%==0 (
    ECHO You need to use a subcommand of either 'build' or 'test'
    GOTO end
)

SET option=%1
SHIFT

IF %option%==build GOTO optionBuild
IF %option%==test GOTO optionTest

ECHO Subcommand not recognized, you must choose either 'build' or 'test'

@REM Start the build subcommand part
:optionBuild

IF %argC%==1 GOTO optionBuildHelp

@REM Get command line options
:optionBuildLoop

IF NOT "%1"=="" (
    IF "%1"=="--build" GOTO :optionBuildBuild
    IF "%1"=="-b" GOTO optionBuildBuild
:optionBuildLoopAfterBuild
    IF "%1"=="--compiler" GOTO optionBuildCompiler 
    IF "%1"=="-c" GOTO optionBuildCompiler 
:optionBuildLoopAfterCompiler
    IF "%1"=="--target" GOTO optionBuildTarget
    IF "%1"=="-t" GOTO optionBuildTarget
:optionBuildLoopAfterTarget
    IF "%1"=="--help" GOTO optionBuildHelp 
    IF "%1"=="-h" GOTO optionBuildHelp 
    SHIFT
    GOTO optionBuildLoop
)

ECHO %system%-%compiler%-%build%-%target%

IF %compiler%==intel CALL "C:\\Program Files (x86)\\Intel\\oneAPI\\setvars.bat"

%CMAKE% --preset="%system%-%compiler%-%build%" -S "."
%CMAKE% --build --preset="%system%-%compiler%-%build%" --target %target%

GOTO end

:optionBuildBuild
SET build=%2
SHIFT
GOTO optionBuildLoopAfterBuild

:optionBuildCompiler
SET compiler=%2
SHIFT
GOTO optionBuildLoopAfterCompiler

:optionBuildTarget
SET target=%2
SHIFT
GOTO optionBuildLoopAfterTarget

:optionBuildHelp
ECHO CLI subcommand "build" to compile the source codes.
ECHO.
ECHO Usage:
ECHO The command will automatically detect your system, you just need to tell the 
ECHO cli which compiler, build type, and target you want to use
ECHO.
ECHO     "cli.bat build -c|--compiler -b|--build -t|--target [-h|--help]"
ECHO.
ECHO Example:
ECHO.
ECHO     "cli.bat build -c intel -b Debug"
ECHO     "cli.bat build -c intel -b Debug -t all"
ECHO     "cli.bat build -c intel -b Debug -t clean"

GOTO end

@REM Start the test subcommand part
:optionTest

IF %argC%==1 GOTO optionTestHelp

@REM Get command line options
:optionTestLoop
IF NOT "%1"=="" (
    IF "%1"=="--build" GOTO optionTestBuild
    IF "%1"=="-b" GOTO optionTestBuild
:optionTestLoopAfterBuild
    IF "%1"=="--compiler" GOTO optionTestCompiler 
    IF "%1"=="-c" GOTO optionTestCompiler 
:optionTestLoopAfterCompiler
    IF "%1"=="--help" GOTO optionTestHelp 
    IF "%1"=="-h" GOTO optionTestHelp 
    SHIFT
    GOTO :optionTestLoop
)

ECHO %system%-%compiler%-%build%-%target%

IF %compiler%==intel CALL "C:\\Program Files (x86)\\Intel\\oneAPI\\setvars.bat"

%CMAKE% --preset="%system%-%compiler%-%build%" -S "."
%CMAKE% --build --preset="%system%-%compiler%-%build%" --target %target%
%CTEST% --preset="%system%-%compiler%-%build%" 

GOTO end

:optionTestBuild
SET build=%2
SHIFT
GOTO optionTestLoopAfterBuild

:optionTestCompiler
SET compiler=%2
SHIFT
GOTO optionTestLoopAfterCompiler

:optionTestdHelp
ECHO CLI subcommand "test" to compile and test the source codes.
ECHO.
ECHO Usage:
ECHO The command will automatically detect your system, you just need to tell the 
ECHO cli which compiler, build type, target, and the testing type, either unit or benchmark
ECHO.
ECHO     "cli.bat test -c|--compiler -b|--build -t|--target -u|--unit -m|--benchmark [-h|--help]"
ECHO.
ECHO Example:
ECHO.
ECHO     "cli.bat test -c intel -b Debug"
ECHO     "cli.bat test -c intel -b Debug -u"
ECHO     "cli.bat test -c intel -b Debug -m -u"
GOTO end

:end
