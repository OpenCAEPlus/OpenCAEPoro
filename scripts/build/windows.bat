SET "system=windows"
SET "compiler=%1"
IF "%compiler%"=="" SET "compiler=intel"
SET "build=%2"
IF "%build%"=="" SET "build=Debug"
SET "target=%3"
IF "%target%"=="" SET "target=all"
ECHO %system%-%compiler%-%build%-%target%


SET CMAKE="C:\\PROGRAM FILES (X86)\\MICROSOFT VISUAL STUDIO\\2019\\COMMUNITY\\COMMON7\\IDE\\COMMONEXTENSIONS\\MICROSOFT\\CMAKE\\CMake\\bin\\cmake.exe"

IF "%compiler%"=="intel" call "C:\\Program Files (x86)\\Intel\\oneAPI\\setvars.bat"

%CMAKE% --preset="%system%-%compiler%-%build%" -S "."
%CMAKE% --build --preset="%system%-%compiler%-%build%" --target %target%