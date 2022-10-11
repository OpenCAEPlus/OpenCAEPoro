SET "system=windows"
SET "compiler=%1"
IF "%compiler%"=="" SET "compiler=intel"
SET "build=%2"
IF "%build%"=="" SET "build=Release"
ECHO "------Testing %system%-%compiler%-%build%"

SET CMAKE="C:\\PROGRAM FILES (X86)\\MICROSOFT VISUAL STUDIO\\2019\\COMMUNITY\\COMMON7\\IDE\\COMMONEXTENSIONS\\MICROSOFT\\CMAKE\\CMake\\bin\\cmake.exe"
SET CTEST="C:\\PROGRAM FILES (X86)\\MICROSOFT VISUAL STUDIO\\2019\\COMMUNITY\\COMMON7\\IDE\\COMMONEXTENSIONS\\MICROSOFT\\CMAKE\\CMake\\bin\\ctest.exe"


call scripts/build/windows.bat %compiler% %build%
%CTEST% --preset="%system%-%compiler%-%build%" 
