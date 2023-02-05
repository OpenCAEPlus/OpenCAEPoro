REM Regression test for OpenCAEPoro project

cls
@echo off
set mydate=%date%
set mytime=%time%
echo Regression test at: %mydate% %mytime%


set dd=%date:~8,2%
set mm=%date:~5,2%
set yy=%date:~0,4%
set Tmm=%time:~3,2%
set Thh=%time:~0,2%
set myfolder="%yy%-%mm%-%dd%-%Thh%-%Tmm%\"

set mylog="%yy%-%mm%-%dd%-%Thh%-%Tmm%-log.dat"

set OCP=OpenCAEPoro.exe
set OCP_Log=log.out
set OCP_SUMMARY=SUMMARY.out
set OCP_FastReview=FastReview.out

set myFlag=Y

set SPE1a_flag=Y
set SPE1a_dir=..\examples\spe1a\
set SPE1a_test=%SPE1a_dir%spe1a.data  	method=FIM dtInit=1 dtMax=100 dtMin=0.1

set SPE1b_flag=Y
set SPE1b_dir=..\examples\spe1b\
set SPE1b_test=%SPE1b_dir%spe1b.data  	method=FIM dtInit=1 dtMax=100 dtMin=0.1

set SPE9_flag=Y
set SPE9_dir=..\examples\spe9\
set SPE9_test=%SPE9_dir%spe9_FIM.data

set CP_flag=Y
set CP_dir=..\examples\cornerpoint\
set CP_test=%CP_dir%CP.data           	method=FIM dtInit=1 dtMax=100 dtMin=0.1

set SPE3_flag=Y
set SPE3_dir=..\examples\spe3\
set SPE3_test=%SPE3_dir%spe3.data    	method=FIM dtInit=1 dtMax=100 dtMin=0.1

set SPE5_flag=Y
set SPE5_dir=..\examples\spe5\
set SPE5_test=%SPE5_dir%spe5.data    	method=FIM dtInit=1 dtMax=50  dtMin=0.1

set ZBCO2_flag=Y
set ZBCO2_dir=..\examples\3.ZB-CO2\OCP\
set ZBCO2_test=%ZBCO2_dir%ZB-CO2.data  	method=FIM dtInit=1 dtMax=10  dtMin=0.1 pl=2


copy ..\..\x64\Release\%OCP% %OCP%


REM Run SPE1a test
if /i %SPE1a_flag% == %myFlag% (	
	echo.
	echo Test problem %SPE1a_test%  >  %mylog%
	%OCP%   %SPE1a_test%   2>&1  | tee %SPE1a_dir%%OCP_Log%
	
	md .\SPE1a\%myfolder%
	
	copy %SPE1a_dir%%OCP_Log%  			.\SPE1a\%myfolder%%OCP_Log%
	copy %SPE1a_dir%%OCP_SUMMARY% 		.\SPE1a\%myfolder%%OCP_SUMMARY%
	copy %SPE1a_dir%%OCP_FastReview% 	.\SPE1a\%myfolder%%OCP_FastReview%
	
	fc   .\SPE1a\%myfolder%%OCP_Log%  		  .\SPE1a\current\%OCP_Log% 		> 	.\SPE1a\CompareLog.out
	fc   .\SPE1a\%myfolder%%OCP_SUMMARY%  	  .\SPE1a\current\%OCP_SUMMARY% 	> 	.\SPE1a\CompareSummary.out
	fc   .\SPE1a\%myfolder%%OCP_FastReview%   .\SPE1a\current\%OCP_FastReview% 	> 	.\SPE1a\CompareFastRevive.out
	
	fc   .\SPE1a\%myfolder%%OCP_Log%  		  .\SPE1a\current\%OCP_Log% 		>> 	%mylog%
	fc   .\SPE1a\%myfolder%%OCP_SUMMARY%  	  .\SPE1a\current\%OCP_SUMMARY% 	>> 	%mylog%
	fc   .\SPE1a\%myfolder%%OCP_FastReview%   .\SPE1a\current\%OCP_FastReview% 	>> 	%mylog%
)




REM Run SPE1b test
if /i %SPE1b_flag% == %myFlag% (
	echo.
	echo Test problem %SPE1b_test%  >>	%mylog%
	%OCP%   %SPE1b_test%   2>&1  | tee %SPE1b_dir%%OCP_Log%
	
	md .\SPE1b\%myfolder%
	
	copy %SPE1b_dir%%OCP_Log%  			.\SPE1b\%myfolder%%OCP_Log%
	copy %SPE1b_dir%%OCP_SUMMARY% 		.\SPE1b\%myfolder%%OCP_SUMMARY%
	copy %SPE1b_dir%%OCP_FastReview% 	.\SPE1b\%myfolder%%OCP_FastReview%
	
	fc   .\SPE1b\%myfolder%%OCP_Log%          .\SPE1b\current\%OCP_Log%          >   .\SPE1b\CompareLog.out
	fc   .\SPE1b\%myfolder%%OCP_SUMMARY%      .\SPE1b\current\%OCP_SUMMARY%      >   .\SPE1b\CompareSummary.out
	fc   .\SPE1b\%myfolder%%OCP_FastReview%   .\SPE1b\current\%OCP_FastReview%   >   .\SPE1b\CompareFastRevive.out
	
	fc   .\SPE1b\%myfolder%%OCP_Log%          .\SPE1b\current\%OCP_Log%          >>   %mylog%
	fc   .\SPE1b\%myfolder%%OCP_SUMMARY%      .\SPE1b\current\%OCP_SUMMARY%      >>   %mylog%
	fc   .\SPE1b\%myfolder%%OCP_FastReview%   .\SPE1b\current\%OCP_FastReview%   >>   %mylog%
)




REM Run SPE9 test
if /i %SPE9_flag% == %myFlag% (
	echo.
	echo Test problem %SPE9_test%  >>	%mylog%
	%OCP%   %SPE9_test%   2>&1  | tee %SPE9_dir%%OCP_Log%
	
	md .\SPE9\%myfolder%
	
	copy %SPE9_dir%%OCP_Log%  			.\SPE9\%myfolder%%OCP_Log%
	copy %SPE9_dir%%OCP_SUMMARY% 		.\SPE9\%myfolder%%OCP_SUMMARY%
	copy %SPE9_dir%%OCP_FastReview% 	.\SPE9\%myfolder%%OCP_FastReview%
	
	fc   .\SPE9\%myfolder%%OCP_Log%          .\SPE9\current\%OCP_Log%          >   .\SPE9\CompareLog.out
	fc   .\SPE9\%myfolder%%OCP_SUMMARY%      .\SPE9\current\%OCP_SUMMARY%      >   .\SPE9\CompareSummary.out
	fc   .\SPE9\%myfolder%%OCP_FastReview%   .\SPE9\current\%OCP_FastReview%   >   .\SPE9\CompareFastRevive.out
	
	fc   .\SPE9\%myfolder%%OCP_Log%          .\SPE9\current\%OCP_Log%          >>   %mylog%
	fc   .\SPE9\%myfolder%%OCP_SUMMARY%      .\SPE9\current\%OCP_SUMMARY%      >>   %mylog%
	fc   .\SPE9\%myfolder%%OCP_FastReview%   .\SPE9\current\%OCP_FastReview%   >>   %mylog%
)




REM Run CP test
if /i %CP_flag% == %myFlag% (
	echo.
	echo Test problem %CP_test%  >>	%mylog%
	%OCP%   %CP_test%   2>&1  | tee %CP_dir%%OCP_Log%
	
	md .\CP\%myfolder%
	
	copy %CP_dir%%OCP_Log%  			.\CP\%myfolder%%OCP_Log%
	copy %CP_dir%%OCP_SUMMARY% 			.\CP\%myfolder%%OCP_SUMMARY%
	copy %CP_dir%%OCP_FastReview% 		.\CP\%myfolder%%OCP_FastReview%
	
	fc   .\CP\%myfolder%%OCP_Log%          .\CP\current\%OCP_Log%          >   .\CP\CompareLog.out
	fc   .\CP\%myfolder%%OCP_SUMMARY%      .\CP\current\%OCP_SUMMARY%      >   .\CP\CompareSummary.out
	fc   .\CP\%myfolder%%OCP_FastReview%   .\CP\current\%OCP_FastReview%   >   .\CP\CompareFastRevive.out
	
	fc   .\CP\%myfolder%%OCP_Log%          .\CP\current\%OCP_Log%          >>   %mylog%
	fc   .\CP\%myfolder%%OCP_SUMMARY%      .\CP\current\%OCP_SUMMARY%      >>   %mylog%
	fc   .\CP\%myfolder%%OCP_FastReview%   .\CP\current\%OCP_FastReview%   >>   %mylog%
)




REM Run SPE3 test
if /i %SPE3_flag% == %myFlag% (
	echo.
	echo Test problem %SPE3_test%	>>	%mylog%
	%OCP%   %SPE3_test%   2>&1  | tee %SPE3_dir%%OCP_Log%
	
	md .\SPE3\%myfolder%
	
	copy %SPE3_dir%%OCP_Log%  				.\SPE3\%myfolder%%OCP_Log%
	copy %SPE3_dir%%OCP_SUMMARY% 			.\SPE3\%myfolder%%OCP_SUMMARY%
	copy %SPE3_dir%%OCP_FastReview% 		.\SPE3\%myfolder%%OCP_FastReview%
	
	fc   .\SPE3\%myfolder%%OCP_Log%          .\SPE3\current\%OCP_Log%          >   .\SPE3\CompareLog.out
	fc   .\SPE3\%myfolder%%OCP_SUMMARY%      .\SPE3\current\%OCP_SUMMARY%      >   .\SPE3\CompareSummary.out
	fc   .\SPE3\%myfolder%%OCP_FastReview%   .\SPE3\current\%OCP_FastReview%   >   .\SPE3\CompareFastRevive.out
	
	fc   .\SPE3\%myfolder%%OCP_Log%          .\SPE3\current\%OCP_Log%          >>   %mylog%
	fc   .\SPE3\%myfolder%%OCP_SUMMARY%      .\SPE3\current\%OCP_SUMMARY%      >>   %mylog%
	fc   .\SPE3\%myfolder%%OCP_FastReview%   .\SPE3\current\%OCP_FastReview%   >>   %mylog%
)



REM Run SPE5 test
if /i %SPE5_flag% == %myFlag% (
	echo.
	echo Test problem %SPE5_test%  >>	%mylog%
	%OCP%   %SPE5_test%   2>&1  | tee %SPE5_dir%%OCP_Log%
	
	md .\SPE5\%myfolder%
	
	copy 	%SPE5_dir%%OCP_Log%  				.\SPE5\%myfolder%%OCP_Log%
	copy 	%SPE5_dir%%OCP_SUMMARY% 			.\SPE5\%myfolder%%OCP_SUMMARY%
	copy 	%SPE5_dir%%OCP_FastReview% 			.\SPE5\%myfolder%%OCP_FastReview%
	
	fc   	.\SPE5\%myfolder%%OCP_Log%          .\SPE5\current\%OCP_Log%          >   .\SPE5\CompareLog.out
	fc   	.\SPE5\%myfolder%%OCP_SUMMARY%      .\SPE5\current\%OCP_SUMMARY%      >   .\SPE5\CompareSummary.out
	fc   	.\SPE5\%myfolder%%OCP_FastReview%   .\SPE5\current\%OCP_FastReview%   >   .\SPE5\CompareFastRevive.out
	
	fc   	.\SPE5\%myfolder%%OCP_Log%          .\SPE5\current\%OCP_Log%          >>   %mylog%
	fc   	.\SPE5\%myfolder%%OCP_SUMMARY%      .\SPE5\current\%OCP_SUMMARY%      >>   %mylog%
	fc   	.\SPE5\%myfolder%%OCP_FastReview%   .\SPE5\current\%OCP_FastReview%   >>   %mylog%
)




REM Run ZBCO2 test
if /i %ZBCO2_flag% == %myFlag% (
	echo.
	echo Test problem %ZBCO2_test%  >>	%mylog%
	%OCP%   %ZBCO2_test%	2>&1	|	tee %ZBCO2_dir%%OCP_Log%
	
	md 		.\ZBCO2\%myfolder%
	
	copy 	%ZBCO2_dir%%OCP_Log%				.\ZBCO2\%myfolder%%OCP_Log%
	copy 	%ZBCO2_dir%%OCP_SUMMARY% 			.\ZBCO2\%myfolder%%OCP_SUMMARY%
	copy 	%ZBCO2_dir%%OCP_FastReview% 		.\ZBCO2\%myfolder%%OCP_FastReview%
	
	fc   	.\ZBCO2\%myfolder%%OCP_Log%			.\ZBCO2\current\%OCP_Log%       	> 	.\ZBCO2\CompareLog.out
	fc   	.\ZBCO2\%myfolder%%OCP_SUMMARY% 	.\ZBCO2\current\%OCP_SUMMARY%    	> 	.\ZBCO2\CompareSummary.out
	fc   	.\ZBCO2\%myfolder%%OCP_FastReview% 	.\ZBCO2\current\%OCP_FastReview%	>	.\ZBCO2\CompareFastRevive.out
	
	fc   	.\ZBCO2\%myfolder%%OCP_Log%			.\ZBCO2\current\%OCP_Log%       	>>   %mylog%
	fc   	.\ZBCO2\%myfolder%%OCP_SUMMARY% 	.\ZBCO2\current\%OCP_SUMMARY%    	>>   %mylog%
	fc   	.\ZBCO2\%myfolder%%OCP_FastReview% 	.\ZBCO2\current\%OCP_FastReview%	>>   %mylog%
)





set /p Success="If accept this result? if then cover <current> directory with it. [Y/N]: "
if /I %Success% == %myFlag% (
	if /I %SPE1a_flag% == %myFlag% 	( xcopy .\SPE1a\%myfolder%  .\SPE1a\current\  /Y )
	if /I %SPE1b_flag% == %myFlag% 	( xcopy .\SPE1b\%myfolder%  .\SPE1b\current\  /Y )
	if /I %SPE9_flag% == %myFlag% 	( xcopy .\SPE9\%myfolder%  	.\SPE9\current\  /Y )
	if /I %CP_flag% == %myFlag% 	( xcopy .\CP\%myfolder%  	.\CP\current\  /Y )
	if /I %SPE3_flag% == %myFlag% 	( xcopy .\SPE3\%myfolder%  	.\SPE3\current\  /Y )
	if /I %SPE5_flag% == %myFlag% 	( xcopy .\SPE5\%myfolder%  	.\SPE5\current\  /Y )
	if /I %ZBCO2_flag% == %myFlag% 	( xcopy .\ZBCO2\%myfolder%  .\ZBCO2\current\  /Y )
	echo "<current> directories have been updated!"
) else (
	echo "<current> directories have not been updated!"
)