echo off
set LOCALHOST=%COMPUTERNAME%
set KILL_CMD="F:\ANASYS\ANSYS Inc\v202\fluent/ntbin/win64/winkill.exe"

"F:\ANASYS\ANSYS Inc\v202\fluent\ntbin\win64\tell.exe" wenzheshang 23615 CLEANUP_EXITING
if /i "%LOCALHOST%"=="wenzheshang" (%KILL_CMD% 102080) 
if /i "%LOCALHOST%"=="wenzheshang" (%KILL_CMD% 102916) 
if /i "%LOCALHOST%"=="wenzheshang" (%KILL_CMD% 107468)
del "f:\Thinking\program\BELab1.0\Workdata\Fluent_Python\2023-03-05_19-38\cleanup-fluent-wenzheshang-102916.bat"
