echo off
set LOCALHOST=%COMPUTERNAME%
set KILL_CMD="F:\ANASYS\ANSYS Inc\v202\fluent/ntbin/win64/winkill.exe"

"F:\ANASYS\ANSYS Inc\v202\fluent\ntbin\win64\tell.exe" wenzheshang 3090 CLEANUP_EXITING
if /i "%LOCALHOST%"=="wenzheshang" (%KILL_CMD% 101632) 
if /i "%LOCALHOST%"=="wenzheshang" (%KILL_CMD% 105480) 
if /i "%LOCALHOST%"=="wenzheshang" (%KILL_CMD% 99016)
del "f:\Thinking\program\BELab1.0\Workdata\Fluent_Python\2023-03-05_19-08\cleanup-fluent-wenzheshang-105480.bat"
