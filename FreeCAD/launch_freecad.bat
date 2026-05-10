@echo off
REM Launch FreeCAD 0.20 with clean PATH and software OpenGL fallback.
REM Error 87 on Windows 11 is typically caused by OpenGL driver incompatibility
REM with FreeCAD's bundled Qt, or by Anaconda DLLs shadowing FreeCAD's own.

set FREECAD_DIR=C:\Program Files\FreeCAD 0.20\bin

REM Strip PATH to essentials + FreeCAD only
set PATH=%FREECAD_DIR%;%SystemRoot%\system32;%SystemRoot%;%SystemRoot%\System32\Wbem;%SystemRoot%\System32\WindowsPowerShell\v1.0

REM Force Qt to use software OpenGL (bypasses GPU driver issues on Windows 11)
set QT_OPENGL=software

REM Fallback: also try ANGLE (DirectX-based) if software still fails — swap the
REM comment below and re-run if FreeCAD still crashes.
REM set QT_OPENGL=angle

start "" "%FREECAD_DIR%\FreeCAD.exe" %*