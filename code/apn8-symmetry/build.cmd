@echo off
rem build.cmd [program ...]   Windows, MSVC (Visual Studio 2022 or its Build Tools, C++ workload).
rem Run from an "x64 Native Tools Command Prompt for VS 2022", or call vcvars64.bat first.
rem Builds every program (or the ones named) into bin\ (or %BIN% if set).
setlocal
cd /d "%~dp0"
where cl >nul 2>nul || (echo cl.exe not found: open an x64 Native Tools Command Prompt, or call vcvars64.bat & exit /b 1)
if "%BIN%"=="" set BIN=bin
if not exist "%BIN%" mkdir "%BIN%"
rem review_odd_engine includes c3k4.cpp's Engine (lines 14-205) as a separate file
powershell -NoProfile -Command "(Get-Content src\c3k4.cpp)[13..204] | Set-Content -Encoding ascii src\c3k4_engine.inc" || exit /b 1
set PROGS=%*
if "%PROGS%"=="" set PROGS=orbit5 fix5check symp c3k4 inv2 aut8 classify8 review_engine review_inv2 review_odd_brute review_odd_engine
for %%p in (%PROGS%) do (
    echo %%p
    cl /nologo /O2 /openmp /EHsc /std:c++17 /arch:AVX2 src\%%p.cpp /Fe"%BIN%\%%p.exe" /Fo"%BIN%\%%p.obj" >"%BIN%\%%p.buildlog" || (type "%BIN%\%%p.buildlog" & exit /b 1)
)
