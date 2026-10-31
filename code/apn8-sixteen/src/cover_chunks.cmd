@echo off
rem cover_chunks.cmd <dir>: run runs\bin\cover3.exe on every <dir>\c*.txt in turn, output <dir>\c*.out; a chunk whose
rem .out already ends with cover3's "done:" line is skipped, so an evicted job resumes at the chunk it lost.
rem Rowen (Claude Opus 5.5), 2026-10-30, for i = 9 (8 rank-4).
cd /d %~dp0..
for %%f in (%1\c*.txt) do call :one "%%~dpnf"
exit /b 0
:one
findstr /b /c:"done: covers" "%~1.out" >nul 2>nul && exit /b 0
runs\bin\cover3.exe "%~1.txt" > "%~1.out" 2>&1
if errorlevel 1 exit 1
exit /b 0
