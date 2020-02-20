@echo off
setlocal EnableDelayedExpansion

rem Create vector with names of days
set /a i=0
for %%d in (31 28 31 30 31 30 31 31 30 31 30 31) do (
   set /A i=i+1
   set day[!i!]=%%d
)	

set /a i=0
for %%d in (01 02 03 04 05 06 07 08 09 10 11 12) do (
   set /A i=i+1
   set month[!i!]=%%d
)	

set /a l_limit=1
set /a u_limit=l_limit + 20
set /a i=1
for /l %%m in (1,1,12) do (
	for /l %%d in (1,1,!day[%%m]!) do (
		if !i! GEQ %l_limit% (
			if %%d LSS 10 (set day_text=0%%d)
			if %%d GEQ 10 (set day_text=%%d)
			(echo retrieve,
			echo   class   = od,
			echo   stream  = enfo,
			echo   expver  = 1,
			echo   date    = 2018!month[%%m]!!day_text!,
			echo   time    = 00/06/12/18,
			echo   step    = 0/1/2/3/4/5/6/7/8/9/10/11/12/13/14/15/16/17/18/19/20/21/22/23/24,
			echo   type    = cf,
			echo   levtype = sfc,
			echo   param   = ssrd,
			echo   grid    = 0.2/0.2,
			echo   area	  = 48.4/-116.2/34.2/-77.8,
			echo   target  = "north_america_ctrl_!month[%%m]!!day_text!.grib") > north_america_ctrl_!month[%%m]!!day_text!.req

			start python C:\Users\kdoubled\mars\bin\mars north_america_ctrl_!month[%%m]!!day_text!.req

			rem rm north_america_ctrl_!month[%%m]!!day_text!.req
		)

		set /A i=i+1
		if !i! GEQ %u_limit% (exit /b)
	)
) 


