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
set /a u_limit=l_limit + 30
set /a i=1
for /l %%m in (1,1,12) do (
	for /l %%d in (1,1,!day[%%m]!) do (
		if !i! GEQ %l_limit% (
			if %%d LSS 10 (set day_text=0%%d)
			if %%d GEQ 10 (set day_text=%%d)

			start grib_to_netcdf -o Raw\north_america_!month[%%m]!!day_text!.nc -T -R 20180101 "GRIB files\north_america_!month[%%m]!!day_text!.grib"

		)

		set /A i=i+1
		if !i! GEQ %u_limit% (exit /b)
	)
) 


