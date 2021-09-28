;*************************************************************
;Aim: cal monthly q2 from daily q2
;Mission:
;   1. Read data from NCEP-DOE data
;   2. Post-post process
;     2.1 Calculate q2
;   3. Output to yearly files
;                       Zhenning Li
;                       2015-05-06
;Last Edited:
;   2015-03-01  Create
;   2015-03-03  Change to fit HEAT Experiment
;   2015-04-05  Change to fit Team Server
;   2015-04-08  Change to fit ECMWF MONLTM MEAN 
;   2015-04-18  Change to caculate q2
;   2015-04-18  Change to caculate daily q2 in NCEP-DOE
;   2015-05-06  Change to Monthly q2
;*************************************************************

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin


;global variables
    ;Define parameters

    ;for directory 
    g_strt_yr=1979
    g_lst_yr=2016
    g_data_path = "/home/yangsong3/data-observation/ERA-interim/pressure/q2_1979-2016/"
    g_data_out_path = "/home/yangsong3/data-observation/ERA-interim/monthly/q2/"




;*************************For F_Compset Ctrl***********************************

; ==============================================================
; Mission 1. Open the file: Read only the user specified period
; ==============================================================
    do nyr = g_strt_yr,g_lst_yr
        print(nyr+"Started!")
        ctrl_in     = addfile (g_data_path+"q2."+nyr+".nc", "r")   ; note the "s" of addfile
        var1 = ctrl_in->q2
        mon_var1 = calculate_monthly_values(var1,"ave",0,False)
        mon_var1@long_name = "Monthly q2 Calculated from daily on Pressure Levels"
        system("rm "+g_data_out_path+"q2."+nyr+".mon.mean.nc")
        ncdf = addfile(g_data_out_path+"q2."+nyr+".mon.mean.nc" ,"c")  ; open output netCDF file
        ; output variables directly
        ncdf->q2=mon_var1
        print(nyr+"Done!")
        delete([/var1, mon_var1/])
    end do
end

