;;*************************************************
;; draw circulation pattern during monsoon year 
;;*************************************************************
function split_12(data)
begin
    dim = dimsizes(data) 
    yearnum = dim(0)/12

    data_12 = new((/12,yearnum/),float)
    do i=0,11   
    data_12(i,:) = data(i::12)
    end do 
return data_12
end 
function bw_bandpass(var)
 begin
    ca    = 9.0*12.0        ; band start (longer period)
    cb    = 4.0        ; band end

    fca   = 1.0/ca      ; 'left'  frequency
    fcb   = 1.0/cb      ; 'right' frequency

    opt   = False       ; use default options (time series of filtered
                        ; values will be returned)
    var = (/bw_bandpass_filter (var,fca,fcb,opt,0)/)  
 return var
 end 
;;*************************************************
;; main program 
;;*************************************************************
begin
    startmon = 101
    endmon = 99912
    startyear = 1
    endyear = 999

    ; startmon = 101
    ; endmon = 120012
    ; startyear = 1
    ; endyear = 1200

    lat1 = -40
    lat2 = 40
    lon1 = -180
    lon2 = 180

; ;;;读入海温数据
    diri="/home/ys17-23/CMIP6/2.5x2.5/CESM2/piControl/Omon/"
    
    fils1=systemfunc("ls "+diri+"thetao_Omon_CESM2_piControl_r1i1p1f1_0-300mean.nc")

    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)

    lat := f1->lat
    lon := f1->lon

    timeind:=ind(date.le.endmon.and.date.ge.startmon)

    sst=lonFlip(f1->thetao(timeind,{lat1:lat2},{lon1:lon2}))
    printVarSummary(sst)
    
    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)

    ;;detrend 
    ssta = (/dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)/)
    ;;;
    delete(date)
    delete(timeind)

    ssta = (/bw_bandpass(ssta)/)
   
    print("************** finish SST data*****************")
;;;;; 

    system("rm highpass_bw_thetao_0-300mean_picontrol.nc")
    ncout = addfile("highpass_bw_thetao_0-300mean_picontrol.nc", "c")
    ncout->thetao = ssta
;;;;;
end 
