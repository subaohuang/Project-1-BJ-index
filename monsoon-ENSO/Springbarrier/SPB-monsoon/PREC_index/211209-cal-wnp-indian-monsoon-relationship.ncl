;;*************************************************
;; main program 
;;*************************************************************
begin
    startmon = 190001
    endmon = 201912

    startyear = 1900
    endyear = 2019
;;
;;AISMR DATA 

    ncols = 18 ;;;":YEAR JAN FEB   MAR   APR   MAY   JUN   JUL   AUG   SEP   OCT   NOV   DEC    JF   MAM  JJAS   OND   ANN"
    nrows=146
    ; data = readAsciiTable("All-indian-monsoon-rainfall-index.txt", ncol, "float", 2)
    data = asciiread("All-indian-monsoon-rainfall-index.txt", (/nrows,ncols/), "float")
    jjas = data(:,15) ;;sum of jjas precipitation  1871-2016
    
    jjas!0 = "time"
    jjas&time = ispan(1871,2016,1)

    pre_jjas = jjas({startyear:2012}) 

    pre_jjas = pre_jjas*0.1 ;;convert to mm
    clm = dim_avg_n_Wrap(pre_jjas, 0)
    pre_jjasa = dim_rmvmean(pre_jjas)
    pre_jjasa = (/dtrend_msg_n(ispan(1,2012-startyear+1,1),pre_jjasa,False,True,0)/)
;;;
    ISMR = pre_jjasa
;;;
;;;读入ERA wind data

   ; fils1 = "/home/yangsong3/data-observation/ERA5-monthly/pressure/u_mon_1x1_1979_2020.nc"
   fils1 = "/home/yangsong3/data-observation/20thCRV2/uwnd.mon.mean.nc"
   f1 = addfile(fils1, "r")
   time1 = f1->time
   date1=cd_calendar(time1, 1)

   ; endmon = endmon - 100  ;;; 2019 is normal year 
   timeind1=ind(date1.ge.startmon.and.date1.le.endmon)
 
  
;;;;U,V
    ; u_200  = short2flt(f1->u(timeind1,{200},:,:))
    ; u_850  = short2flt(f1->u(timeind1,{850},:,:))

    u_200  = f1->uwnd(timeind1,{200},:,:)
    u_850  = f1->uwnd(timeind1,{850},:,:)

    u_850 = (/rmMonAnnCycTLL(u_850)/)
    u_200 = (/rmMonAnnCycTLL(u_200)/)
;;;

   ; fils2 = "/home/yangsong3/data-observation/ERA5-monthly/pressure/v_mon_1x1_1979_2020.nc"
   ; f2 = addfile(fils2, "r")
   ; time2 = f2->time
   ; date2 = cd_calendar(time2, 1)
   ; timeind2 = ind(date2.ge.197901.and.date2.le.201912)

;;;;U,V
    ; v_200  = short2flt(f2->v(timeind2,{200},:,:))
    ; v_850  = short2flt(f2->v(timeind2,{850},:,:))
    ; v_850 = (/rmMonAnnCycTLL(v_850)/)
    ; v_200 = (/rmMonAnnCycTLL(v_200)/)

     u_850 = (/dtrend_msg_n(ispan(1,dimsizes(u_850&time),1),u_850,False,True,0)/)
;;;


;;;;;EAESM指数 : 850hPA ZONAL wind diff 110-140 22.5-32.5 - 90-130,5-15
    DMI_ctrl = u_850(:,1,1)

    DMI_ctrl  =  (/dim_avg_n_Wrap(u_850(:,{22.5:32.5},{110:140}), (/1,2/)) - \
                  dim_avg_n_Wrap(u_850(:,{5:15},{90:130}), (/1,2/))/)

    DMI_ctrl_JJA = month_to_season(DMI_ctrl, "JJA")  ;;;jja 少一年，从25601开始
;;;


;;;ISM index  850hPA ZONAL wind diff 5-15N 40-80E - 20-30N,70-90E
    ISM = u_850(:,1,1)

    ISM = (/dim_avg_n_Wrap(u_850(:,{5:15},{40:80}), (/1,2/)) - \
                  dim_avg_n_Wrap(u_850(:,{20:30},{70:90}), (/1,2/))/)

    ISM_JJA = month_to_season(ISM, "JJA")  ;;;jja 少一年，从25601开始
;;;;


; corr = escorc_n(DMI_ctrl_JJA, ISM_JJA, 0, 0)
corr = escorc_n(DMI_ctrl_JJA, ISMR, 0, 0)
; corr = escorc_n(ISM_JJA, ISMR, 0, 0)


print(corr)

end 