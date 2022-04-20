;;*************************************************
;; main program 
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

begin
    ; startmon = 101
    ; endmon = 89712
    ; startyear = 1
    ; endyear = 897

    startmon = 101
    endmon = 99712
    startyear = 1
    endyear = 997
;;;读入海温数据
    diri="/home/yangsong3/CMIP6/1.0x1.0/CESM2/piControl/Omon/"
    fils1=systemfunc("ls "+diri+"tos_Omon_CESM2_piControl_r1i1p1f1_gn_0{0..9}"+"*")
  
    f1=addfiles(fils1, "r")
    date:=cd_calendar(f1[:]->time, 1)

    lat = f1[0]->lat
    lon = f1[0]->lon

    timeind:=ind(date.le.endmon+200.and.date.ge.startmon)

    sst=lonFlip(f1[:]->tos(timeind,:,:))
    printVarSummary(sst)
    
    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)

    ;;detrend 
    ssta = (/dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)/)
    ;;;
    delete(date)
    delete(timeind)
    
    nino34 =  dim_avg_n_Wrap(ssta(:,{5:-5},{190:240}), (/1,2/))
    nino34_djf = month_to_season(nino34, "DJF")

    printVarSummary(nino34_djf)

    nino34_djf(0:dimsizes(nino34_djf)-2) = nino34_djf(1:) 
    nino34_djf(dimsizes(nino34_djf)-1) = nino34(dimsizes(nino34)-1) 

    nino34_12_r =  split_12(nino34)

    nino34_djf = (/dim_avg_n_Wrap(nino34_12_r(5:8,:), 0)/)


    std_nino = dim_stddev(nino34_djf)

    sst_total_year = ispan(startyear, endyear+2, 1) ;;;0-199

    ; nino_year = sst_total_year(ind(nino34_djf.ge.std_nino))
    ; lanina_year = sst_total_year(ind(nino34_djf.le.-1*std_nino))

    nino_year = sst_total_year(ind(nino34_djf.gt.0.43*std_nino))
    lanina_year = sst_total_year(ind(nino34_djf.lt.-0.43*std_nino))
    normal_year = sst_total_year(ind(nino34_djf.gt.-0.43*std_nino.and.nino34_djf.lt.0.43*std_nino))

    print("************** finish SST data*****************")
;;;


;;读入ua
    diri2 = "/home/yangsong3/CMIP6/2.5x2.5/CESM2/piControl/"
    fils2 = systemfunc("ls "+diri2+"ua_Amon_CESM2_piControl_r1i1p1f1_gn_0{0..9}"+"*")
  
    f2 = addfiles(fils2, "r")
    date:=cd_calendar(f2[:]->time, 1)

    lat := f2[0]->lat
    lon := f2[0]->lon

    timeind:=ind(date.le.endmon.and.date.ge.startmon)

    u_850 = lonFlip(f2[:]->ua(timeind,{85000},:,:))
    printVarSummary(u_850)
 
;;;;


    u_850 = (/rmMonAnnCycTLL(u_850)/)
    u_850 = (/dtrend_msg_n(ispan(1,dimsizes(u_850&time),1),u_850,False,True,0)/)


;;;;;EAESM指数 : 850hPA ZONAL wind diff 110-140 22.5-32.5 - 90-130,5-15
    DMI_ctrl = u_850(:,1,1)

    DMI_ctrl  =  (/dim_avg_n_Wrap(u_850(:,{22.5:32.5},{110:140}), (/1,2/)) - \
                  dim_avg_n_Wrap(u_850(:,{5:15},{90:130}), (/1,2/))/)

    DMI_ctrl_JJA = month_to_season(DMI_ctrl, "JJA")  ;;;jja 少一年，从25601开始
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


corr = escorc_n(DMI_ctrl_JJA, ISM_JJA, 0, 0)

print(corr)

end 