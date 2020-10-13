
; function  draw_corr(mode,wks)
begin
; ipindex1=asciiread("/mnt/f/share/MRE_msf_mean_cold_locationxi_index.ascii", 55, "float")
; ipindex2=asciiread("/mnt/f/share/MRE_msf_mean_cold_locationdong_index.ascii", 55, "float")

ipindex1=asciiread("/mnt/f/share/JRA-55_msf_cold_locationxi_index.ascii", 55, "float")
ipindex2=asciiread("/mnt/f/share/JRA-55_msf_cold_locationdong_index.ascii", 55, "float")

ipindex3=ipindex2- ipindex1


ipindex1=dim_rmvmean(ipindex1)
ipindex2=dim_rmvmean(ipindex2)
ipindex3=dim_rmvmean(ipindex3)

;;;读入HadiSST海温数据
  diri="/mnt/f/share/ersst/"
  fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
  f1=addfile(fils1, "r")
  date=cd_calendar(f1->time, 1)

  timeind=ind(date.le.201312.and.date.ge.195801)
  

  sst=lonFlip(f1->sst(timeind,{30:-30},:))
  ssta=(/rmMonAnnCycTLL(sst)/)
  copy_VarCoords(sst, ssta)
  ssta_de=dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)
  copy_VarCoords(ssta, ssta_de)
  ssta_son=month_to_season(ssta_de, "SON")
  ssta_djf=month_to_season(ssta_de, "DJF")

  ssta_s=(/ssta_son(0:54,:,:)+ssta_djf(1:,:,:)/)*0.5
  copy_VarCoords(ssta_son(0:54,:,:), ssta_s)
    delete(date)
    delete(timeind)
    copy_VarCoords(ssta_son(0:54,:,:), ssta_s)

    emi=dim_avg_n(ssta_s(:,{10:-10},{165:220}),(/1,2/))-0.5*dim_avg_n(ssta_s(:,{20:-10},{125:145}), (/1,2/))-0.5*dim_avg_n(ssta_s(:,{5:-15},{250:290}),(/1,2/))
    
    emi_center=dim_avg_n(ssta_s(:,{10:-10},{165:220}),(/1,2/))

    emi_std=dim_standardize_n(emi, 1, 0)
    std_emi=dim_stddev_n(emi, 0)

    nino3=dim_avg_n(ssta_s(:,{5:-5},{210:270}), (/1,2/))
    nino3_std=dim_standardize_n(nino3, 1, 0)
    std_nino3=dim_stddev_n(nino3, 0)
     
    nino3_4=dim_avg_n(ssta_s(:,{5:-5},{190:240}), (/1,2/))  
;;; 
    x=ipindex1
    y=ipindex2
   
     mxlag    = 9
     x_Lead_y = esccr(x,y,mxlag)
     y_Lead_x = esccr(y,x,mxlag)    ; switch the order of the series

     ccr = new ( 2*mxlag+1, float)    
     ccr(0:mxlag-1) = y_Lead_x(1:mxlag:-1)  ; "negative lag", -1 reverses order
     ccr(mxlag:)    = x_Lead_y(0:mxlag)     ; "positive lag"
  
wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/JRA-55&hadisst超前滞后相关")

    res=True
    plot=gsn_csm_xy(wks, fspan(-9,9,19), ccr, res)
end