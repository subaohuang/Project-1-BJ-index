;;计算偏相关，id1和id2,去掉id3的影响
function partial_corr(id1,id2,id3)
begin
    r12 = escorc_n(id1, id2, 0, 0)
    r13 = escorc_n(id1, id3, 0, 0)
    r23 = escorc_n(id2, id3, 0, 0)

    r12_3 = (/r12 - r13*r23/)/(sqrt((1 - r13^2)*(1 - r23^2)))
return r12_3
end 
begin
fils1  =  "/home/ys17-19/lsh/data/era/era_pre_global.nc"
; fils1 = "/home/yangsong3/data-observation/NCEP_NCAR/monthly/uwnd.mon.mean.nc"
; fils1  =  "/home/ys17-19/lsh/data/era/era_sur_global.nc"


    f1 = addfile(fils1, "r")
    date  =  cd_calendar(f1->time, 1)

    timeind  =  ind(date.le.201212.and.date.ge.197901)
    delete(date)


    u_heat = short2flt(f1->u(timeind,{200},:,:))
    ; u_heat = short2flt(f1->u10(timeind,:,:))

    ; v_heat = short2flt(f1->v(timeind,{1000},:,:))
    ; u_heat = f1->uwnd(timeind,{1000},:,:)

    u_heat = (/dtrend_msg_n(ispan(1,dimsizes(u_heat&time),1),u_heat,False,True,0)/)
    delete(timeind)
    
   
   ; ;;;超前滞后相关
   ;  x = dim_avg_n_Wrap(u_heat(:,{-5:5},{150:210}), (/1,2/))
   ;  y = dim_avg_n_Wrap(u_heat(:,{-5:5},{70:100}), (/1,2/))

   ;  mxlag    = 12
   ;  x_Lead_y = esccr(x,y,mxlag)
   ;  y_Lead_x = esccr(y,x,mxlag)    ; switch the order of the series

   ;  ccr = new ( 2*mxlag+1, float)    
   ;  ccr(0:mxlag-1) = y_Lead_x(1:mxlag:-1)  ; "negative lag", -1 reverses order
   ;  ccr(mxlag:)    = x_Lead_y(0:mxlag)     ; "positive lag"

   ;  wks=gsn_open_wks("pdf", "NCEP_NCAR-IWC与PWC的超前滞后相关")

   ;  res=True
   ;  plot=gsn_csm_xy(wks, fspan(-12,12,25), ccr, res)

   ;  delete(plot)
   ;  delete(wks)
   ;  exit()
   ;  ;;;

    ; ;;;
    ; season = (/"MAM","JJA","SON","DJF"/)
    ; dim = dimsizes(u_heat)
    ; u_heat_mon = new((/12,34,dim(1),dim(2)/),typeof(u_heat))
    ; do i = 0,11,1
    ; u_heat_mon(i,:,:,:) = u_heat(i:dim(0)-1:12,:,:
    ; end do
    ; ;;;计算每个季节的相关系数
    ; do i = 0,3   
    ; u_heat_mam = month_to_season(u_heat, season(i))

    ; v_heat_mam = month_to_season(v_heat, season(i))


    ; u_heat_mam_pa = dim_avg_n_Wrap(u_heat_mam(:,{-5:5},{150:210}), (/1,2/))

    ; u_heat_mam_id = dim_avg_n_Wrap(u_heat_mam(:,{-5:5},{70:100}), (/1,2/))

    ; corr_heat = escorc_n(u_heat_mam_pa, u_heat_mam_id, 0, 0)
    ; print("corr_heat="+corr_heat)
    ; end do 
    ;;;

    ;;;冬季的partial corrlation
      ;;读入Nino3.4
      f2=addfile("/home/ys17-19/lsh/data/sst/HadISST_sst.nc", "r")
      date=cd_calendar(f2->time, 1)
      timeind=ind(date.ge.197901.and.date.le.201212)
      sst=lonFlip(f2->sst(timeind,{20:-20},:))
      ssta=(/rmMonAnnCycTLL(sst)/)
      copy_VarCoords(sst, ssta)
      ssta_de=dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)
      copy_VarCoords(ssta, ssta_de)
      nino3_4=dim_avg_n_Wrap(ssta_de(:,{5:-5},{190:240}), (/1,2/))   
      nino3_4_anu=(/month_to_season(nino3_4, "DJF")/)
      ;;

      ;;计算偏相关
        u_heat_djf = month_to_season(u_heat, "DJF")

        pa_djf = dim_avg_n_Wrap(u_heat_djf(:,{-5:5},{150:210}), (/1,2/))

        id_djf = dim_avg_n_Wrap(u_heat_djf(:,{-5:5},{70:100}), (/1,2/))
       
        corr_djf = partial_corr(pa_djf,id_djf,nino3_4_anu)

        print(corr_djf)
        exit()
       ;;




    ;;;计算每个月的相关系数

    u_heat_mon_pa = dim_avg_n_Wrap(u_heat_mon(:,:,{-5:5},{150:210}), (/2,3/))

    u_heat_mon_id = dim_avg_n_Wrap(u_heat_mon(:,:,{-5:5},{70:100}), (/2,3/))

    corr_heat_mon = escorc_n(u_heat_mon_id, u_heat_mon_pa, 1, 1)

    printVarSummary(corr_heat_mon)
    ;;;


    ;;;画图
    
    wks = gsn_open_wks("pdf", "era-mon-corr-iwc-pwc")
    res =True
    res@tmXBMode = "Manual"
    res@tmXBTickStartF = 1
    res@tmXBTickEndF = 12
    res@xyMonoMarkLineMode = True
    res@xyMarkLineMode  = (/"MarkLines"/) 
    
    res@xyLineColors      = (/"black"/)          ; change line color
    
    
    plot = gsn_csm_xy(wks, ispan(1, 12, 1), corr_heat_mon, res)


end    





;;;
