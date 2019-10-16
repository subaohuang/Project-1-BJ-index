;;;************************************************************
;; investigate and draw the annual cycle of WC using 10-m wind ,slp, msf,
;；随时间的变化情况
;; linshuheng 2019.10.14
;;;*********************************************************

function annual_cycle(u10_inten)
local u10_inten_month, nmos,nmo
begin
    u10_inten_month = new(12,typeof(u10_inten))
    
    nmos = 12
    ntim = 456
    do nmo=0,nmos-1
    u10_inten_month(nmo) = dim_avg( u10_inten(nmo:ntim-1:nmos) )
    end do
return u10_inten_month
end

begin 
;;=========read location data=================
   ;  locationxi = asciiread("/mnt/f/share/rea_msf/era_1.5wlocationxi.ascii", 456, "float")
   ;  locationxi = (/dtrend_msg_n(ispan(1,456,1),locationxi,False,True,0)/)
   ;  locationxi!0 = "time"
   ;  locationxi&time = ispan(1,456,1)



   ;  locationxi_month = new(12,typeof(locationxi))
   ;  nmos = 12
   ;  ntim = 456

   ; do nmo=0,nmos-1
   ;    locationxi_month(nmo) = dim_avg( locationxi(nmo:ntim-1:nmos) )
   ; end do
;;;

;;=======read sst and calculate correlation==============
    ; f2=addfile("/mnt/f/share/ersst/HadISST_sst.nc", "r")
    ; time = f2->time
    ; date = cd_calendar(time, 1)
    ; timeind=ind(date.ge.197901.and.date.le.201612)
    ; sst=lonFlip(f2->sst(timeind,{30:-30},:))

    ; sst = (/dtrend_msg_n(ispan(1,456,1),sst,False,True,0)/)
    
    ; sst_clm   = clmMonTLL(sst)
    ; nino3 = dim_avg_n_Wrap( sst_clm(:,{-5:5},{210:270}) , (/1,2/))
 
    ; ; nino3  =  dim_avg_n_Wrap(sst_clm(:,{-5:5},{80:160}), (/1,2/))-dim_avg_n_Wrap(sst_clm(:,{-5:5},{100:280}), (/1,2/)) 
    ; delete([/time,date,timeind/])
;;;


;;=========read msf data=============================
    fils1 = systemfunc("ls /home/ys17-19/lsh/data/rea_msf/msf_time_era_global.nc ")
    f1 = addfile(fils1, "r")

    time = f1->time
    date = cd_calendar(time, 1)

    tind = ind(date.le.201612.and.date.ge.197901)
    msf_time = f1->msf_time(:,:,tind)

   ;;;去趋势
    msf_time = (/dtrend_msg_n(ispan(1,dimsizes(msf_time&time),1),msf_time,False,True,2)/)
   ;;;
   
    msf_inten = dim_avg_n_Wrap(msf_time(:,{150:240},:), (/0,1/))
    msf_inten_month = annual_cycle(msf_inten)


    ud_time = f1->ud_time({1000},{150:240},tind)
    ud_time = (/dtrend_msg_n(ispan(1,dimsizes(ud_time&time),1),ud_time,False,True,1)/)

    ud_inten = dim_avg_n_Wrap(ud_time, 0)
    ud_inten = (/-1*ud_inten/)
    ud_inten_month = annual_cycle(ud_inten)
    
    w_time = f1->w_time({500},:,tind)
    w_time = (/dtrend_msg_n(ispan(1,dimsizes(w_time&time),1),w_time,False,True,1)/)
    
    ;;intensity of PWC defined by wang et al using (5S– 5N, 160–120W) and the equatorial western Pacific
    ; (5S–5N, 120–160E) as an index
    w_inten = dim_avg_n_Wrap(w_time({120:160},:),0) - \            
                  dim_avg_n_Wrap(w_time({200:240},:),0)
    
    w_inten_month = annual_cycle(w_inten)
    printVarSummary(ud_inten_month)

    delete([/time,date,tind/])

;;---------------------------------------
;;read  wind and slp data
;----------------------------------------
    fils = systemfunc("ls /home/ys17-19/lsh/data/era/sur_1.5/era_interim_sur1.5_"+"*"+".nc")
    f3 = addfiles(fils, "r")

    time=f3[:]->time
    date=cd_calendar(time, 1)
    tind=ind(date.le.201612.and.date.ge.197901)
    u10 = short2flt(f3[:]->u10(tind,:,:))
    slp = short2flt(f3[:]->msl(tind,:,:))
    printVarSummary(u10)     

    u10 = (/dtrend_msg_n(ispan(1,dimsizes(u10&time),1),u10,False,True,0)/)
    slp = (/dtrend_msg_n(ispan(1,dimsizes(slp&time),1),slp,False,True,0)/)

    u10_inten = dim_avg_n_Wrap(u10(:,{-5:5},{150:210}),(/1,2/))
    ; u10_inten = dim_avg_n_Wrap(u10(:,{-5:5},{120:170}),(/1,2/))
    u10_inten = (/-1*u10_inten/)
    u10_inten_month = annual_cycle(u10_inten)

    slp_inten = dim_avg_n_Wrap(slp(:,{-5:5},{200:280}),(/1,2/)) - \
                  dim_avg_n_Wrap(slp(:,{-5:5},{80:160}),(/1,2/))
    copy_VarCoords(slp(:,1,1), slp_inten)
                
    slp_inten_month = annual_cycle(slp_inten)
;;;

;;;===============read soi data=======================================
   ;;read in data ncols = 9 nrows = 3 
   ;;这个soi数据的格式是https://crudata.uea.ac.uk/cru/data/soi/ 上说的，year, 12 * monthly value，annual value
   ;;data time span from 186601-201912
   ; soi_r = asciiread ("/home/ys17-19/lsh/data/soi.dat", (/154,14/), "float")
   dar_r = asciiread ("/home/ys17-19/lsh/data/soi_dar.dat", (/165,13/), "float")
   tah_r = asciiread ("/home/ys17-19/lsh/data/soi_tah.dat", (/165,13/), "float")
   soi_r = (/tah_r - dar_r/)/10
   

   soi_2d = soi_r(:,1:12)
   soi = ndtooned(soi_2d)

   ;;;select time  
   tistart = (1979-1866)*12
   tiend = (2016-1866)*12+11
   ;tistart = (1979-1855)*12
   ;tiend = (2016-1855)*12+11
   soi_inten = soi(tistart:tiend)
   soi_inten_month = annual_cycle(soi_inten)
       


;;===========标准化强度指数==========================

    msf_inten_month_std = dim_standardize(msf_inten_month, 1)

    u10_inten_month_std = dim_standardize(u10_inten_month, 1)

    slp_inten_month_std = dim_standardize(slp_inten_month, 1)

    ud_inten_month_std = dim_standardize(ud_inten_month, 1)

    soi_inten_month_std = dim_standardize(soi_inten_month, 1)

    w_inten_month_std = dim_standardize(w_inten_month, 1)    
   
;;;

    
 
;; corr 
   
; ;*************************************************** 超前滞后相关
    mxlag    = 12
    x_Lead_y = esccr(u10_inten,slp_inten, mxlag)
    y_Lead_x = esccr(slp_inten,u10_inten, mxlag)    ; switch the order of the series

    ccr = new ( 2*mxlag+1, double)    
    ccr(0:mxlag-1) = y_Lead_x(1:mxlag:-1)  ; "negative lag", -1 reverses order
    ccr(mxlag:)    = x_Lead_y(0:mxlag)     ; "positive lag"
    x1 = ispan(-12,12,1)     
    wks1=gsn_open_wks("pdf", "u10定义的沃克环流强度与dslp的超前滞后相关")
    res=True
    res@gsnXRefLine = 0
    plot1=gsn_csm_xy(wks1,x1,ccr,res)
;;;

x = ispan(1,12,1)

;;---------------------------------------
;;draw_piture
;----------------------------------------

;;;
    wks=gsn_open_wks("pdf", "沃克环流强度的年循环_allindex")

    resxy=True
    resxy@gsnDraw=False
    resxy@gsnFrame=False

    resxy@xyLineThicknesses= (/ 2.0/)          ; make second line thicker
    resxy@xyLineColors      = (/"red","pink","orange","green","blue","Purple"/)          ; change line color
    resxy@xyMarkerColors = (/"red","pink","orange","green","blue","Purple"/) 
    resxy@xyMonoMarkLineMode = True
    resxy@xyMarkLineMode  = (/"MarkLines"/) 
    resxy@xyMarker         =  (/16/)                    ; choose type of marker  
    resxy@xyMarkerSizeF     = 0.01        ; Marker size (default 0.01)
    resxy@xyMonoDashPattern = True
    resxy@xyDashPattern = 0 
  
    
    resxy@gsnLeftStringFontHeightF=0.01
    resxy@gsnRightStringFontHeightF=0.01

    resxy@trXMaxF = 12
    resxy@trXMinF = 1
    
    resxy@tmXBMode      = "Explicit"                ; Define own tick mark labels.
    resxy@tmXBValues=ispan(1,12,1)
    resxy@tmXBLabels = ""+x
   
    
    resxy@tmXBLabelFontAspectF=1.0
    resxy@tmXBLabelFontHeightF=0.008
    resxy@tmXBLabelFont=30     ;黑体
    resxy@tmXBLabelFontThicknessF=2.0
    resxy@tmXBMajorLengthF=0.005        ;調整刻度的长度
    resxy@tmXBMajorOutwardLengthF=0.005  ;调整刻度的长度，注意2个命令要同时使用
  

    resxy@tmXBMinorLengthF=0.003
    resxy@tmXBMinorOutwardLengthF=0.003

    resxy@tmYLMajorLengthF=0.005        ;調整刻度的长度
    resxy@tmYLMajorOutwardLengthF=0.005  ;调整刻度的长度，注意2个命令要同时使用
    resxy@tmYLLabelFont=30
    resxy@tmYLLabelFontHeightF=0.008
    resxy@tmYLMinorOn=False
    resxy@tmXTOn=False
    
    

    resxy@pmLegendDisplayMode    = "Always"            ; turn on legend
    resxy@pmLegendSide           = "Top"               ; Change location of 
    resxy@pmLegendParallelPosF   = .3                  ; move units right
    resxy@pmLegendOrthogonalPosF = -0.24                ; move units down
    resxy@pmLegendWidthF         = 0.15                ; Change width and
    resxy@pmLegendHeightF        = 0.1                ; height of legend.
    resxy@lgPerimOn              = False               ; turn off box around
    resxy@lgLabelFontHeightF     = .01                ; label font height
    resxy@lgItemOrder =(/5,4,3,2,1,0/)
    resxy@xyExplicitLegendLabels = (/"MSF","Ud","domega","U10","DSLP","SOI"/)         ; create explicit labels


    resxy@vpHeightF=0.4
    resxy@vpWidthF=0.6
    ; resxy@vpXF=0.15
    ; resxy@vpYF=0.95
    
    ; res2= resxy
    ; res2@xyLineColor      = (/"black"/)          ; change line color
    ; res2@xyMarkerColor = "black"
    ; res2@xyExplicitLegendLabels = (/"WL"/)         ; create explicit labels
    
    ; plot = gsn_csm_xy2(wks, ispan(1,12, 1), nino3, locationxi_month, resxy, res2)
     
     plot = gsn_csm_xy(wks, ispan(1,12, 1), (/msf_inten_month_std,ud_inten_month_std,w_inten_month_std,\
                                             u10_inten_month_std,slp_inten_month_std, soi_inten_month_std/), resxy)

     ; wks2=gsn_open_wks("pdf", "tahit年循环")
     ; resxy@tmYLValues = ispan(10040,10160,10)
     ; plot2 = gsn_csm_xy(wks2, ispan(1,12, 1), (/soi_inten_month_std/), resxy)
       
     ; draw(plot2)
     ; frame(wks2)


     draw(plot)
     frame(wks)
   end