load "/mnt/d/ncl/Dot-Matrix-master/WCHAR.ncl" ; All related function tar into this lib 
;----------------------------------------------------------------------
; Procedure to attach a horizontal labelbar to the bottom of plot.
;----------------------------------------------------------------------
undef("add_labelbar")
procedure add_labelbar(wks,plot,colors,labels)
local vph, vpw, nboxes, lbres, lbid, amres, annoid
begin
  getvalues plot                      ; Get plot size for use in
    "vpHeightF" : vph                 ; creating labelbar.
    "vpWidthF"  : vpw
  end getvalues

  nboxes = dimsizes(colors)

  lbres                    = True          ; labelbar only resources
  lbres@lbAutoManage       = True          ; Necessary to control sizes
  lbres@vpWidthF           = 0.15 * vpw    ; labelbar width
  lbres@vpHeightF          = 0.95 * vph    ; labelbar height
  lbres@lbFillColors       = colors        ; labelbar colors
  lbres@lbMonoFillPattern  = True          ; Solid fill pattern
  lbres@lbLabelFontHeightF = 0.1          ; font height. default is small

  lbres@lbLabelAlignment   = "BoxCenters"  ; center of box
  ; lbres@lbLabelAlignment   = "BoxCenters"  ; center of box
  lbres@lbLabelAlignment   = "ExternalEdges"  ; bottom of box
  lbres@lbOrientation      = "Vertical"
  lbres@lbPerimOn          = False

  lbid = gsn_create_labelbar(wks,nboxes,labels,lbres)
   
;
; Now, create some annotation resources indicating how we want to
; attach the labelbar to the plot. Here, we are using the top right
; corner of the labelbar as the point which we are going to position
; it, and then we use amParallelPosF and amOrthogonalPosF to indicate
; where we want to place it.
;
; amParallelPosF/amOrthogonalPosF                                       ;
;    0.0/ 0.0  -  annotation in dead center of plot
;    0.5/ 0.5  -  annotation at bottom right of plot
;    0.5/-0.5  -  annotation at top right of plot
;   -0.5/-0.5  -  annotation at top left of plot
;   -0.5/ 0.5  -  annotation at bottom left of plot
;                                                                 
  amres                  = True
  amres@amJust           = "CenterLeft"
  amres@amParallelPosF   =  0.52
  amres@amOrthogonalPosF =  0.0
  plot@annoid = gsn_add_annotation(plot,lbid,amres)


end

undef("rm_enso")
function rm_enso(data,nino3_4_anu_std)
begin
  ;;;去掉ENSO信号
  rc_nino=regCoef_n(nino3_4_anu_std,data,0,0)
  data_fit = rc_nino*nino3_4_anu_std + rc_nino@yintercept
  data_resi = data - data_fit
  ;;;

return (data_resi)
end 

;;***********************************************************************
;;主程序
;;;**********************************************************************
begin
    diri="/mnt/f/share/rea_msf/"
    f1=addfile(diri+"intension_rea_r1i1p1_djf.nc", "r")
    date=f1->time
    timeind=ind(date.ge.1979.and.date.le.2012)
    
    inten_rea_pa=f1->intension_pa(:,timeind)
    inten_rea_id=f1->intension_id(:,timeind)
    inten_rea_alt=f1->intension_alt(:,timeind)

    inten_rea_id = (/-1*inten_rea_id/)
   
    delete(timeind)
    delete(date)
    delete(f1)
    delete(diri)

    
    diri="/mnt/f/share/rea_msf/"
    f1=addfile(diri+"intension_rea_r1i1p1_jja.nc", "r")
    date=f1->time
    timeind=ind(date.ge.1979.and.date.le.2012)
    
    inten_rea_pa_jja=f1->intension_pa(:,timeind)
    inten_rea_id_jja=f1->intension_id(:,timeind)
    inten_rea_alt_jja=f1->intension_alt(:,timeind)
    
    inten_rea_id_jja = (/-1*inten_rea_id_jja/)

    delete(timeind)
    delete(date)
    delete(f1)
    delete(diri)
;;;;

;;;;load data olr

   diri="/mnt/f/share/NOAA/"
   fils1=systemfunc("ls "+diri+"olr.mon.mean.nc")
    ; fils1=systemfunc("ls "+diri+"/MERRA/MERRA/msf_time_MERRA_global.nc")
    ; fils1=systemfunc("ls "+diri+"MRE/msf_time_mre_global.nc")
    year1=1979
    year2=2012
    yearend=201212
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201212)
    olr_r=f1->olr(timeind,:,:)
    olr=(/olr_r@add_offset+olr_r@scale_factor*olr_r/)
    copy_VarCoords(olr_r, olr)
    copy_VarAtts(olr_r, olr)
    delete(olr&time)
    olr&time=date(timeind)
    printVarSummary(olr)
    olr_pi=lonPivot(olr, 30)
   ;;;先得到秋季和冬季eninio年的沃克环流情况SONDJF
    olr_djf_r=month_to_season(olr_pi(time|:,lat|:,lon|:), "DJF")
    olr_mam_r=month_to_season(olr_pi(time|:,lat|:,lon|:), "MAM")

    olr_jja_r=month_to_season(olr_pi(time|:,lat|:,lon|:), "JJA")
    olr_son_r=month_to_season(olr_pi(time|:,lat|:,lon|:), "SON")
  
    ; olr_djf_r=(/month_to_annual(olr_pi(time|:,lat|:,lon|:), 1)/)
    olr_djf=olr_djf_r
    olr_djf=(/(olr_djf_r+olr_mam_r)/2.0/)
    olr_jja = olr_jja_r
    olr_jja = (/(olr_jja_r+olr_son_r)/2.0/)

    olr_djf_pool=dim_avg_n_Wrap(olr_djf(:,{10:10},{80:150}), (/1,2/))
    olr_jja_pool=dim_avg_n_Wrap(olr_jja(:,{10:10},{80:150}), (/1,2/))

    delete(date)
    delete(timeind)
    ;;;

    f2=addfile("/mnt/f/share/GPCP/precip.mon.mean.nc", "r")
    ; f2=addfile("/mnt/f/share/GPCC/precip.mon.total.v7.nc", "r")
    date=cd_calendar(f2->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201212)
    p=f2->precip(timeind,:,:)

    delete(p&time)
    p&time=date(timeind)
    printVarSummary(p)
    p_pi=lonPivot(p, 30)
   ;;;先得到秋季和冬季eninio年的沃克环流情况SONDJF
    p_mam_r=month_to_season(p_pi(time|:,lat|:,lon|:), "MAM")
    p_djf_r=month_to_season(p_pi(time|:,lat|:,lon|:), "DJF")
    
    p_jja_r = month_to_season(p_pi(time|:,lat|:,lon|:), "JJA")
    p_son_r = month_to_season(p_pi(time|:,lat|:,lon|:), "SON")
    
    p_jja = p_jja_r
    p_jja = (/(p_jja_r+p_son_r)/2.0/)

    p_djf = p_djf_r
    p_djf = (/(p_djf_r+p_mam_r)/2.0/)
    
    
    p_djf_pool=dim_avg_n_Wrap(p_djf(:,{10:10},{80:150}), (/1,2/))
    p_jja_pool=dim_avg_n_Wrap(p_jja(:,{10:10},{80:150}), (/1,2/))


    delete(date)
    delete(timeind)
    delete(p)
    delete(p_pi)
    ;;;
    
    ;;;load sst dataset
    ;;;读入nino3.4 inde
      f2=addfile("/mnt/f/share/ersst/HadISST_sst.nc", "r")
      date=cd_calendar(f2->time, 1)
      timeind=ind(date.ge.197901.and.date.le.201312)
      sst=lonFlip(f2->sst(timeind,{20:-20},:))
      ssta=(/rmMonAnnCycTLL(sst)/)
      copy_VarCoords(sst, ssta)
      ssta_de=dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)
      copy_VarCoords(ssta, ssta_de)
      nino3_4=dim_avg_n_Wrap(ssta_de(:,{5:-5},{190:240}), (/1,2/))   
      ; nino3_4=(/dtrend(nino3_4, True)/)
      ; nino3_4_anu=month_to_annual(nino3_4, 1)
      nino3_4_anu=(/month_to_season(nino3_4, "DJF")/)

      nino3_4_djf=nino3_4_anu(0:33)   ;;;注意要JJA的相关，Nino指数要用前期的事件，研究Pacific这个问题
      ; nino3_4_djf=nino3_4_anu(1:)
      nino3_4_anu_std=dim_standardize(nino3_4_djf, 1)
      nino3_4_anu_std!0="time" 
      nino3_4_anu_std&time=ispan(1979, 2012, 1)


    ;;;

      ;;;;强度与Enino指数的相关性

      r_id_jja =  escorc_n(inten_rea_id_jja, nino3_4_anu_std,1,0)
      r_pa_jja =  escorc_n(inten_rea_pa_jja, nino3_4_anu(1:34),1,0)

      r_id_djf =  escorc_n(inten_rea_id, nino3_4_anu_std,1,0)
      r_pa_djf =  escorc_n(inten_rea_pa, nino3_4_anu_std,1,0)

      print("r_id_jja="+r_id_jja)
      print("r_pa_jja="+r_pa_jja)
      print("r_id_djf="+r_id_djf)
      print("r_pa_djf="+r_pa_djf)

      ;;;;

;;;求相关性


  ; ;; 去掉做相关前的ENSO信号
  ;   do i = 0,6
  ;   inten_rea_id(i,:) =(/rm_enso(inten_rea_id(i,:),nino3_4_anu_std)/)
  ;   inten_rea_pa(i,:)  = (/rm_enso(inten_rea_pa(i,:),nino3_4_anu_std)/)
  ;   inten_rea_alt(i,:)  = (/rm_enso(inten_rea_alt(i,:),nino3_4_anu_std)/)
    
  ;   inten_rea_id_jja(i,:) = (/rm_enso(inten_rea_id_jja(i,:),nino3_4_anu_std)/)
  ;   inten_rea_pa_jja(i,:)  = (/rm_enso(inten_rea_pa_jja(i,:),nino3_4_anu_std)/)
  ;   inten_rea_pa_jja(i,:)  = (/rm_enso(inten_rea_pa_jja(i,:),nino3_4_anu(1:34))/)
  ;   ; inten_rea_alt_jja(i,:)  = (/rm_enso(inten_rea_alt_jja(i,:),nino3_4_anu_std)/)
  ;   end do
  ;   olr_djf_pool  = (/rm_enso(olr_djf_pool,nino3_4_anu_std)/)
  ;   olr_jja_pool  = (/rm_enso(olr_jja_pool,nino3_4_anu(1:34))/)
  ;   ; olr_jja_pool  = (/rm_enso(olr_jja_pool,nino3_4_anu_std)/)
    
  ;   p_djf_pool  = (/rm_enso(p_djf_pool,nino3_4_anu_std)/)
  ;   p_jja_pool  = (/rm_enso(p_jja_pool,nino3_4_anu(1:34))/)
  ;   p_jja_pool  = (/rm_enso(p_jja_pool,nino3_4_anu_std)/)
  ; ;; 
    

;;; 重新计算相关性 
   cor_jja = escorc_n(inten_rea_pa_jja, inten_rea_id_jja,1,1)
   cor_djf = escorc_n(inten_rea_id, inten_rea_pa,1,1)
   printVarSummary(cor_djf)

   ;;; OLR与强度的相关性
   
   cor_jja_pa_olr = escorc_n(inten_rea_pa_jja,olr_jja_pool ,1,0)
   cor_jja_id_olr = escorc_n(inten_rea_id_jja,olr_jja_pool ,1,0)
   
   cor_djf_pa_olr = escorc_n(inten_rea_pa,olr_djf_pool ,1,0)
   cor_djf_id_olr = escorc_n(inten_rea_id,olr_djf_pool ,1,0)
 
   ;;; precipitation correlation 

   cor_jja_pa_p = escorc_n(inten_rea_pa_jja,p_jja_pool ,1,0)
   cor_jja_id_p = escorc_n(inten_rea_id_jja,p_jja_pool ,1,0)
   
   cor_djf_pa_p = escorc_n(inten_rea_pa,p_djf_pool ,1,0)
   cor_djf_id_p = escorc_n(inten_rea_id,p_djf_pool ,1,0)
  
   ;;;显著性检验
   ;;;自由度
     
     df=34
     p  = 0.975
     T = cdft_t(p, df)    ; T =  2.12
     print("T=="+T)
    t=T
    r_test=sqrt(t^2/(34-2+t^2))                     ;0.05显著性检验的t值为2.101
    print(r_test)
;;;;






; Read the labels
; ncol_labels = (/"SI","WI","SPO","SIO","WPO","WIO","SPP","SIP","WPP","WIP"/)
ncol_labels = (/"WI","WPO","WIO","WPP","WIP","SI","SPO","SIO", "SPP","SIP"/)
nrow_labels =(/"20CR","ERAIM","JRA-55","MERRA","NCEP-1","NCEP-2","ENS"/)


values=new((/7,10/),typeof(cor_jja))

values(:,0)=cor_djf


values(:,1) = cor_djf_pa_olr  
values(:,2) = cor_djf_id_olr

values(:,3) = cor_djf_pa_p  
values(:,4) = cor_djf_id_p  

values(:,5)=cor_jja

values(:,6) = cor_jja_pa_olr  
values(:,7) = cor_jja_id_olr

values(:,8) = cor_jja_pa_p  
values(:,9) = cor_jja_id_p  


text_values =  sprintf("%3.2f", values) + ""

ncols=10
nrows=7
ncr    = dimsizes(values)
col_width  = 1./ncols
row_width  = 1./nrows
col_width2 = col_width/2.
row_width2 = row_width/2.


;;***********************************************************************
;;画table图
;;;**********************************************************************

; Define colors to use for the table. Be sure to include white and black for background
; and foreground.

; ;
  colors = (/"White","Black","Navy","MediumBlue","MediumTurquoise","SeaGreen1",\
             "azure1","azure1","Orange","Firebrick1","Red","Brown"/)

;---Set levels and colors to use
       
        ; rgba   = read_colormap_file("BkBlAqGrYeOrReViWh200")
        ; colors = span_color_rgba(rgba(25:185,:),10)
        ; colors(0,3) = 0.0   ; make the first color transparent

  labels = new(dimsizes(colors)-1,string)    ; for labelbar

; Assign a color to a particular range of data.
  cols_for_vals = new(ncr,"string")
 
  cols_for_txt = new(ncr,"string")

  rng_beg=(/-1.0,-0.8,-0.6,-0.4,-0.35,0.0,0.35,0.4,0.6,0.8/)
  rng_end=(/-0.8,-0.6,-0.4,-0.35,0.0,0.35,0.4,0.6,0.8,1.0/)

  label_beg=(/-1.0,-0.8,-0.6,-0.4,-0.35,0.0,0.35,0.4,0.6,0.8,1.0/)

  do i=0,9,1
    
   
    cols_for_vals = where(rng_beg(i).le.values.and.values.lt.rng_end(i),colors(i+2),\
                          cols_for_vals)
;---Create formatted labels for the labelbar
    ; labels(i) = sprintf("%3.2f", rng_beg(i)) + " - " + \
    ;               sprintf("%3.2f",rng_end(i))            
  end do
   labels= ""+sprintf("%3.2f",label_beg)

  cols_for_txt = where(cols_for_vals.eq."Navy".or.cols_for_vals.eq."MediumBlue","white","black")

; Start the graphics.

   wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/相关系数的table图")
  gsn_define_colormap(wks,colors)

;
; Create a blank plot so we can get some tickmarks.
;
; Do the blank plot before the table, so we can maximize
; the blank plot. Later we'll make sure the table goes
; in the same location as the blank plot.
;
  resb                         = True
  resb@gsnDraw                 = False
  resb@gsnFrame                = False

  resb@vpYF                    = 0.92
  resb@vpXF                    = 0.15
  resb@vpHeightF               = 0.8
  resb@vpWidthF                = 0.7

; Explicitly label X axis. The blank plot goes from 0 to 1, by default.
  resb@tmXBMode                = "Explicit"
  resb@tmXBValues              = fspan(col_width2,1.-col_width2,ncols)
  resb@tmXBLabels              = ncol_labels
  resb@tmXBLabelFontHeightF    = 0.015            ; Make font a little smaller

; Explicitly label Y axis.
  resb@tmYLMode                = "Explicit"
  resb@tmYLValues              = fspan(row_width2,1.-row_width2,nrows)
;
; The Y labels need to be reversed, because the label in the first
; row in the ASCII file is associated with the top of the plot.
;
  resb@tmYLLabels              = nrow_labels(::-1)
  resb@tmYLLabelFontHeightF    = 0.015            ; Make font a little smaller

; Rotate the X axis labels 90 degrees.
  resb@tmXBLabelAngleF         = 90.
  resb@tmXBLabelJust           = "CenterRight"

; Turn off right and top tickmarks
  resb@tmYROn                  = False
  resb@tmXTOn                  = False

; Main title
  resb@tiMainString            = ""

  blank = gsn_csm_blank_plot(wks,resb)
  add_labelbar(wks,blank,colors(2:),labels)   ; Attach labelbar
;
; Get position and size of the blank plot so we can
; be sure to draw the table in same location.
;
  getvalues blank
    "vpXF"      : vpx
    "vpYF"      : vpy
    "vpWidthF"  : vpw
    "vpHeightF" : vph
  end getvalues

; Start and end coordinates for table.
  x = (/vpx,vpx+vpw/)
  y = (/vpy-vph,vpy/)

; Set up resources for table
  rest               = True

  rest@gsLineColor   = -1                ; No border lines
  rest@gsFillColor   = cols_for_vals     ; Each cell of table
  rest@txFontColor   = cols_for_txt

  rest@txFontHeightF = 0.015             ; Size of centered text

  gsn_table(wks,ncr,x,y,text_values,rest)    ; Draw table
  draw(blank)                                ; Draw tickmarks and labelbar

; Advance the frame.
  frame(wks)

  






   end 






                     



    

