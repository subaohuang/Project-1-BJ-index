  ;;********************************************************
;;计算海洋平流项的时间变化的比较，分为动力项和热力项的作用,只关注
;;南北两侧的水平平流的动力项与热力项以及垂直扩散项
;;                                        
;;                                     linshuheng 2020.03.11
;;
;;
;**********************************************************
function  rc_test(ssta_rc)
begin
  tval = onedtond(ssta_rc@tval , dimsizes(ssta_rc))
   df   = onedtond(ssta_rc@nptxy, dimsizes(ssta_rc)) - 2
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)    ; prob(nlat,nlon)
   copy_VarCoords(ssta_rc, prob)

return prob
end
function  wgt_annual(total_month)
begin
    weight1 = ispan(1,12,1)*1.0
    weight1 = (2*weight1-1)/24

    weight2 = ispan(1,12,1)*1.0
    weight2 = (25-2*weight2)/24
    total_annual_test = new(34,"float")
    total_annual_test@_FillValue = default_fillvalue("float")
    printVarSummary(total_annual_test)

    do i =0,32,1
      total_annual_test(i+1) = sum(total_month(i*12:(i+1)*12-1)*weight1) + sum(total_month((i+1)*12:(i+2)*12-1)*weight2)
    end do    
return total_annual_test
end    
    


;;**************************************************************
;;;main program
;;***************************************************************
begin
 ;;;time span
    startmon =198001
    endmon = 201512

    startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 4 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1
  ;;
;;;region
;;;region
   lat1 = -5
   lat2 = 20
   lon1 = 100
   lon2 = 170
;;;

   
;;;read soda data and then calculate ocean dynamic : bar(v*σT)， bar(w*σT)，bar denotes 
;;;averaged in mixed layer
    
    f2=addfile("/home/yangsong3/data-observation/SODA/3.4.2/soda3.4.2_u_mn_1980-2015_chazhi.nc", "r")   
    f3=addfile("/home/yangsong3/data-observation/SODA/3.4.2/soda3.4.2_v_mn_1980-2015_chazhi.nc", "r")   
    f4=addfile("/home/yangsong3/data-observation/SODA/3.4.2/soda3.4.2_wt_mn_1980-2015_chazhi.nc", "r")   
    f5=addfile("/home/yangsong3/data-observation/SODA/3.4.2/soda3.4.2_temp_mn_1980-2015_chazhi.nc", "r")
    f6=addfile("/home/ys17-19/lsh/Project/SCS-rain/extend_to_2015/SODA3.4.2_mld_density_criterion_0.8_2015.nc", "r") 
  


    indyear = (2015-1980+1)*12-1
    u_r = f2->u(0:indyear,{0:330},{lat1:lat2},{lon1:lon2})
    v_r = f3->v(0:indyear,{0:330},{lat1:lat2},{lon1:lon2})
    wt_r = f4->wt(0:indyear,{0:330},{lat1:lat2},{lon1:lon2})
    T_r= f5->temp(0:indyear,{0:330},{lat1:lat2},{lon1:lon2})
    mlt = f6->mld(0:indyear,{lat1:lat2},{lon1:lon2})
  


    u = tofloat(u_r)
    v = tofloat(v_r)
    wt = tofloat(wt_r)
    T = tofloat(T_r)

    copy_VarCoords(u_r, u)
    copy_VarCoords(v_r, v)
    copy_VarCoords(wt_r, wt)
    copy_VarCoords(T_r, T)


    delete(u_r)
    delete(v_r)
    delete(wt_r)
    delete(T_r)
  

    depth = T&depth
printVarSummary(wt)



;;;calculate lateral  and bottom boundaries of the domain Heat transport, bottom is 45m , lateral is 
  ;;; four horizatal boundaries

    ;;parameter 
     pi = get_pi("float")
     dx = u&lat
     dx = 6378388.*cos(0.0174533*u&lat)
     dx!0 = "lat"
     copy_VarCoords(u(1,1,:,1), dx)

     dy = (pi/180)*1*6378388.
   


      dz_3d_test = new(dimsizes(u),float)
      do i = 0,dimsizes(mlt&time)-1
        dz_3d_test(i,:,:,:) = dpres_plevel(depth, mlt(i,:,:), 0, 0)
      end do 
      copy_VarCoords(u, dz_3d_test)
      printVarSummary(dz_3d_test)
  
   
     dx_3d_test = conform_dims(dimsizes(T(:,:,{0:15},{125:160})) ,dx({0:15}) , 2)


     dx_3d_test  = where(ismissing(T(:,:,{0:15},{125:160})), 0.0, dx_3d_test)
     dz_3d_test  = where(ismissing(T), 0.0, dz_3d_test)
     
     volume = dim_sum_n(dx_3d_test*dz_3d_test(:,:,{0:15},{125:160})*dy, (/1,2,3/))


     Tbox =dim_sum_n_Wrap(T(:,:,{0:15},{125:160})*dy*dx_3d_test*dz_3d_test(:,:,{0:15},{125:160}), (/1,2,3/))/volume
     copy_VarCoords(T(:,1,1,1), Tbox)

     ; volume = sum(45*dy*45.58*dx({0:15})) 

     dT = center_finite_diff_n(Tbox, 1, False, 0, 0)
     copy_VarCoords(Tbox, dT)
     dT_annual  = month_to_annual(dT, 1) 
    
     Tbox_annual = month_to_annual(Tbox, 1)


  
     dT_test = center_finite_diff_n(Tbox_annual, 1, False, 0, 0)
     dT_wgt = dT_test
     dT_wgt(1:35) = Tbox_annual(1:35) - Tbox_annual(0:34)
   
     

    ;; lateral heat transport north and south v*(T -Tbox)dxdz,south: 0N ,110-155,north: 15N ,110-155
       x1 = 0
       x2 = 15
       y1  = 125.5
       y2 = 160
       z1  = 45.6
    
       dim =dimsizes(v(:,:,{x1},{y1:y2}))
    
       dz_3d1 = dz_3d_test(:,:,{x1},{y1:y2})
       dz_3d2 = dz_3d_test(:,:,{x2},{y1:y2})
       Tbox_3d = conform_dims(dim ,Tbox , 0)


    ;;;;;;;垂直扩散变化项
     
    ;;; lateral heat transport vertical w*(T -Tbox)dxdy 
       
       ; dim:= dimsizes(wt(:,4,{x1:x2},{y1:y2}))
       ; dx_3d :=conform_dims(dim ,dx({0:15}) , 1)
       ; Tbox_3d :=conform_dims(dim, Tbox, 0)       

       ; wtT = wt(:,4,{x1:x2},{y1:y2})
       ; mlt_scs = mlt(:,{x1:x2},{y1:y2})

       ; printVarSummary(mlt_scs)
      
       ; mlt_scs = where(ismissing(mlt_scs), avg(mlt_scs), mlt_scs)
      
       ; printMinMax(mlt_scs, 1)
       
;        ;;;将wt进行垂直方向差值，为了得到更加精确的结果
;            print("开始差值")
;             xi  = wt&depth 
;             xo = ispan(5,330,5)

;             wt_chazhi = linint1_n_Wrap(xi,wt(:,:,{x1:x2},{y1:y2}),False,xo,0,1)
;             T_chazhi = linint1_n_Wrap(xi,T(:,:,{x1:x2},{y1:y2}),False,xo,0,1)

;             wt_mlt = wt_chazhi(:,1,:,:)
;             T_mlt = T_chazhi(:,1,:,:)
;             T_mlt_diff = T_chazhi(:,1,:,:)

;             delete(wt)
;             delete(T)
        
;         print("start triple loops")
;          do i = 0,dim(0)-1
;             do j = 0,dim(1)-1
;                do k = 0,dim(2)-1
;                   wt_mlt(i,:,:) = wt_chazhi(i,{mlt_scs(i,j,k)},{x1:x2},{y1:y2})
;                   T_mlt(i,:,:) = T_chazhi(i,{mlt_scs(i,j,k)},{x1:x2},{y1:y2}) 
;                   T_mlt_diff(i,:,:) = T_chazhi(i,{mlt_scs(i,j,k)},{x1:x2},{y1:y2}) - Tbox(i)
;                end do 
;               end do
;          end do
;        ;;;
; ;;;;;;

; ;;;;写nc文件,简单的方法  
;      system("/bin/rm -f SODA3.4.2_wt_T_in_MLD.nc")   ; remove any pre-existing file
      
;     ncdf = addfile("SODA3.4.2_wt_T_in_MLD.nc" ,"c")  ; open output netCDF file

;     ;===================================================================
;     ; create global attributes of the file (optional)
;     ;===================================================================
;        fAtt               = True            ; assign file attributes
;        fAtt@title         = "SODA3.4.2 wt and T in mld"
;        fAtt@Conventions   = "None"
;        fAtt@creation_date = systemfunc ("date")
;        fileattdef( ncdf, fAtt )            ; copy file attributes

;     ;===================================================================
;     ; make time an UNLIMITED dimension; recommended  for most applications
;     ;===================================================================
;        filedimdef(ncdf,"time",-1,True) 

;     ;===================================================================
;     ; output variables directly; NCL will call appropriate functions
;     ; to write the meta data associated with each variable
;     ;===================================================================
;        ncdf->wt_mlt  = wt_mlt                          ; 4D               
;        ncdf->T_mlt = T_mlt  
;        ncdf->T_mlt_diff = T_mlt_diff                    ; 3D 
; exit()
;  ;;;;;;
  
  
   ;;δ(WT)
    f1 = addfile("/home/ys17-19/lsh/Project/SCS-rain/extend_to_2015/SODA3.4.2_wt_T_in_MLD_density_defined_0.8_2015.nc", "r")

    wt_mlt = f1->wt_mlt(:,{lat1:lat2},{lon1:lon2})
    T_mlt = f1->T_mlt_diff(:,{lat1:lat2},{lon1:lon2})

    wtT = wt_mlt
    wtT = wt_mlt*T_mlt
    


   wtT&lat@units = "degrees_north"
   wtT&lon@units = "degrees_east"


    dwt  = wtT(1,:,:)
    dwt = (/dim_avg_n_Wrap(wtT(12*18:,:,:), 0) - dim_avg_n_Wrap(wtT(0:12*18-1,:,:), 0)/)
   dwt&lat@units = "degrees_north"
   dwt&lon@units = "degrees_east"

    printVarSummary(dwt)

   ;;

   ;;δ(W)T
    tdw = wtT(1,:,:)
    T_clm  = dim_avg_n_Wrap(T_mlt(0:12*18-1,:,:), 0)
    w_delta = dim_avg_n_Wrap(wt_mlt(12*18:,:,:), 0) - dim_avg_n_Wrap(wt_mlt(0:12*18-1,:,:), 0)
    
    tdw = (/T_clm*w_delta/)
   
    ;;;

    ;;δ(T)w
    wdt = wtT(1,:,:)
    w_clm  = dim_avg_n_Wrap(wt_mlt(0:12*18-1,:,:), 0)
    T_delta = dim_avg_n_Wrap(T_mlt(12*18:,:,:), 0) - dim_avg_n_Wrap(T_mlt(0:12*18-1,:,:), 0)
    
    wdt = (/T_delta*w_clm/)


    ;;

   ;;δ(T)δ(w)
    dwdt = wtT(1,:,:)
    dwdt = (/w_delta*T_delta/)

   ;;

   ;;eddy 
   eddy = wtT(1,:,:)
   eddy = (/dwt - wdt - tdw - dwdt/)
   eddy&lat@units = "degrees_north"
   eddy&lon@units = "degrees_east"

   ;;
  
  dwt = (/dwt*10^5/)
  wdt = (/wdt*10^5/)
  tdw = (/tdw*10^5/)
  dwdt = (/dwdt*10^5/)
  eddy = (/eddy*10^5/)






;;;;draw


    
;;**************************************************************
;画图
;;***************************************************************
  wks   = gsn_open_wks("pdf", "WT的分解-philipian")

  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnLeftStringFontHeightF=20.0
  res@cnFillPalette="BlRe"
  res@cnFillOn=True
  res@cnLinesOn=False
  res@gsnAddCyclic = False

  res@cnLineLabelsOn=False
  res@lbLabelBarOn=False
  res@lbLabelStride=2
  res@lbOrientation="Vertical"
  res@mpLimitMode           = "LatLon"
  ; res@mpMinLatF             = -5
  ; res@mpMaxLatF             = 45
  ; res@mpMinLonF             =90
  ; res@mpMaxLonF             = 160

  res@mpMinLatF             = 0
  res@mpMaxLatF             = 15
  res@mpMinLonF             = 125
  res@mpMaxLonF             = 160



  res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

  res@tmXTOn=False
  res@tmYROn=False
  res@tmXBOn=True
  res@tmYROn=False
  res@tmXBMode="Explicit"
  res@tmYLMode="Explicit"

  res@tmXBLabelFontHeightF=0.015
  res@tmYLLabelFontHeightF=0.015
  res@tmXBMajorThicknessF=3.0
  res@tmYLMajorThicknessF=3.0
  res@tmBorderThicknessF=4.0

  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   -0.5
  res@cnMaxLevelValF       =   0.5
  res@cnLevelSpacingF      =  0.1
 

  plot = new(5,graphic)

  res@gsnLeftString = "Total"
  plot(0) = gsn_csm_contour_map(wks, dwt, res)
  res@gsnLeftString = "dW.Tclm"
  plot(1) = gsn_csm_contour_map(wks, tdw, res)
  res@gsnLeftString = "dT.Wclm"
  plot(2) = gsn_csm_contour_map(wks, wdt, res)
  res@gsnLeftString = "dT.dW"
  plot(3) = gsn_csm_contour_map(wks, dwdt, res)
  res@gsnLeftString = "EDDY"
  plot(4) = gsn_csm_contour_map(wks, eddy, res)
  



  pres = True
  pres@gsnPanelLabelBar  = True

  gsn_panel(wks, plot, (/3,2/), pres)
 


end







   ;;;

;;;;;





