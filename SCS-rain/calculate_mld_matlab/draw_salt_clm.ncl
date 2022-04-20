

;;**************************************************************
;;;main program
;;***************************************************************
begin


;;;读入HadiSST海温数据
    diri="/mnt/f/share/SODA/soda3.4.2_salt_mn_1980-2015_chazhi_new.nc"
    
    f1=addfile(diri, "r")
    salt_r = f1->salt
    lat =f1->lat
    lon =f1->lon
    time = f1->time
    

    date = cd_calendar(time, 1)
    salt = tofloat(salt_r)
    copy_VarMeta(salt_r, salt)
    printVarSummary(salt)
  


    ; salt  = where(salt.eq.salt@_FillValue, 10^9, salt)
    ; salt@_FillValue = 10^9
    ; salt@missing_value = 10^9
    ; salt = where(salt.eq.0, salt@_FillValue,salt)
     

    ; fils2 = "/mnt/f/share/landsea.nc"
    ; f2 = addfile(fils2, "r")
    ; lsdata = f2->LSMASK
    ; lsm = landsea_mask(lsdata,salt&lat,salt&lon) 
       
    ; salt_mask = mask(salt, lsm.eq.1, False)
    ; copy_VarMeta(salt, salt_mask)
    ; printVarSummary(salt_mask)

      ; ;;;;写nc文件,简单的方法  
      ;      system("/bin/rm -f soda3.4.2_salt_mn_1980-2015_chazhi_new.nc")   ; remove any pre-existing file
            
      ;     ncdf = addfile("soda3.4.2_salt_mn_1980-2015_chazhi_new.nc" ,"c")  ; open output netCDF file

      ;     ;===================================================================
      ;     ; create global attributes of the file (optional)
      ;     ;===================================================================
      ;        fAtt               = True            ; assign file attributes
      ;        fAtt@title         = "SODA3.4.2_soda at 5m,198001-201312"
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
      ;        ncdf->salt  = salt                       ; 4D      
        
      ;     ;;;;
      ;     exit()

    salt_mask = salt
    salt_mask&lat@units = "degrees_north"
    salt_mask&lon@units = "degrees_east"


    density_clm = dim_avg_n_Wrap(salt_mask, 0)
    
    ; print(salt(1,:,:))
    printVarSummary(density_clm)
   
    ; printMinMax(salt(400,:,:),1)
 
    


  wks = gsn_open_wks("pdf", "/mnt/f/share/SODA/盐度的气候态")
  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnAddCyclic = False
  
  ; res@cnFillPalette="MPL_RdBu"
  res@cnFillOn=True
  res@cnLinesOn=False
  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   30.0
  res@cnMaxLevelValF       =   37.0
  res@cnLevelSpacingF      =  1.0
  res@cnLineLabelsOn=False
  res@lbLabelBarOn=True
  res@lbLabelStride=2
  res@lbOrientation="Vertical"
  res@mpLimitMode           = "LatLon"
  res@mpMinLatF             = min(lat)
  res@mpMaxLatF             = max(lat)
  res@mpMinLonF             = min(lon)
  res@mpMaxLonF             = max(lon)
  
  res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

  res@tmXTOn=False
  res@tmYROn=False
  res@tmXBOn=True
  res@tmYROn=False

  res@tmXBLabelFontHeightF=0.015
  res@tmYLLabelFontHeightF=0.015
  res@tmXBMajorThicknessF=3.0
  res@tmYLMajorThicknessF=3.0
  res@tmBorderThicknessF=4.0

  res@mpFillOn = False

  plot =gsn_csm_contour(wks,density_clm,res)
  draw(plot)
  frame(wks)


end 

;;;