begin
    ;;read  eva data 
    year  = ispan(1958,2018,1) + ""

    fils =  "/home/yangsong3/data-observation/OAflux/monthly/evaporation/evapr_oaflux_"+year+".nc"

    f1 = addfiles(fils, "r")
    ; date = cd_calendar(f1[:]->time, 1)
    ; timeind = ind(date.ge.197901.and.date.le.201312)
    
    evp_r1 = short2flt(f1[:]->evapr(:,:,:))
     
    evp_r1@missing_value = 1.0e+36
    evp_r1@_FillValue =   1.0e+36
    evp_r1 = where(evp_r1.eq.3276.6, 1.0e+36, evp_r1)
    

    year  := ispan(2019,2020,1) + ""
    fils :=  "/home/yangsong3/data-observation/OAflux/monthly/evaporation/evapr_oaflux_"+year+".nc"
    f1 := addfiles(fils, "r")
    evp_r2 = dble2flt(f1[:]->evapr(:,:,:))
    
    evp_r2&time := toint(evp_r2&time) 
    evp  := array_append_record(evp_r1, evp_r2, 0)

    ; evp =  evp_r*10 / 365.0   ;;convert t0 day
    ; copy_VarMeta(evp_r, evp)
    ; delete(evp_r)
    printVarSummary(evp)
 
    ;;;处理数据
    ;;; 水平插值 lat:-88.75：88.75 ； lon:1.25 2.5×2.5度

    yo = fspan(-88.75, 88.75, 72)
    xo = fspan(1.25, 358.75, 144)
    evp_new = linint2_Wrap(evp&lon, evp&lat, evp, True, xo, yo, 0)

    delete(evp_new&LAT)
    delete(evp_new&LON)
    evp_new!1 = "lat"
    evp_new!2 = "lon"
    evp_new&lat = yo
    evp_new&lon = xo

   
   f2=addfile("/home/yangsong3/data-observation/linshh/data/sst/HadISST_sst.nc", "r")
   time = f2->time
   date:=cd_calendar(f2->time, 1)
   timeind:=ind(date.ge.195801.and.date.le.202012)

   time_sel = time(timeind)
   printVarSummary(time_sel)

   evp_new&time := time_sel



     f3 = addfile("evapr_oaflux_chazhi.nc", "c")
     f3->evapr = evp_new
     end 
    ;;;;


      
