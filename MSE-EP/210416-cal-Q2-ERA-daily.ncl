;*************************************************************
;Aim: cal daily q2 from ERA Datasets, Output 
; daily q2, file by year
;Mission:
;   1. Read data from ERA data
;   2. Post-post process
;     2.1 Calculate q2
;   3. Output to yearly files
;                      shuheng lin
;                   

;*************************************************************
begin


;======================================================
;Variables List:
;   t--------------Air Temperature (K, TLLL)
;   u--------------U Wind (m/s, TLLL)
;   v--------------V Wind (m/s, TLLL)
;   w--------------Vertical Velocity (Pa/s, TLLL)
;   z--------------Geopotential Height (m^2/s^2, TLLL)
;   slp------------Mean Surface Pressure (Pa, TLLL)
;======================================================
;global constants
    G_CPD   =   1004.07     ;Cpd=1004.07 J*kg^-1*K^-1
    G_RD    =   287.05      ;Rd=287.05 J*kg^-1*K^-1
    G_G     =   9.8         ;g=9.8 m*s^-2
;global variables
    ;Define parameters

    ;for directory 
    g_strt_yr=1979
    g_lst_yr=2016
    g_data_dir   =(/"shum","uwnd","vwnd","omega"/)
    g_data_path = "/home/yangsong3/data-observation/ERA-interim/pressure/"+g_data_dir+"/"+g_data_dir

    g_data_out_path = "/home/yangsong3/data-observation/ERA-interim/pressure/q2_1979-2016/"



    g_latS=-10.
    g_latN=20.
    g_lonW=100.
    g_lonE=170.
    ; create an array of desired pressure levels:

    ;g_fig_title   =   "1979-2014 "+(/"Mar","Apr","May","Jun","Jul","Aug"/)+" q2 Trend"

;*************************For F_Compset Ctrl***********************************

; ==============================================================
; Mission 1. Open the file: Read only the user specified period
; ==============================================================
    do nyr = g_strt_yr, g_lst_yr
        ctrl_in1     = addfile (g_data_path(0)+".interim."+nyr+".nc", "r")   ; note the "s" of addfile
        q  =ctrl_in1->shum(:,:,{g_latS:g_latN},{g_lonW:g_lonE})
        ctrl_in2     = addfile (g_data_path(1)+".interim."+nyr+".nc", "r")   ; note the "s" of addfile
        u  = ctrl_in2->uwnd(:,:,{g_latS:g_latN},{g_lonW:g_lonE})
        ctrl_in3     = addfile (g_data_path(2)+".interim."+nyr+".nc", "r")   ; note the "s" of addfile
        v  =ctrl_in3->vwnd(:,:,{g_latS:g_latN},{g_lonW:g_lonE})
        ctrl_in4     = addfile (g_data_path(3)+".interim."+nyr+".nc", "r")   ; note the "s" of addfile
        w  = ctrl_in4->omega(:,:,{g_latS:g_latN},{g_lonW:g_lonE})

        printMinMax(q, 1)
        printMinMax(u, 1)
        printMinMax(v, 1)
        printMinMax(w, 1)

        print(nyr+" Read done!")
    ;   Prepare some global variables

        g_lvl   =   q&level*100 
        g_lat   =   q&lat
        g_lon   =   q&lon
        g_nlat  =   dimsizes(g_lat)
        g_nlon  =   dimsizes(g_lon)
        g_nlvl =   dimsizes(g_lvl)
    ; ==============================================================
    ; Mission 2. post-post process data to deserved results
    ; ==============================================================

    ;-----------------------------------------------------------
    ;q2 Caculation: (partial diff replaced by complete diff)
    ;    
    ;
    ;   q2  =  -L*( d(q)/dt + V*grad(q) + w*dq/dp)
    ;               [1]         [2]         [3]
    ;
    
    ;-----------------------------------------------------------
    
    ;CALCULATE Q2         
    ;*****************TERM [1], easy********************
    dqdt    =   q
    dqdt    =   center_finite_diff_n(q,86400,False,0,0)
    ;***************************************************


    ;***********TERM [2], this is a little tricky*******
    long_name = "moisture advection"
    units     = "kg/(kg-s)"     ; (m/s)*(g/kg)*(1/m)
    cyclic = False
    opt_adv = 0
    advq      = advect_variable_cfd(u,v,q,g_lat,g_lon, cyclic , long_name,units,opt_adv) 
   
    ;*****************TERM [3], easy********************
   
    dqdp = center_finite_diff_n(q,g_lvl,False,1,1)
    copy_VarCoords(q, dqdp)
    dqdp@longname = "Specific Humidity: Local Vertical Derivative"
    dqdp@units = "kg/(kg-Pa)"   ; (g/kg)/Pa 
    printVarSummary(dqdp)
    printMinMax(dqdp,0)
    print("-----")

    dqdp  = w*dqdp               ; overwrite ... convenience
    dqdp@longname = "Specific Humidity: Vertical Moisture Transport"
    dqdp@units    = "kg/(kg-s)"       ; [(Pa/s)(g/kg)/Pa)]  
    printVarSummary(dqdp)
    printMinMax(dqdp,0)
    print("-----")

    ;***************************************************
        Lc           = 2.26e6    ; [J/kg]=[m2/s2]  Latent Heat of Condensation/Vaporization
        Lc@long_name = "Latent Heat of Condensation/Vaporization"
        Lc@units     = "J/kg"  ; ==> "m2/s2"


        q2      =   dqdp ;copy metadata
        q2      =   dqdt + advq + dqdp
        q2      =   q2*Lc

        q2@long_name    =   "Daily Q2 on Pressure Levels"
        q2@units    =   "m2/s or J/kg/s"
        q2@var_desc    =   "Q2"
        q2!1        =   "lev"
        print(nyr+" Calculate done!")
    ;************************************************
    ; High Efficiency to Write the file
    ;************************************************

    ;Get dimsize
    ;    dims=dimsizes(q2(0,0,:,:))
    ;    nlat=dims(0)
    ;    nlon=dims(1)
        
    ;Set fileoption    
        system("rm "+g_data_out_path+"q2."+nyr+".nc")
        fout = addfile( g_data_out_path+"q2."+nyr+".nc","c")  ; open output netCDF file
        setfileoption(fout,"DefineMode",True)

    ;Set All field attribution
        fileAtt     = True
        fileAtt@creation_date=systemfunc("date")
        fileattdef(fout,fileAtt)

    ;Define Coordinate
        dimNames    =(/"time","lev","lat","lon"/)
        dimSizes    =(/-1,g_nlvl,g_nlat,g_nlon/)
        dimUnlim    =(/True,False,False,False/)
        filedimdef(fout,dimNames,dimSizes,dimUnlim)

    ;Define var, type and dim
        filevardef(fout,"time",typeof(q2&time),getvardims(q2&time))
        filevardef(fout,"lev",typeof(q2&lev),getvardims(q2&lev))
        filevardef(fout,"lat",typeof(q2&lat),getvardims(q2&lat))
        filevardef(fout,"lon",typeof(q2&lon),getvardims(q2&lon))
        filevardef(fout,"q2",typeof(q2),getvardims(q2))

    ;Define Attribute
        filevarattdef(fout,"q2",q2)
        filevarattdef(fout,"time",q2&time)
        filevarattdef(fout,"lev",q2&lev)
        filevarattdef(fout,"lat",q2&lat)
        filevarattdef(fout,"lon",q2&lon)
        
        fout->time=(/q2&time/)
        fout->lev=(/q2&lev/)
        fout->lat=(/q2&lat/)
        fout->lon=(/q2&lon/)
        fout->q2=(/q2/)

        print(nyr+" Write done!")
        delete([/q2, dqdp, dqdt, advq/])
        delete([/q, u, v, w/])
        print("Variable Deleted!")
    end do
end

