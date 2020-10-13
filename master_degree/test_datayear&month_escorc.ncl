function clm_ano(x[*]:numeric)  ;去除一维时间序列的年循环
begin

   dimx  = dimsizes (x)
   ntim  = dimx(0)

   nmos  = 12
   modCheck ("clm", ntim, nmos)   ; error check; Compute all 12 monthly averages first. 
   aveMonth = new((/nmos/),typeof(x) \
                                    ,getFillValue(x))

   do nmo=0,nmos-1
      aveMonth(nmo) = dim_avg( x(nmo:ntim-1:nmos) )
   end do
   copy_VarAtts (x,aveMonth)
   aveMonth!0   = "month"           ; create a "month" named dim
   aveMonth&month = ispan(0,nmos-1,1) ; create a month coord var
   xAnom=x
   do yr=0,ntim-1,nmos
      xAnom(yr:yr+nmos-1) = (/ x(yr:yr+nmos-1)- aveMonth /)
   end do
   return (xAnom)
end

begin
;;;读入指数数据
    diri="/mnt/f/share/era/data/"
    ipindex=asciiread("/mnt/f/share/era/data/ip_index_1950_2010.ascii", 732, "float")

   
    ipindex_mon_std=dim_standardize_n(clm_ano(ipindex),0,0)
    ipindex_mon_std!0="time"
    ipindex_mon_std&time=ispan(1,732,1)
    ipindex_mon_std_de=dtrend(ipindex_mon_std, True)
    copy_VarCoords(ipindex_mon_std, ipindex_mon_std_de)
     
    ipindex_mon_std_de=(/runave_n(ipindex_mon_std_de, 3, 0, 0)/)


    ipindex_annual=dim_standardize_n(month_to_annual(ipindex, 1), 0, 0)
    ipindex_annual_std_de=dtrend(ipindex_annual, True)
    
    ipindex_season_std_de=month_to_seasonN(ipindex_mon_std_de, "DJF")

;;;

    fils2=systemfunc("ls /mnt/f/share/ersst/HadISST_sst.nc") 
    f2=addfile(fils2, "r")
    date=cd_calendar(f2->time, 1)
    timeind=ind(date.le.201012.and.date.ge.195001)
    sst=lonFlip(f2->sst(timeind,{20:-20},:))
    delete(sst&time)

    sst&time=date(timeind)
    
    ssta=rmMonAnnCycTLL(sst)

    sst_de=dtrend_msg_n(ispan(0,dimsizes(ssta&time)-1,1), ssta, True, False, 0)
    copy_VarCoords(ssta, sst_de)
   

;;;读入nino指数数据
    nino3_4=asciiread(diri+"Hadley_nino4_1950_2010_nodtrend.ascii", 732, "float")
    nino3_4=(/dtrend(nino3_4, True)/)
    nino3_4_mon_std=dim_standardize(clm_ano(nino3_4), 1)
    nino3_4_mon_std!0="time" 
    nino3_4_mon_std&time=ispan(1, 732, 1)

    nino3_4_mon_std=(/runave_n(nino3_4_mon_std, 3, 0, 0)/)
;;;
    
   nino3_4_mon_std_season=month_to_seasonN(nino3_4_mon_std, "DJF")
   
   nino3_4_annual_std=month_to_annual(nino3_4_mon_std, 1)

   printVarSummary(nino3_4_mon_std_season)
  
 



;;;计算相关性

    r1=escorc(nino3_4_mon_std, ipindex_mon_std_de)
    print("r1="+r1)

    r2=escorc(nino3_4_mon_std_season, ipindex_annual_std_de)

    print("r2="+r2)
   

    r3=escorc(nino3_4_annual_std, ipindex_annual_std_de)

    print("r3="+r3)
    

    r4=escorc(nino3_4_mon_std_season, ipindex_season_std_de)

    print("r4="+r4)

end 


   

   













