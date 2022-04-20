begin


    fils1=systemfunc("ls /mnt/i/CMIP5/history/ua/CSIRO-Mk3-6-0/r1i1p1/ua_Amon_CSIRO-Mk3-6-0_historical_r1i1p1_199001-200512.nc")
    f1=addfile(fils1, "r")    
    delete(fils1)
    yearstart=199001
    yearend=200512

    time=f1->time    
    printVarSummary(time)
    date=cd_calendar(time, 1)
  
    print("******************日期*********************")
    timeind=ind(date.le.yearend.and.date.ge.yearstart)
    print(date(timeind))
    delete(date)
    print("读入风场")
    z0=(/925.0/)*100.0
    levind=get1Dindex(f1->plev,z0)

    uwnd_r=f1->ua(timeind,levind,:,:)  ;;;195801-201212
    print(uwnd_r)
    exit()
    end 