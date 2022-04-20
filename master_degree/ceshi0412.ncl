    begin
    xi_model="CMCC-CM"
    fils1=systemfunc("ls /mnt/j/CMIP5/history/uv/uv_"+xi_model+"_history_r1i1p1_"+"*"+".nc")
    
    yearstart=197901
    yearend=200512
    print(fils1)
    f1=addfile(fils1, "r")

    time=f1->time    
    printVarSummary(time)
    
    date=cd_calendar(time, 1)
    print("******************日期*********************")
    timeind=ind(date.le.yearend.and.date.ge.yearstart)
    print(date(timeind))
    end