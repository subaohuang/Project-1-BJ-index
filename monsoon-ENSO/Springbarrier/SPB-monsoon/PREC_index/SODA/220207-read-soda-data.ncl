

begin   
f1 = addfile("http://iridl.ldeo.columbia.edu/SOURCES/.CARTON-GIESE/.SODA/.v2p2p4/temp/dods","r")
time = f1->time 
date = cd_calendar(time, 1)

temp = f1->temp(:,:,{-5:5},{100:120})

print(temp)
exit()
end