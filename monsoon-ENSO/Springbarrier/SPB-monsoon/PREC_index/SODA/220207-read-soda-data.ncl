

begin   
f1 = addfile("http://iridl.ldeo.columbia.edu/SOURCES/.CARTON-GIESE/.SODA/.v2p2p4/w/dods","r")
time = f1->time 
date = cd_calendar(time, 1)

temp = f1->w(1,:,{-5:5},{100:120})

print(temp&depth)
exit()
end