close;
clear;
clc;

//Extracting Radiation Data
[fd,SST,Sheetnames,Sheetpos] = xls_open('Model Inputs.xls')
[Value,TextInd] = xls_read(fd,Sheetpos(5))
Rs = Value(2:8761,7)
mclose()

//Extracting Azimuth Angle
[fd,SST,Sheetnames,Sheetpos] = xls_open('Orientation.xls')
[Value,TextInd] = xls_read(fd,Sheetpos(2))
Az = Value(2:366,2:25)'
mclose()
Az = Az(:)

//For a north-facing wall - Azimuth between 90 & 270 deg
for i = 1:length(Az)
    if (Az(i) > 180 & Az(i) < 360) 
    Rs(i) = Rs(i)
    else 
    Rs(i) = 0
    end
end

//Accounting for various seasons and getting seasonal average
//Summer - December to end of February
Rsummer = [Rs(1:1416); Rs(8017:8760)]
//Autumn - March to end of May
Rautumn = Rs(1417:3624)
//Winter - June to end of August
Rwinter = Rs(3625:5832)
//Spring - September to end of November
Rspring = Rs(5833:8016)

//Number of 24 hr periods
np_Rsummer = length(Rsummer)/24
np_Rautumn = length(Rautumn)/24
np_Rwinter = length(Rwinter)/24
np_Rspring = length(Rspring)/24

j = 1
k = 1
for l = 1:length(Rsummer)
    Sum24(k,j) = Rsummer(l)
    k = k +1
    if k == 25
        k = 1
        j = j+1
    end 
end

for m = 1:24
    SumAve(m) = mean(Sum24(m,1:$))
end

j = 1
k = 1
for l = 1:length(Rautumn)
    Aut24(k,j) = Rautumn(l)
    k = k +1
    if k == 25
        k = 1
        j = j+1
    end 
end

for m = 1:24
    AutAve(m) = mean(Aut24(m,1:$))
end

j = 1
k = 1
for l = 1:length(Rwinter)
    Wint24(k,j) = Rwinter(l)
    k = k +1
    if k == 25
        k = 1
        j = j+1
    end 
end

for m = 1:24
    WintAve(m) = mean(Wint24(m,1:$))
end

j = 1
k = 1
for l = 1:length(Rspring)
    Spr24(k,j) = Rspring(l)
    k = k +1
    if k == 25
        k = 1
        j = j+1
    end 
end

for m = 1:24
    SprAve(m) = mean(Spr24(m,1:$))
end

Ave = [SumAve AutAve WintAve SprAve]

scf(1)
clf()
plot((0:23),Ave)
plot(0,1000)
title('West-facing')
xlabel('Time of day [hr]')
ylabel('Solar radiation [W/m^2]')
legend('Summer','Autumn','Winter','Spring')
