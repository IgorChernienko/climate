#Загружаем необходимые пакеты
require(ncdf4)
require(raster)
require(reshape)
require(maptools)
mask<-readShapePoly('masks/Sakhalin_east.shp') #Читаем шейп-файл с полигонами
###########################################
year<-c(1982,2016) #Задаем диапазон лет
month<-8           #Задаем месяц 
day<-7             #Задаем день 
prm<-'temp'#'ise'#  Выбираем тип данных -- ТПМ или лёд
###########################################
result<-matrix(nrow=year[2]-year[1]+1,ncol=nrow(mask)+1) #формируем матрицу
rownames(result)<-year[1]:year[2]                        #в которую поместим результаты
colnames(result)<-c(paste0('polygon_',1:nrow(mask)),'total_average')

# Запускаем цикл чтения netCDF по годам
for(ye in year[1]:year[2]){
  #формируем дату и диапазон координат для запроса
  tm.dia<-paste(paste(ye,month,sep='-'),day,sep='-')
  lat.dia<-as.vector(bbox(mask)[2,])#
  lon.dia<-as.vector(bbox(mask)[1,])
  vari<-ifelse(prm=='temp','sst','icec')
  fil.name<-paste0(paste0(paste0(vari,'.day.mean.'),ye),'.v2.nc')
  fil.name<-ifelse(prm=='temp',paste0('SST/',fil.name),paste0('ISE/',fil.name))
  #Открываем файл
  fid<-nc_open(fil.name)
  #Извлекаем значения точек измерений
  time.v<-as.vector(fid$dim$time$vals)
  lat.v<-as.vector(fid$dim$lat$vals)
  lon.v<-as.vector(fid$dim$lon$vals)
  #Приводим диапазоны извлечения данных к требованиям функции
  #по времени
  time.range<-as.Date(time.v,origin = '1800-01-01')==as.Date(tm.dia)
  tm.select<-c(
    length(time.v[as.Date(time.v,origin = '1800-01-01')<=as.Date(tm.dia,origin = '1900-01-01')]),
    1
  )
  #по широте
  lat.range<-lat.v>=lat.dia[1]&lat.v<=lat.dia[2]
  lat.select<-c(
    length(lat.v[lat.v<=min(lat.v[lat.range])]),
    length(lat.v[lat.range])
  )
  #по долготе
  lon.range<-lon.v>=lon.dia[1]&lon.v<=lon.dia[2]
  lon.select<-c(
    length(lon.v[lon.v<=min(lon.v[lon.range])]),
    length(lon.v[lon.range])
  )
  #Получаем массив со значениями температуры,[Долгота,Широта,Дата]
  q.res<-ncvar_get(nc = fid,varid = vari,
                   start = c(lon.select[1],lat.select[1],tm.select[1]),
                   count = c(lon.select[2],lat.select[2],tm.select[2]))
  rm(fid)
  dimnames(q.res)<-list(lon.v[lon.range],
                        lat.v[lat.range])
  sect.melt<-melt(q.res)
  rm(q.res)
  colnames(sect.melt)<-c('lon','lat','var')
  sect.sp<-sect.melt
  rm(sect.melt)
  coordinates(sect.sp)<-~lon+lat
  #Запускаем цикл по полигонам
  for(pg in 1:nrow(mask)){
    mask.sp<-mask[pg,]
    sect.over<-sect.sp%over%mask.sp
    sect.over<-as.vector(as.matrix(sect.over))
    sect.over<-sect.sp[!is.na(sect.over),]
    result[as.numeric(rownames(result))==ye,pg]<-mean(sect.over$var,na.rm=T)
    rm(sect.over)
  }
  #Вычисляем общее среднее для всех полигонов
  sect.over<-sect.sp%over%mask
  sect.over<-as.vector(as.matrix(sect.over))
  sect.over<-sect.sp[!is.na(sect.over),]
  result[as.numeric(rownames(result))==ye,nrow(mask)+1]<-mean(sect.over$var,na.rm=T) 
  rm(sect.over)
}
#Просмотр результата
fix(result)
#Графики временных рядов
plot(result[,1]~c(year[1]:year[2]),type='l',lwd=2,ylim=range(result,na.rm = T))
for(i in 2:nrow(mask)){
  lines(result[,i]~c(year[1]:year[2]),col=i,lwd=2)
}
lines(result[,nrow(mask)+1]~c(year[1]:year[2]),col='darkgray',lty=2,lwd=3)
grid(lty=3)

#Запись на диск графика временных рядов
win.metafile(filename = 'ncdf_time_ser.wmf')

plot(result[,1]~c(year[1]:year[2]),type='l',lwd=2,ylim=range(result,na.rm = T),
     xlab='Год',ylab='ТПМ')
for(i in 2:nrow(mask)){
  lines(result[,i]~c(year[1]:year[2]),col=i,lwd=2)
}
lines(result[,nrow(mask)+1]~c(year[1]:year[2]),col='darkgray',lty=2,lwd=3)
grid(lty=3)
dev.off()

png(filename = 'ncdf_time_ser.png',width = 800,height = 400,pointsize = 18)
plot(result[,1]~c(year[1]:year[2]),type='l',lwd=2,ylim=range(result,na.rm = T),
     xlab='Год',ylab='ТПМ')
for(i in 2:nrow(mask)){
  lines(result[,i]~c(year[1]:year[2]),col=i,lwd=2)
}
lines(result[,nrow(mask)+1]~c(year[1]:year[2]),col='darkgray',lty=2,lwd=3)
grid(lty=3)
dev.off()

#Запись на диск таблицы с результатами
write.csv2(x = result,file = 'ncdf_time_ser.csv')
