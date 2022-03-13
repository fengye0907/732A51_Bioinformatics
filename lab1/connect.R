path <- "C:\\Users\\fengy\\Desktop\\asd.txt"
data <- read.table(file = path)
string <- as.character(data[1,1])
for (i in 2:95){
  string <- paste0(string,as.character(data[i,1]))
}
write.table(string, file = "C:\\Users\\fengy\\Desktop\\sd.txt", row.names = F, quote = F) 
