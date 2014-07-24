# Last Change: Mon Apr 21 10:45:23 AM 2014 CDT

#options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
#print(args)
# trailingOnly=TRUE means that only your arguments are returned, check:
# print(commandsArgs(trailingOnly=FALSE))

file <- args[1]
output <- args[2]

data <- read.table(file,sep="\t",colClasses=c("character","numeric","character","character","numeric","character","numeric","character","numeric","character","numeric","character","numeric"))

for (i in 1:nrow(data)) {
	counter <- 0; alpha <- c(NA,NA); num <- c(0,0)
#	print(data[i,])
	for (j in seq(5,11,2)) {
		if (data[i,j] > 5) {
		counter <- counter + 1
		alpha[counter] <- data[i,j-1]
		num[counter] <- data[i,j]
#		print(counter); print(i); print(j); print(alpha); print(num)
		}
	}
	if (num[1] >= 12 & num[2] == 0) {
#		print(data[i,])
#		cat("single", "alpha = ", alpha, "num = ", num, "\n\n")
		data[i,14] <- alpha[1]
	}
#	if (counter == 2 & (num[1]>=10 || num[2]>=10)) {
	if (counter == 2 & (num[1] + num[2]) >= 12) {
		total <- num[1] + num[2]; half <- round((num[1] + num[2])/2)
#		cat("total = ", total, "total/2 = ", half, "\n")
#		print(data[i,])
		p1 <- fisher.test(rbind(c(num[1],num[2]),c(total,0)), alternative="two.sided", conf.level=0.95)$p.value
#		cat("double1", "p1 (", total, ": 0 ) = ", p1, "\n")
		p2 <- fisher.test(rbind(c(num[1],num[2]),c(0,total)), alternative="two.sided", conf.level=0.95)$p.value
#		cat("double2", "p2 ( 0 :", total, ") = ", p2, "\n")
		p3 <- fisher.test(rbind(c(num[1],num[2]),c(half,half)), alternative="two.sided", conf.level=0.95)$p.value
#		cat("double3", "p3 (", half, ":", half, ") = ", p3, "\n")
#		if (p1 >= 0.05 & p1 >= p2 & p1 >= p3) {
		if (p1 >= p2 & p1 >= p3) {
#			cat("selected double1", "p1 (", total, ": 0 ) = ", p1, "\n\n")
			data[i,14] <- alpha[1]
		}
#		if (p2 >= 0.05 & p2 >= p1 & p2 >= p3) {
		if (p2 >= p1 & p2 >= p3) {
#			cat("selected double2", "p2 ( 0 :", total, ") = ", p2, "\n\n")
			data[i,14] <- alpha[2]
		}
#		if (p3 >= 0.05 & p3 >= p1 & p3 >= p2) {
		if (p3 >= p1 & p3 >= p2) {
#			cat("selected double3", "p3 (", half, ":", half, ") = ", p3, "\n\n")
			data[i,14] <- paste(alpha[1],alpha[2],sep="/")
		}
	}
#	if (counter > 2) { print("more than 2") }
}

#print(data); cat("\n")

write.table(data,file=output,row.names=F,col.names=F,quote=F,sep="\t")



