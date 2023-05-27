##
food <- read.csv(file="06.csv", + head=FALSE,sep=",")
type <- food[,2:ncol(food)]
head(type)

foodwithnames <- read.csv(file="06.csv", head=TRUE,sep=",")
##two way to normalize
normalized1_M <- apply(type, 2, function(x) (x - mean(x)) / sd(x))
normalized2_M <- apply(type, 2, function(x) (x - min(x)) / (max(x) - min(x)))
##distance inside each way
dis1 = dist(normalized1_M , method="minkowski",p=2)
dis2 = dist(normalized2_M , method="minkowski",p=2)
dis3 = dist(type , method="minkowski",p=2)
##MDS 
model1_1 <- cmdscale(dis1,k=1)
model1_2 <- cmdscale(dis1,k=2)
model1_3 <- cmdscale(dis1,k=3)
model1_4 <- cmdscale(dis1,k=4)
model2_1 <- cmdscale(dis2,k=1)
model2_2 <- cmdscale(dis2,k=2)
model2_3 <- cmdscale(dis2,k=3)
model2_4 <- cmdscale(dis2,k=4)
model3_1 <- cmdscale(dis3,k=1)
model3_2 <- cmdscale(dis3,k=2)
model3_3 <- cmdscale(dis3,k=3)
model3_4 <- cmdscale(dis3,k=4)
##MDS model give us similarity between points (only for 1D and 2D)
par(mfrow = c(1, 2))
plot(model1_1, rep(0, length(model1_1)),asp=1, ylim=c(-2,2),xlim=c(-3,4), ylab="",xlab="")
plot(model1_2,asp=1, ylim=c(-2,2),xlim=c(-3,4), ylab="",xlab="")
par(mfrow = c(1, 2))
plot(model2_1, rep(0, length(model2_1)),asp=1, ylim=c(-2,2),xlim=c(-1.5,1.5), ylab="",xlab="")
plot(model2_2,asp=1, ylim=c(-2,2),xlim=c(-1.5,1.5), ylab="",xlab="")
par(mfrow = c(1, 2))
plot(model3_1, rep(0, length(model3_1)),asp=1, ylim=c(-200,200),xlim=c(-250,200), ylab="",xlab="")
plot(model3_2,asp=1, ylim=c(-200,200),xlim=c(-250,200), ylab="",xlab="")

library(wordcloud)
textplot(model2_2[,1],model2_2[,2], gsub("(ˆ[ˆ\\s]+\\s{1})","",food$V1,perl=TRUE),
         asp=1,xlim=c(-1,1),ylim=c(-1,1),
         cex=0.6)



##what dimention is enough
par(mfrow = c(1, 3))
eigen_values1 <- cmdscale(dis1,k=1,eig=TRUE)$eig
plot(eigen_values1, ylim=c(-20,80),xlim=c(0,8), ylab="eigenvalue",xlab="number of dimensions")
title("Eigenvalue of distance matrix D1")

eigen_values2 <- cmdscale(dis2,k=1,eig=TRUE)$eig
plot(eigen_values2, ylim=c(-1,5),xlim=c(0,8),  ylab="eigenvalue",xlab="number of dimensions")
title("Eigenvalue of distance matrix D2")

eigen_values3 <- cmdscale(dis3,k=1,eig=TRUE)$eig
plot(eigen_values3,xlim=c(0,8),  ylab="eigenvalue",xlab="number of dimensions")
title("Eigenvalue of distance matrix D3")
## GOF
Gof1_1 <- cmdscale(dis1,k=1,eig=TRUE)$GOF
Gof1_2 <- cmdscale(dis1,k=2,eig=TRUE)$GOF
Gof1_3 <- cmdscale(dis1,k=3,eig=TRUE)$GOF
Gof1_4 <- cmdscale(dis1,k=4,eig=TRUE)$GOF

Gof2_1 <- cmdscale(dis2,k=1,eig=TRUE)$GOF
Gof2_2 <- cmdscale(dis2,k=2,eig=TRUE)$GOF
Gof2_3 <- cmdscale(dis2,k=3,eig=TRUE)$GOF
Gof2_4 <- cmdscale(dis2,k=4,eig=TRUE)$GOF

Gof3_1 <- cmdscale(dis3,k=1,eig=TRUE)$GOF
Gof3_2 <- cmdscale(dis3,k=2,eig=TRUE)$GOF
Gof3_3 <- cmdscale(dis3,k=3,eig=TRUE)$GOF
Gof3_4 <- cmdscale(dis3,k=4,eig=TRUE)$GOF

##
par(mfrow = c(3, 4))
plot(dis1, dist(model1_1), ylab="",xlab="")
plot(dis2, dist(model2_1), ylab="",xlab="")
plot(dis3, dist(model3_1), ylab="",xlab="")

plot(dis1, dist(model1_2), ylab="",xlab="")
plot(dis2, dist(model2_2), ylab="",xlab="")
plot(dis3, dist(model3_2), ylab="",xlab="")

plot(dis1, dist(model1_3), ylab="",xlab="")
plot(dis2, dist(model2_3), ylab="",xlab="")
plot(dis3, dist(model3_3), ylab="",xlab="")


plot(dis1, dist(model1_4), ylab="",xlab="")
plot(dis2, dist(model2_4), ylab="",xlab="")
plot(dis3, dist(model3_4), ylab="",xlab="")

m_dis1 <- as.matrix(dis1)
m_dis2 <- as.matrix(dis2)
m_dis3 <- as.matrix(dis3)

M_e1 <- as.matrix(eigen_values1)
M_e2 <- as.matrix(eigen_values2)
M_e3 <- as.matrix(eigen_values3)

par(mfrow = c(2, 3))
plot( model2_2[,1],normalized2_M[,1], ylab="Energy in D2",xlab="Model2_2_X")
plot(model2_2[,1], normalized2_M[,2], ylab="Protein in D2",xlab="Model2_2_X")

plot( model2_2[,1], normalized2_M[,3],ylab="Fat in D2",xlab="Model2_2_X")

plot( model2_2[,1], normalized2_M[,4],ylab="Calcium in D2",xlab="Model2_2_X")

plot( model2_2[,1], normalized2_M[,5],ylab="Iron in D2",xlab="Model2_2_X")

##
par(mfrow = c(2, 3))
plot( model2_2[,2],normalized2_M[,1], ylab="Energy in D2",xlab="Model2_2_Y")
plot(model2_2[,2], normalized2_M[,2], ylab="Protein in D2",xlab="Model2_2_Y")

plot( model2_2[,2], normalized2_M[,3],ylab="Fat in D2",xlab="Model2_2_Y")

plot( model2_2[,2], normalized2_M[,4],ylab="Calcium in D2",xlab="Model2_2_Y")

plot( model2_2[,2], normalized2_M[,5],ylab="Iron in D2",xlab="Model2_2_Y")



max_diffA1 <- max(abs(dis1-dist(model1_1)))
max_diffA2 <- max(abs(dis1-dist(model1_2)))
max_diffA3 <- max(abs(dis1-dist(model1_3)))
max_diffA4 <- max(abs(dis1-dist(model1_4)))
max_diffB1 <- max(abs(dis2-dist(model2_1)))
max_diffB2 <- max(abs(dis2-dist(model2_2)))
max_diffB3 <- max(abs(dis2-dist(model2_3)))
max_diffB4 <- max(abs(dis2-dist(model2_4)))

mean_diffA1 <- mean(abs(dis1-dist(model1_1)))
mean_diffA2 <- mean(abs(dis1-dist(model1_2)))
mean_diffA3 <- mean(abs(dis1-dist(model1_3)))
mean_diffA4 <- mean(abs(dis1-dist(model1_4)))
mean_diffB1 <- mean(abs(dis2-dist(model2_1)))
mean_diffB2 <- mean(abs(dis2-dist(model2_2)))
mean_diffB3 <- mean(abs(dis2-dist(model2_3)))
mean_diffB4 <- mean(abs(dis2-dist(model2_4)))

meanabs_diffA1 <- mean(abs((dis1-dist(model1_1))/dis1))
meanabs_diffA2 <- mean(abs((dis1-dist(model1_2))/dis1))
meanabs_diffA3 <- mean(abs((dis1-dist(model1_3))/dis1))
meanabs_diffA4 <- mean(abs((dis1-dist(model1_4))/dis1))
meanabs_diffB1 <- mean(abs((dis2-dist(model2_1))/dis2))
meanabs_diffB2 <- mean(abs((dis2-dist(model2_2))/dis2))
meanabs_diffB3 <- mean(abs((dis2-dist(model2_3))/dis2))
meanabs_diffB4 <- mean(abs((dis2-dist(model2_4))/dis2))

corE <- cor(model2_2[,2],normalized2_M[,1])
corP <- cor(model2_2[,2],normalized2_M[,2])
corF <- cor(model2_2[,2],normalized2_M[,3])
corC <- cor(model2_2[,2],normalized2_M[,4])
corI <- cor(model2_2[,2],normalized2_M[,5])


