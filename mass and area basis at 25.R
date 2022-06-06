
library(car); library(lubridate); library(lme4);
library(MuMIn);library(car); library(multcomp); library(effects)
library(lme4); library(MuMIn); library(lmerTest); library(doBy); library(plotrix);library(sjPlot)

setwd("E:/WETFERT Analysis/Cleaned/Analysis")
aci<-read.csv("One at a time values_WETFERT.csv")
# aci$Date<-as.Date(aci$date.a, "%m/%d/%Y")
# aci$MoYr<-as.factor(floor_date(aci$Date, "month"))

setwd("E:/WETFERT Analysis/WETFERT Respiration/Analysis")
df<-read.csv("LPR WETFERT R Model.csv")
# df$Date<-as.Date(df$Date, "%m/%d/%Y")


# DF3<-merge(aci, df,by=c("date.c","UserIDs_in"),all=TRUE)
# levels(DF1$Code);levels(DF2$Code)
# write.csv(DF3,"photorespFERT.csv")
# levels(as.factor(DF1$Sampling))
# 
# levels(DF1$Sampling)
# str(DF1)

df<-read.csv("photorespFERT.csv")

# dum<-subset(df, Leaf.age == "new" & R_calc == "R.area" & MyTemp == "25" & keep == "Y")
# 
# m1<- lm(Asat ~ rdref*treatment,dum)
# Anova(m1)
# summary(m1)
# plot(allEffects(m1))
df$jmax.25mass<-as.numeric(df$jmax.25mass)
df$vcmax.25mass<-as.numeric(df$vcmax.25mass)

sdf<-summaryBy(rdref + q10 + LMA_g.m2 + jmax + vcmax + vcmax.25 + vcmax.25mass + jmax.25 + jmax.25mass + Asat + gs ~ MyTemp + Date + keep + R_calc + Leaf.age * Treatment, FUN = c(mean,std.error), na.rm = T, df)

sdf<-summaryBy(rdref + q10 + LMA_g.m2 + vcmax.25 + vcmax.25mass + jmax.25 + jmax.25mass + Asat + gs ~ MyTemp + Date + keep + R_calc + Leaf.age * Treatment, FUN = c(mean,std.error), na.rm = T, df)




sdf<-summaryBy(vcmax.25mass + jmax.25mass ~ MyTemp + Date + keep + R_calc + Leaf.age * Treatment, FUN = c(mean,std.error), na.rm = T, df)




tiff(file = "mass and are basis at 25.tiff", height = 8, width = 4, res = 600, units = "in", compression = "zip+p")
par(mfrow = c(2,1), omi = c(1, 0.5, 0.25, 0.1), mar = c(2.5,6,2.5,0.5))

#########################################################################

dum<-subset(df,R_calc == "R.mass" & Treatment == "c" & MyTemp == "25" & Leaf.age == "new" & keep == "Y" & Jkeep == "Y")
plot(jmax.25mass ~ vcmax.25mass, dum, pch = 16, col= "blue", ylim = c(0, 1200), xlim = c(0,500), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
# rect(xleft = as.Date("2020-08-17"), xright = as.Date("2020-08-18"), ybottom = -999, ytop = 999, col = "grey80", bty = "n", border = F)
# box()
dum<-subset(df,R_calc == "R.mass"  & Treatment == "c" & MyTemp == "25" & Leaf.age == "new" & keep == "Y" & Jkeep == "Y")
points(jmax.25mass ~ vcmax.25mass, dum, pch = 16, col= "black")
dum<-subset(df,R_calc == "R.mass"  & Treatment == "n" & MyTemp == "25" & Leaf.age == "new" & keep == "Y" & Jkeep == "Y" )
points(jmax.25mass ~ vcmax.25mass, dum, pch = 16, col= "dodgerblue")
dum<-subset(df,R_calc == "R.mass"  & Treatment == "p" & MyTemp == "25" & Leaf.age == "new" & keep == "Y" & Jkeep == "Y") 
points(jmax.25mass ~ vcmax.25mass, dum, pch = 16, col= "forestgreen")



dum$jmax.25mass<-as.factor(dum$jmax.25mass)
dum$vcmax.25mass<-as.factor(dum$vcmax.25mass)
dum<-subset(df,R_calc == "R.mass" & MyTemp == "25" & Leaf.age == "new" & keep == "Y" & Jkeep == "Y")
m1<-lm(jmax.25mass ~ vcmax.25mass, dum)
anova(m1)
summary(m1)
p1<-plot_model(m1, type= c("pred"), terms= c("vcmax.25mass"))
new<-as.data.frame(p1$data)
m2<-lm(predicted~x, new)
ablineclip(m2,x1=min(dum$vcmax.25mass,na.rm = TRUE),x2=max(dum$vcmax.25mass,na.rm = TRUE), lwd=2 )


axis(1, at = seq(0,500,100), las = 1, cex.axis = 0.85, labels =T)
axis(2, at = seq(0,1200,200), las = 2, cex.axis = 1.2)
# mtext(side = 2, expression(italic(R)[mass]*''^25^degree*""^C~(n*mol~g^-1~s^-1)), cex = 1.4, padj = -2, outer= F)

mtext(side = 2, expression(italic(J)[max]*''^25^degree*""^C~(n*mol~g^-1~s^-1)), cex = 1.2, padj = -2.2, outer= F)
mtext(side = 2, expression("Mass Basis"), cex = 1.4, padj = -8, outer= F)
# legend("topright", c("vcmax.25*  ","T  ","V x T  "), bty = "n", cex = 1)
legend("topright", c(expression(italic(V)[cmax^25]~'***',"Treatment", italic(V)[cmax^25]~' x Trt' )), bty = "n", cex = 0.8)
mtext(side = 1, expression(italic(V)[cmax]*''^25^degree*""^C~(n*mol~g^-1~s^-1)), cex = 1.2, padj = 2, outer= F)
text(100, 60, expression(italic(r)^2~'= 0.63'),cex=1.2)
legend("bottomright", c("Control","Nitrogen","Phosphorus"), col=c("black", "dodgerblue", "forestgreen"), pch=16, cex = 0.8, horiz = F, bty='n')

legend("topleft", c("(a)"), bty = "n", cex = 1.2)

###########################################################################

dum<-subset(df,R_calc == "R.mass" & Treatment == "c" & MyTemp == "25" & Leaf.age == "new" & keep == "Y" & Jkeep == "Y")
plot(jmax.25 ~ vcmax.25, dum, pch = 16, col= "blue", ylim = c(40, 200), xlim = c(0,100), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
# rect(xleft = as.Date("2020-08-17"), xright = as.Date("2020-08-18"), ybottom = -999, ytop = 999, col = "grey80", bty = "n", border = F)
# box()
dum<-subset(df,R_calc == "R.mass"  & Treatment == "c" & MyTemp == "25" & Leaf.age == "new" & keep == "Y" & Jkeep == "Y")
points(jmax.25 ~ vcmax.25, dum, pch = 16, col= "black")
dum<-subset(df,R_calc == "R.mass"  & Treatment == "n" & MyTemp == "25" & Leaf.age == "new" & keep == "Y" & Jkeep == "Y")
points(jmax.25 ~ vcmax.25, dum, pch = 16, col= "dodgerblue")
dum<-subset(df,R_calc == "R.mass"  & Treatment == "p" & MyTemp == "25" & Leaf.age == "new" & keep == "Y" & Jkeep == "Y") 
points(jmax.25 ~ vcmax.25, dum, pch = 16, col= "forestgreen")


dum<-subset(df,R_calc == "R.mass" & MyTemp == "25" & Leaf.age == "new" & keep == "Y" & Jkeep == "Y")
m1<-lm(jmax.25 ~ vcmax.25*treatment, dum)
anova(m1)
summary(m1)
p1<-plot_model(m1, type= c("pred"), terms= c("vcmax.25"))
new<-as.data.frame(p1$data)
m2<-lm(predicted~x, new)
ablineclip(m2,x1=min(dum$vcmax.25,na.rm = TRUE),x2=max(dum$vcmax.25,na.rm = TRUE), lwd=2 )

mtext(side = 1, expression(italic(V)[cmax]*''^25^degree*""^C~(mu*mol~m^-2~s^-1)), cex = 1.2, padj = 2, outer= F)
mtext(side = 2, expression("Area Basis"), cex = 1.4, padj = -8, outer= F)
mtext(side = 2, expression(italic(J)[max]*''^25^degree*""^C~(mu*mol~m^-2~s^-1)), cex = 1.2, padj = -2, outer= F)
axis(1, at = seq(0,100,20), las = 1, cex.axis = 0.85, labels =T)
axis(2, at = seq(40,200,40), las = 2, cex.axis = 1.2, labels = T)

legend("topright", c(expression(italic(V)[cmax^25]~'***',"Treatment", italic(V)[cmax^25]~' x Trt' )), bty = "n", cex = 0.8)
text(60, 50, expression(italic(r)^2~'= 0.61'),cex=1.2)
# mtext(side = 2, expression(italic(R)[mass]*''^25^degree*""^C~(n*mol~g^-1~s^-1)), cex = 1.4, padj = -2, outer= F)

# mtext(side = 2, expression(italic(J)[max]~(mu*mol~m^-2~s^-1)), cex = 1.4, padj = -2, outer= F)

# legend("topright", c("vcmax.25*  ","T  ","V x T  "), bty = "n", cex = 1)

# legend("bottomleft", c("Control","Nitrogen","Phosphorus"), col=c("black", "dodgerblue", "forestgreen"), pch=1, lty = 1, cex = 1.05, horiz = F, bty='n')

legend("topleft", c("(b)"), bty = "n", cex = 1.2)


##############################################################################

dev.off() 
