## Question 1
setwd("D:/研一下/生物统计/生统考试")
q1 <- read.table("Q1_gene_data.txt")
str(q1)
G3 <- as.numeric(q1[3,])
plot(density(G3),xlab = "Sample",ylab = "Gene expression")
min(G3);max(G3);mean(G3);var(G3);

boxplot(q1)
q1_clean <- q1
q1_clean[q1_clean>100] <- NA
boxplot(q1_clean)


## Question 2
q2 <- read.table("Q2_skin_test.txt",header = T,row.names = 1)
attach(q2)
shapiro.test(right);shapiro.test(left)
var.test(right,left,alternative = "two.sided")
wilcox.test(right,left,paired = T,alternative = "two.sided",exact = F)
t.test(right,left,paired = T,var.equal = T)


# Question 3
q3 <- read.csv("Q3_SBP.csv")
attach(q3)
shapiro.test(follow_up_1)
shapiro.test(follow_up_2)
var.test(follow_up_1,follow_up_2,alternative = "two.sided")
wilcox.test(follow_up_1,follow_up_2,paired = T,alternative = "two.sided",exact = F)
t.test(follow_up_1,follow_up_2,paired = T,var.equal = T)


# Question 4
q4 <- read.csv("Q4_antibiotic.csv")
attach(q4)
q4_aov <- aov(percent~antibiotic)
summary(q4_aov)
TukeyHSD(q4_aov)

#Calculate the ANOVA Manually
total_SS <- sum((percent - mean(percent))^2)
group_var <- aggregate(percent,by = list(antibiotic),FUN = var)
group_num <- aggregate(percent,by = list(antibiotic),FUN = length)
within_SS <- sum(group_var$x * (group_num$x-1))
between_SS <- total_SS - within_SS

df_between <- length(group_var$Group.1) - 1
df_within <- length(percent) - length(group_var$Group.1)

between_MS <- between_SS/df_between
within_MS <- within_SS/df_within
F_ratio <- between_MS/within_MS
p_value <- pf(F_ratio,df_between,df_within,lower.tail = FALSE)
