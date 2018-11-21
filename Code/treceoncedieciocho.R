d <- read.table("Data_FG/FG.txt", header = TRUE, sep = "\t")

levels(d$time)
d$time <- as.factor(d$time)
levels(d$time)

levels(d$date)
levels(d$date) <- c("Jun12", "Jul13", "May13", "May14", "Oct13", "Sep12")
d$date <- factor(d$date, levels = c("Jun12", "Sep12", "May13", "Jul13", "Oct13", "May14"))
levels(d$date)

levels(d$block)
d$block <- as.factor(d$block)
levels(d$block)

levels(d$treatment)
d$treatment <- factor(d$treatment, levels = c("control", "Cb", "Cd", "Ci", "Ee", "Pg", "Pp", "Cb.Cd.Pg",
                                                "Cb.Ci.Pg", "Cd.Pp.Ee","Ci.Ee.Pg", "Ci.Pp.Ee", "Ci.Pp.Pg", "6sp"))
levels(d$treatment)

levels(d$datapoint)
d$datapoint <- as.factor(d$datapoint)
levels(d$datapoint)

levels(d$Focal)
boxplot(d$Cover)

d$total <- apply(d[, seq(15, 91, by = 2)], 1, function (x) sum(x))
boxplot(d$total)

write.table(d, file = "Data_FG/FG.txt", sep = "\t", row.names = FALSE)

