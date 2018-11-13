df <- read.table("Data_FG/FG.txt", header = TRUE, sep = "\t")

levels(df$time)
df$time <- as.factor(df$time)
levels(df$time)

levels(df$date)
levels(df$date) <- c("Jun12", "Jul13", "May13", "May14", "Oct13", "Sep12")
df$date <- factor(df$date, levels = c("Jun12", "Sep12", "May13", "Jul13", "Oct13", "May14"))
levels(df$date)

levels(df$block)
df$block <- as.factor(df$block)
levels(df$block)

levels(df$treatment)
df$treatment <- factor(df$treatment, levels = c("control", "Cb", "Cd", "Ci", "Ee", "Pg", "Pp", "Cb.Cd.Pg",
                                                "Cb.Ci.Pg", "Cd.Pp.Ee","Ci.Ee.Pg", "Ci.Pp.Ee", "Ci.Pp.Pg", "6sp"))
levels(df$treatment)

levels(df$datapoint)
df$datapoint <- as.factor(df$datapoint)
levels(df$datapoint)

levels(df$Focal)
boxplot(df$Cover)

df$total <- apply(df[, seq(15, 91, by = 2)], 1, function (x) sum(x))
df$total

colnames(df)
