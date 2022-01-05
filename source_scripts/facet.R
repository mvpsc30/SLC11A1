data$group=ifelse(data$SLC11A1>median(data$SLC11A1),"high","low")
data$age=ifelse(data$Age>40,">40","≤40")
library(survminer)
fit<- survfit(Surv(survival, status) ~ group,
              data = data)
colnames(data)
grade=ggsurvplot_facet(fit,data,
                            facet.by = c("Grade"),
                            risk.table = F,       # show risk table.
                            pval = TRUE, 
                            pval.method = F,# show p-value of log-rank test.
                            conf.int = TRUE,         # show confidence intervals for 
                            # point estimates of survival curves.
                            palette = c( "#E7B800","#00AFBB"),
                            title="Grade",
                            xlab = "Time in months",   # customize X axis label.
                            break.time.by = 50,     # break X axis in time intervals by 500.
                            ggtheme = theme_bw(), # customize plot and risk table with a theme.
                            conf.int.style = "step",  # customize style of confidence intervals
                            #surv.median.line = "hv",  # add the median survival pointer.
                            legend.labs = c("High", "Low")    # change legend labels.)
                           )

IDH_status=ggsurvplot_facet(fit,data,
                            facet.by = c("IDH_status"),
                            risk.table = F,       # show risk table.
                            pval = TRUE, 
                            pval.method = F,# show p-value of log-rank test.
                            conf.int = TRUE,         # show confidence intervals for 
                            # point estimates of survival curves.
                            palette = c( "#E7B800","#00AFBB"),
                            title="IDH status",
                            xlab = "Time in months",   # customize X axis label.
                            break.time.by = 50,     # break X axis in time intervals by 500.
                            ggtheme = theme_bw(), # customize plot and risk table with a theme.
                            conf.int.style = "step",  # customize style of confidence intervals
                            #surv.median.line = "hv",  # add the median survival pointer.
                            legend.labs = c("High", "Low")    # change legend labels.)
                            )

co_del_1p19q=ggsurvplot_facet(fit,data,
                            facet.by = c("co_del_1p19q"),
                            risk.table = F,       # show risk table.
                            pval = TRUE, 
                            pval.method = F,# show p-value of log-rank test.
                            conf.int = TRUE,         # show confidence intervals for 
                            # point estimates of survival curves.
                            palette = c( "#E7B800","#00AFBB"),
                            title="1p/19q codeletion",
                            xlab = "Time in months",   # customize X axis label.
                            break.time.by = 50,     # break X axis in time intervals by 500.
                            ggtheme = theme_bw(), # customize plot and risk table with a theme.
                            conf.int.style = "step",  # customize style of confidence intervals
                            #surv.median.line = "hv",  # add the median survival pointer.
                            legend.labs = c("High", "Low")    # change legend labels.)
                            )

age=ggsurvplot_facet(fit,data,
                     facet.by = c("age"),
                     risk.table = F,       # show risk table.
                     pval = TRUE, 
                     pval.method = F,# show p-value of log-rank test.
                     conf.int = TRUE,         # show confidence intervals for 
                     # point estimates of survival curves.
                     palette = c( "#E7B800","#00AFBB"),
                     title="Age",
                     xlab = "Time in months",   # customize X axis label.
                     break.time.by = 50,     # break X axis in time intervals by 500.
                     ggtheme = theme_bw(), # customize plot and risk table with a theme.
                     conf.int.style = "step",  # customize style of confidence intervals
                     #surv.median.line = "hv",  # add the median survival pointer.
                     legend.labs = c("High", "Low")    # change legend labels.)
)


(grade+IDH_status+plot_layout(widths = c(3,2)))/
(co_del_1p19q+age+plot_layout(widths = c(1,1)))+plot_annotation(tag_levels = "A")

ggsave("plot3.pdf",height =8,width = 12)
