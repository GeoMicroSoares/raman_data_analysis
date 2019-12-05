fitpar_treated=fitpar %>% 
  filter(str_detect(param, "FWHM")) %>% 
  separate(name, c("strain","conc","rep"),
           sep = "_")

fitpar_treated$conc <- factor(fitpar_treated$conc,
                             levels = c("05mM", "1mM", "5mM", "10mM"))

ggplot(fitpar_treated,
       aes(conc,fit,
           colour=rep)) +
  geom_jitter(size=4,alpha=.6,
              width = .3, height = .3) +
  facet_wrap(param~strain, 
             nrow = 3, scales = "free_y") +
  labs(y="FWHM (nm)",
       colour="Replicate") +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"),
        strip.text.x = element_text(face = "bold"),
        panel.grid.major.y = element_line(size=1.5),
        panel.grid.minor.y = element_line(size=.5)) +
  ggtitle("Full-width half maxima for PBH481 and PBS1",
          subtitle = expression(paste("Experiment 1: Increasing concentrations of ",
                                      Na[2],SeO[3])))
ggsave("fwhm_pbhps_exp1.png",
       height = 6, width = 6)
