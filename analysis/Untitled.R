#=========================================================================================
#========== Variance partitioning and correlations
#=========================================================================================

### Pelagic

# pel_par_names <- c(names(pelagic_rds$fit)[str_detect(names(pelagic_rds$fit), "b_lz")],
#                    names(pelagic_rds$fit)[str_detect(names(pelagic_rds$fit), "o_lz")],
#                    names(pelagic_rds$fit)[str_detect(names(pelagic_rds$fit), "r_lz")])
# 
# pel_pars_full <- pelagic_rds$fit %>%
#   rstan::extract(pars = pel_par_names) %>%
#   parallel::mclapply(as_tibble) %>%
#   bind_cols() %>%
#   set_names(pel_par_names) %>%
#   mutate(step = row_number()) %>%
#   gather(var, val, -step) %>%
#   mutate(par = strsplit(var, "\\[|\\]|,") %>% map_chr(~as.character(.x[1])),
#          id = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2]))) %>%
#   left_join(pelagic_rds$dd %>%
#               select(date, site, event) %>%
#               unique() %>%
#               mutate(id = row_number()))
# 
# 
# pel_pars_full %>%
#   ggplot(aes(val,
#              color = var))+
#   facet_wrap(~par)+
#   geom_density()+
#   scale_color_discrete(guide = "none")
# 
# pel_cor <- pel_pars_full %>%
#   split(.$step) %>%
#   parallel::mclapply(function(x){
#     x %>%
#       select(par, val, event) %>%
#       ungroup() %>%
#       arrange(par) %>%
#       pivot_wider(names_from = par,
#                   values_from = val) %>%
#       select(-event) %>%
#       cor() %>%
#       as_tibble() %>%
#       mutate(step = unique(x$step),
#              row = row_number())
#   }) %>%
#   bind_rows
# 
# pel_cor <- pel_cor %>%
#   pivot_longer(cols = c(b_lz, o_lz, r_lz)) %>%
#   group_by(row, name) %>%
#   summarize(lo = quantile(value, probs = 0.16),
#             mi = median(value),
#             hi = quantile(value, probs = 0.84)) %>%
#   ungroup() %>%
#   pivot_wider(names_from = name,
#               values_from = c(lo, mi, hi))
# 
# write_csv(pel_cor, "analysis/output/pel_cor.csv")
# 
# pel_var_part <- pel_pars_full %>%
#   split(.$step) %>%
#   parallel::mclapply(function(x_){
#     x_ %>%
#       split(.$par) %>%
#       lapply(function(xx_){
#         d_ = xx_ %>% select(step,par) %>% unique()
#         d_$f = anova(lm(val ~ site, data = xx_))[1,"F value"]
#         return(d_)}) %>% bind_rows()
#   }) %>% bind_rows()
# write_csv(pel_var_part, "analysis/output/pel_var_part.csv")

pel_cor <- read_csv("analysis/output/pel_cor.csv")
pel_var_part <- read_csv("analysis/output/pel_var_part.csv")

pel_var_part_sum <- pel_var_part  %>%
  group_by(par) %>%
  summarize(lo = quantile(log(f), probs = 0.16),
            mi = median(log(f)),
            hi = quantile(log(f), probs = 0.84))


### Benthic

# ben_par_names <- c(names(benthic_rds$fit)[str_detect(names(benthic_rds$fit), "b_lz")],
#                    names(benthic_rds$fit)[str_detect(names(benthic_rds$fit), "o_lz")],
#                    names(benthic_rds$fit)[str_detect(names(benthic_rds$fit), "r_lz")])
# ben_pars_full <- benthic_rds$fit %>%
#   rstan::extract(pars = ben_par_names) %>%
#   parallel::mclapply(as_tibble) %>%
#   bind_cols() %>%
#   set_names(ben_par_names) %>%
#   mutate(step = row_number()) %>%
#   gather(var, val, -step) %>%
#   mutate(par = strsplit(var, "\\[|\\]|,") %>% map_chr(~as.character(.x[1])),
#          id = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2]))) %>%
#   left_join(benthic_rds$dd %>%
#               select(date, site, event) %>%
#               unique() %>%
#               mutate(id = row_number()))
# 
# ben_pars_full %>%
#   ggplot(aes(val,
#              color = var))+
#   facet_wrap(~par)+
#   geom_density()+
#   scale_color_discrete(guide = "none")
# 
# ben_cor <- ben_pars_full %>%
#   split(.$step) %>%
#   parallel::mclapply(function(x){
#     x %>%
#       select(par, val, event) %>%
#       ungroup() %>%
#       arrange(par) %>%
#       pivot_wider(names_from = par,
#                   values_from = val) %>%
#       select(-event) %>%
#       cor() %>%
#       as_tibble() %>%
#       mutate(step = unique(x$step),
#              row = row_number())
#   }) %>%
#   bind_rows
# 
# ben_cor <- ben_cor %>%
#   pivot_longer(cols = c(b_lz, o_lz, r_lz)) %>%
#   group_by(row, name) %>%
#   summarize(lo = quantile(value, probs = 0.16),
#             mi = median(value),
#             hi = quantile(value, probs = 0.84)) %>%
#   ungroup() %>%
#   pivot_wider(names_from = name,
#               values_from = c(lo, mi, hi))
# 
# write_csv(ben_cor, "analysis/output/ben_cor.csv")
# 
# 
# ben_var_part <- ben_pars_full %>%
#   split(.$step) %>%
#   parallel::mclapply(function(x_){
#     x_ %>%
#       split(.$par) %>%
#       lapply(function(xx_){
#         d_ = xx_ %>% select(step,par) %>% unique()
#         d_$f = anova(lm(val ~ site, data = xx_))[1,"F value"]
#         return(d_)}) %>% bind_rows()
#   }) %>% bind_rows()
# write_csv(ben_var_part, "analysis/output/ben_var_part.csv")

ben_cor <- read_csv("analysis/output/ben_cor.csv")
ben_var_part <- read_csv("analysis/output/ben_var_part.csv")

ben_var_part_sum <- ben_var_part  %>%
  group_by(par) %>%
  summarize(lo = quantile(log(f), probs = 0.16),
            mi = median(log(f)),
            hi = quantile(log(f), probs = 0.84))

# combine
var_part_sum = bind_rows(pel_var_part_sum %>% mutate(type = "Pelagic"),
                         ben_var_part_sum %>% mutate(type = "Benthic")) %>%
  mutate(type = factor(type,
                       levels = c("Pelagic",
                                  "Benthic")),
         name = factor(par,
                       levels = c("o_lz",
                                  "b_lz",
                                  "r_lz"),
                       labels = c("Optimum PAR",
                                  "Max GPP",
                                  "Respiration")))

# plot
fig_var <- ggplot(data = var_part_sum,
                  aes(x = name,
                      y = mi,
                      color = type))+
  geom_hline(yintercept = 0,
             size = 0.25,
             linetype = 1,
             color = "gray30")+
  geom_point(size = 3, position = position_dodge(width = 0.25))+
  geom_errorbar(aes(ymin = lo,
                    ymax = hi),
                position = position_dodge(width = 0.25),
                width = 0.1)+
  scale_y_continuous(Log~"F"~ratio~(between~vs.~within~sites),
                     limits = c(-5, 5),
                     breaks = c(-5, -2.5, 0, 2.5, 5),
                     labels = c(-5, -2.5, "0", 2.5, 5))+
  scale_color_manual("",values = c("black","dodgerblue"))+
  scale_x_discrete("")+
  theme(legend.position = "top",
        panel.border = element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.text = element_blank(),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25))
fig_var

# cairo_pdf(file = "analysis/figures/fig_var.pdf",
#           width = 3.5, height = 2.5, family = "Arial")
# fig_var
# dev.off()

#=========================================================================================
