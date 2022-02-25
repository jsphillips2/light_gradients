#=========================================================================================
#========== Preliminaries
#=========================================================================================

# load packages
library(tidyverse)
library(lemon)
library(cowplot)
library(lubridate)
library(viridisLite)

# import data
pelagic_rds <- read_rds("analysis/output/pelagic.rds")
benthic_rds <- read_rds("analysis/output/benthic.rds")
options(mc.cores = parallel::detectCores()-2)

# create date levels
date_levels <- bind_rows(pelagic_rds$dd,
                         benthic_rds$dd) %>%
  mutate(date_level = paste(month(as_date(date),
                                  label = T),
                            year(as_date(date))) %>% 
           factor(levels = c("Jun 2018",
                             "Jul 2018",
                             "Aug 2018",
                             "Jul 2019",
                             "Aug 2019"),
                  labels = c("Jun '18",
                             "Jul '18",
                             "Aug '18",
                             "Jul '19",
                             "Aug '19"))) %>%
  select(date, date_level) %>%
  unique()

# define date level colors
date_colors <- viridis(5,
                       begin = 0,
                       end = 0.9,
                       option = "viridis")
names(date_colors) <- levels(date_levels$date_level)

# define site levels
site_levels <- bind_rows(pelagic_rds$dd,
                         benthic_rds$dd) %>%
  mutate(site_level = factor(site,
                             levels = c("reyk",
                                        "st33",
                                        "e5"),
                             labels = c("north",
                                        "south",
                                        "east"))) %>%
  select(site, site_level) %>%
  unique()

# define site colors
site_colors <- viridis(3,
                       begin = 0.8,
                       end = 0.1,
                       option = "magma")
names(site_colors) <- levels(site_levels$site_level)

# set theme
theme_set(theme_bw() %+replace%
            theme(panel.grid = element_blank(),
                  strip.background = element_blank(),
                  plot.margin = margin(t = 1,
                                       r = 1,
                                       b = 1,
                                       l = 1),
                  legend.margin = margin(t = 0,
                                         r = 0,
                                         b = 0,
                                         l = -4),
                  legend.text = element_text(size = 8),
                  axis.text = element_text(size = 10, color="black",family = "sans"),
                  axis.title = element_text(size =10),
                  axis.title.y = element_text(angle = 90, margin=margin(0,5,0,0)),
                  axis.title.x = element_text(margin = margin(5,0,0,0)),
                  panel.spacing = unit(0.1, "lines"),
                  axis.ticks = element_line(size = 0.25)))

#=========================================================================================





#=========================================================================================
#========== PI curves: Pelagic
#=========================================================================================

# extract fits
pelagic_nep <- pelagic_rds$fit_summary %>%
  filter(str_detect(.$var, "nep_sum")) %>%
  mutate(name = strsplit(var, "\\[|\\]|,") %>% map_chr(~.x[1]),
         id = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2]))) %>%
  filter(name %in% c("nep_sum")) %>%
  rename(lo = `16%`,
         mi = `50%`,
         hi = `84%`) %>%
  select(id, lo, mi, hi) %>%
  unique() %>%
  left_join(pelagic_rds$dd_sum %>%
              mutate(id = row_number())) %>%
  select(-id) %>%
  unique() %>%
  left_join(date_levels) %>%
  left_join(site_levels)


# plot
fig_nep_pel <- ggplot(data = pelagic_rds$dd %>%
                        left_join(date_levels)  %>%
                        left_join(site_levels),
                      aes(x = par, 
                          y = do_flux, 
                          color = date_level))+
  facet_rep_wrap(~site_level, 
                 nrow = 2,
                 strip.position = "top")+
  geom_hline(yintercept = 0,
             size = 0.25,
             linetype = 1,
             color = "gray30")+
  geom_ribbon(data = pelagic_nep,
              aes(x = par,
                  ymin = lo,
                  ymax = hi,
                  fill = date_level),
              linetype = 0,
              alpha = 0.2,
              inherit.aes = F)+
  geom_point(size = 1,
             alpha = 0.75,
             stroke = 0.25)+
  geom_line(data = pelagic_nep,
            aes(y = mi))+
  scale_color_manual("",
                     values = date_colors)+
  scale_fill_manual(values = date_colors,
                    guide = "none")+
  scale_y_continuous(Net~production~(mg~O[2]~m^{-2}~h^{-h}),
                     breaks = c(-0.06, -0.03, 0, 0.03, 0.06),
                     labels = 1000*c(-0.06, -0.03, 0, 0.03, 0.06),
                     limits = c(-0.063, 0.063))+
  scale_x_continuous(PAR~(mu*mol-photons~m^-{2}~s^{-1}),
                     breaks = c(0, 200, 400, 600, 800),
                     limits = c(0, 800))+
  theme(panel.border = element_blank(),
        panel.spacing = unit(0.75, "lines"),
        strip.text = element_text(margin = margin(t = 0,
                                                  b = 0,
                                                  l = 0,
                                                  r = -1),
                                  vjust = 0.25),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25),
        plot.margin = margin(t = 1,
                             b = 1,
                             l = 0,
                             r = 10),
        legend.position = "top",
        legend.margin = margin(t = 1,
                               b = 10,
                               l = 1,
                               r = 25),
        legend.text = element_text(size = 8,
                                   margin=margin(l = -8)),
        legend.key.size = unit(1, "lines"),
        legend.spacing.y = unit(0, "lines"),
        legend.spacing.x = unit(0.6, "lines"))

# examine
fig_nep_pel

# cairo_pdf(file = "analysis/figures/fig_nep_pel.pdf",
#           width = 3.5, height = 3.9, family = "Arial")
# fig_nep_pel
# dev.off()

#=========================================================================================





#=========================================================================================
#========== PI curves: Pelagic
#=========================================================================================

# extract fits
benthic_nep <- benthic_rds$fit_summary %>%
  filter(str_detect(.$var, "nep_sum")) %>%
  mutate(name = strsplit(var, "\\[|\\]|,") %>% map_chr(~.x[1]),
         id = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2]))) %>%
  filter(name %in% c("nep_sum")) %>%
  rename(lo = `16%`,
         mi = `50%`,
         hi = `84%`) %>%
  select(id, lo, mi, hi) %>%
  unique() %>%
  left_join(benthic_rds$dd_sum %>%
              mutate(id = row_number())) %>%
  select(-id) %>%
  unique() %>%
  left_join(date_levels)  %>%
  left_join(site_levels)

# benthic
fig_nep_ben <- ggplot(data = benthic_rds$dd %>% 
                        left_join(date_levels) %>%
                        left_join(site_levels),
                      aes(x = par, 
                          y = do_flux, 
                          color = date_level))+
  facet_rep_wrap(~site_level, 
                 nrow = 3,
                 strip.position = "top")+
  geom_hline(yintercept = 0,
             size = 0.25,
             linetype = 1,
             color = "gray30")+
  geom_ribbon(data = benthic_nep,
              aes(x = par,
                  ymin = lo,
                  ymax = hi,
                  fill = date_level),
              linetype = 0,
              alpha = 0.2,
              inherit.aes = F)+
  geom_point(size = 1,
             alpha = 0.75,
             stroke = 0.25)+
  geom_line(data = benthic_nep,
            aes(y = mi))+
  scale_color_manual("",
                     values = date_colors)+
  scale_fill_manual(values = date_colors,
                    guide = "none")+
  scale_y_continuous(Net~production~(mg~O[2]~m^{-2}~h^{-h}),
                     breaks = c(-0.28, -0.14, 0, 0.14, 0.28),
                     labels = 1000 * c(-0.28, -0.14, 0, 0.14, 0.28),
                     limits = c(-0.28, 0.28))+
  scale_x_continuous(PAR~(mu*mol-photons~m^-{2}~s^{-1}),
                     breaks = c(0, 250, 500, 750, 1000),
                     limits = c(0, 1000))+
  theme(panel.border = element_blank(),
        panel.spacing = unit(0.75, "lines"),
        strip.text = element_text(margin = margin(t = 0,
                                                  b = 0,
                                                  l = 0,
                                                  r = -1),
                                  vjust = 0.25),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25),
        plot.margin = margin(t = 1,
                             b = 1,
                             l = 0,
                             r = 10),
        legend.position = "top",
        legend.margin = margin(t = 1,
                               b = 10,
                               l = 1,
                               r = 30),
        legend.text = element_text(size = 8,
                                   margin=margin(l = -8)),
        legend.key.size = unit(1, "lines"),
        legend.spacing.y = unit(0, "lines"),
        legend.spacing.x = unit(0.6, "lines"))

# examine
fig_nep_ben

# cairo_pdf(file = "analysis/figures/fig_nep_ben.pdf",
#           width = 3.5, height = 5.8, family = "Arial")
# fig_nep_ben
# dev.off()

#=========================================================================================





#=========================================================================================
#========== Parameters through time
#=========================================================================================

##### Extract parameters

# pelagic
pel_pars <- pelagic_rds$fit_summary %>%
  filter(var %in% c(paste0("b[",1:11,"]"),
                    paste0("o[",1:11,"]"),
                    paste0("r[",1:11,"]"))) %>%
  mutate(name = strsplit(var, "\\[|\\]|,") %>% map_chr(~as.character(.x[1])),
         id = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2]))) %>%
  rename(lo = `16%`,
         mi = `50%`,
         hi = `84%`) %>%
  select(name, id, lo, mi, hi) %>%
  left_join(pelagic_rds$dd %>%
              select(date, site, event) %>%
              unique() %>%
              mutate(id = row_number())) %>%
  select(-id) %>% 
  unique() %>%
  left_join(date_levels) %>%
  left_join(site_levels) %>%
  mutate(id = row_number(),
         id = as.numeric(date_level),
         id = id + ifelse(site_level == "north", -0.1, 
                          ifelse(site_level == "south", 0.1, 0)),
         year = year(as_date(date)))  %>%
  ungroup()


# benthic
ben_pars <-  benthic_rds$fit_summary %>%
  filter(var %in% c(paste0("b[",1:11,"]"),
                    paste0("o[",1:11,"]"),
                    paste0("r[",1:11,"]"))) %>%
  mutate(name = strsplit(var, "\\[|\\]|,") %>% map_chr(~as.character(.x[1])),
         id = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2]))) %>%
  rename(lo = `16%`,
         mi = `50%`,
         hi = `84%`) %>%
  select(name, id, lo, mi, hi) %>%
  left_join(benthic_rds$dd %>%
              select(date, site, event) %>%
              unique() %>%
              mutate(id = row_number())) %>%
  select(-id) %>% 
  unique() %>%
  left_join(date_levels) %>%
  left_join(site_levels) %>%
  group_by(name) %>%
  arrange(site_level, event) %>%
  mutate(id = row_number(),
         id = as.numeric(date_level),
         id = id + ifelse(site_level == "north", -0.1, 
                          ifelse(site_level == "south", 0.1, 0.3)),
         year = year(as_date(date)))  %>%
  ungroup()

# plot function
fig_par_fun <- function(data_,
                        y_scale = 1,
                        y_title_,
                        y_breaks_,
                        y_limits_,
                        x_breaks_,
                        x_label_) {
  ggplot(data = data_,
         aes(x = id, 
             y = y_scale * mi, 
             color = site_level,
             group = interaction(site_level,year)))+
    geom_line()+
    geom_errorbar(aes(ymin = y_scale * lo,
                      ymax = y_scale * hi),
                  width = 0,
                  size = 0.3)+
    geom_point(size = 1.5)+
    scale_color_manual(values = site_colors,
                       guide = "none")+
    scale_fill_manual(values = site_colors,
                      guide = "none")+
    scale_y_continuous(y_title_,
                       breaks = y_breaks_,
                       limits = y_limits_)+
    scale_x_continuous("",
                       breaks = x_breaks_,
                       labels = x_label_)+
    theme(panel.border = element_blank(),
          panel.spacing = unit(0, "lines"),
          axis.line.x = element_line(size = 0.25),
          axis.line.y = element_line(size = 0.25))
}

# plot pelagic panels 

fig_par_pel_o <- fig_par_fun(data_ = pel_pars %>% filter(name == "o"),
                             y_title_ = "Optimum PAR",
                             y_breaks = c(0, 200, 400),
                             y_limits = c(0, 400),
                             x_breaks_ = 1:5,
                             x_label_ = NULL)+
  theme(plot.margin = margin(t = 0,
                             b = 0,
                             l = 1,
                             r = 6))

fig_par_pel_b <- fig_par_fun(data_ = pel_pars %>% filter(name == "b"),
                             y_scale = 1000,
                             y_title_ = "Max GPP",
                             y_breaks = c(0, 40, 80),
                             y_limits = c(0, 80),
                             x_breaks_ = 1:5,
                             x_label_ = NULL)+
  theme(plot.margin = margin(t = 0,
                             b = 0,
                             l = 1,
                             r = 6))

fig_par_pel_r <- fig_par_fun(data_ = pel_pars %>% filter(name == "r"),
                             y_scale = 1000,
                             y_title_ = "Respiration",
                             y_breaks = c(0, 15, 30),
                             y_limits = c(0, 30),
                             x_breaks_ = 1:5,
                             x_label_ = unique(date_levels$date_level))+
  theme(plot.margin = margin(t = 0,
                             b = 0,
                             l = 1,
                             r = 6),
        axis.text.x = element_text(angle = 40, 
                                   vjust = 0.75, 
                                   hjust=0.9))

fig_par_pel <- plot_grid(fig_par_pel_o,
                         fig_par_pel_b,
                         fig_par_pel_r,
                         nrow = 3,
                         align = "v")

# plot benthic panels 

fig_par_ben_o <- fig_par_fun(data_ = ben_pars %>% filter(name == "o"),
                             y_title_ = "",
                             y_breaks = c(0, 500, 1000),
                             y_limits = c(0,  1000),
                             x_breaks_ = 1:5,
                             x_label_ = NULL)+
  theme(plot.margin = margin(t = 0,
                             b = 0,
                             l = -4,
                             r = 5))

fig_par_ben_b <- fig_par_fun(data_ = ben_pars %>% filter(name == "b"),
                             y_scale = 1000,
                             y_title_ = "",
                             y_breaks = c(0, 200, 400),
                             y_limits = c(0, 400),
                             x_breaks_ = 1:5,
                             x_label_ = NULL)+
  theme(plot.margin = margin(t = 0,
                             b = 0,
                             l = -4,
                             r = 5))

fig_par_ben_r <- fig_par_fun(data_ = ben_pars %>% filter(name == "r"),
                             y_scale = 1000,
                             y_title_ = "",
                             y_breaks = c(0, 150, 300),
                             y_limits = c(0, 300),
                             x_breaks_ = 1:5,
                             x_label_ = unique(date_levels$date_level))+
  theme(plot.margin = margin(t = 0,
                             b = 0,
                             l = -4,
                             r = 5),
        axis.text.x = element_text(angle = 40, 
                                   vjust = 0.75, 
                                   hjust=0.9))

fig_par_ben <- plot_grid(fig_par_ben_o,
                         fig_par_ben_b,
                         fig_par_ben_r,
                         nrow = 3,
                         align = "v")

fig_par_comb <- plot_grid(fig_par_pel,
                          fig_par_ben,
                          ncol= 2,
                          align = "v",
                          rel_widths = c(0.95, 1),
                          labels = c("Pelagic",
                                     "Benthic"),
                          label_size = 10,
                          label_fontface = "plain",
                          hjust = c(-2.25, -1.75),
                          vjust = c(-1, -1))

fig_par_comb

# dummy plot for color guides
fig_par_color <- get_legend(ggplot(data = ben_pars  %>%
                                     left_join(site_levels),
                                   aes(x = id, 
                                       y = 1000 * mi, 
                                       color = site_level,
                                       group = interaction(site_level,year)))+
                              geom_line()+
                              geom_errorbar(aes(ymin = 1000 * lo,
                                                ymax = 1000 * hi),
                                            width = 0,
                                            size = 0.3)+
                              geom_point(size = 1.5)+
                              scale_color_manual("",
                                                 values = site_colors)+
                              theme(legend.position = "top",
                                    legend.margin = margin(t = 1,
                                                           b = 10,
                                                           l = 1,
                                                           r = 1),
                                    legend.text = element_text(size = 8,
                                                               margin=margin(l = -8)),
                                    legend.key.size = unit(1, "lines"),
                                    legend.spacing.y = unit(0, "lines"),
                                    legend.spacing.x = unit(0.6, "lines")))

fig_par <- plot_grid(fig_par_color,
                     fig_par_comb,
                     rel_heights = c(0.2, 1),
                     nrow = 2)
fig_par
# cairo_pdf(file = "analysis/figures/fig_par.pdf",
#           width = 3.5, height = 5, family = "Arial")
# fig_par
# dev.off()

#=========================================================================================




#=========================================================================================
#========== Random effects SDs
#=========================================================================================

# pelagic
re_sd_pel <- pelagic_rds$fit_summary %>%
  filter(var %in% c("lo_s[1]",
                    "lb_s[1]",
                    "lr_s[1]")) %>%
  mutate(name = strsplit(var, "\\[|\\]|,") %>% map_chr(~as.character(.x[1])),
         id = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2]))) %>%
  rename(lo = `16%`,
         mi = `50%`,
         hi = `84%`) %>%
  select(name, id, lo, mi, hi) 

# benthic
re_sd_ben <- benthic_rds$fit_summary %>%
  filter(var %in% c("lo_s[1]",
                    "lb_s[1]",
                    "lr_s[1]")) %>%
  mutate(name = strsplit(var, "\\[|\\]|,") %>% map_chr(~as.character(.x[1])),
         id = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2]))) %>%
  rename(lo = `16%`,
         mi = `50%`,
         hi = `84%`) %>%
  select(name, id, lo, mi, hi)

# combine
re_sd <- bind_rows(re_sd_pel %>% mutate(type = "Pelagic"),
                   re_sd_ben %>% mutate(type = "Benthic")) %>%
  mutate(type = factor(type,
                       levels = c("Pelagic",
                                  "Benthic")),
         name = factor(name,
                       levels = c("lo_s",
                                  "lb_s",
                                  "lr_s"),
                       labels = c("Optimum PAR",
                                  "Max GPP",
                                  "Respiration")))

# plot
fig_sd <- ggplot(data = re_sd,
                 aes(x = name,
                     y = mi,
                     color = type))+
  geom_point(size = 3, position = position_dodge(width = 0.25))+
  geom_errorbar(aes(ymin = lo,
                    ymax = hi), 
                position = position_dodge(width = 0.25),
                width = 0.1)+
  scale_y_continuous(Random~effect~SD~(log-scale),
                     limits = c(0, 1.2),
                     breaks = c(0, 0.4, 0.8, 1.2))+
  scale_color_manual("",values = c("black","dodgerblue"))+
  scale_x_discrete("")+
  theme(legend.position = "top",
        panel.border = element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.text = element_blank(),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25))
fig_sd

# cairo_pdf(file = "analysis/figures/fig_sd.pdf",
#           width = 3.5, height = 2.5, family = "Arial")
# fig_sd
# dev.off()
#=========================================================================================





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





#=========================================================================================
#========== Approximate half saturation
#=========================================================================================

### Direct solution

# half saturation would be given by 0.5 = l/w * exp(1 - l/w)
# constraint 0 < l < w 
# solution cannot be expressed in elementary functions:
# half sat = -omega W(-1/2e) where W is product log function
# approximates as 0.23961 * omega



### Approximate with hyperbolic tangent

# half sat for hyperbolic tangent is beta / alpha tanh-1(0.5)
# initial slope of photo-inhibition model is alpha = beta * e / omega
# so, half sat equals omega tanh-1(0.5) / e
# this in turn equals 0.2020784 * omega, which is pretty close to the direct solution



### Approximate with Michaelis-Menton

# initial slope of MM model is alpha = beta / k
# if set alphas equal, obtianing beta * e / omega = beta / k
# which yields k = omega / e
# this equals 0.3678794 * omega, which is a bit further off 



### Set scaling factors to approximat K for the photo-inhibition model
kf_d <- 0.23961
kf_ht <- 0.2020784
kf_mm <- 0.3678794



### Calculate half sat pelagic
pel_omega_mean <- pelagic_rds$data_list$par %>% mean() * 
  exp({rstan::extract(pelagic_rds$fit, pars = "lo_m")}$lo_m)
pel_half_omega_d <- kf_d * mean(pel_omega_mean)
pel_half_omega_ht <- kf_ht * mean(pel_omega_mean)
pel_half_omega_mm <- kf_mm * mean(pel_omega_mean)



### Calculate half sat benthic
ben_omega_mean <- benthic_rds$data_list$par %>% mean() * 
  exp({rstan::extract(benthic_rds$fit, pars = "lo_m")}$lo_m)
ben_half_omega_d <- kf_d * mean(ben_omega_mean)
ben_half_omega_ht <- kf_ht * mean(ben_omega_mean)
ben_half_omega_mm <- kf_mm * mean(ben_omega_mean)


#=========================================================================================





#=========================================================================================
#========== Parameters with temp
#=========================================================================================

##### Extract parameters

# pelagic
pel_pars <- pelagic_rds$fit_summary %>%
  filter(var %in% c(paste0("b[",1:11,"]"),
                    paste0("o[",1:11,"]"),
                    paste0("r[",1:11,"]"))) %>%
  mutate(name = strsplit(var, "\\[|\\]|,") %>% map_chr(~as.character(.x[1])),
         id = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2]))) %>%
  rename(lo = `16%`,
         mi = `50%`,
         hi = `84%`) %>%
  select(name, id, lo, mi, hi) %>%
  left_join(pelagic_rds$dd %>%
              select(date, site, event, temp) %>%
              unique() %>%
              mutate(id = row_number())) %>%
  select(-id) %>% 
  unique() %>%
  left_join(date_levels) %>%
  left_join(site_levels) %>%
  mutate(id = row_number(),
         id = as.numeric(date_level),
         id = id + ifelse(site_level == "north", -0.1, 
                          ifelse(site_level == "south", 0.1, 0)),
         year = year(as_date(date)))  %>%
  ungroup()


# benthic
ben_pars <-  benthic_rds$fit_summary %>%
  filter(var %in% c(paste0("b[",1:11,"]"),
                    paste0("o[",1:11,"]"),
                    paste0("r[",1:11,"]"))) %>%
  mutate(name = strsplit(var, "\\[|\\]|,") %>% map_chr(~as.character(.x[1])),
         id = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2]))) %>%
  rename(lo = `16%`,
         mi = `50%`,
         hi = `84%`) %>%
  select(name, id, lo, mi, hi) %>%
  left_join(benthic_rds$dd %>%
              select(date, site, event, temp) %>%
              unique() %>%
              mutate(id = row_number())) %>%
  select(-id) %>% 
  unique() %>%
  left_join(date_levels) %>%
  left_join(site_levels) %>%
  group_by(name) %>%
  arrange(site_level, event) %>%
  mutate(id = row_number(),
         id = as.numeric(date_level),
         id = id + ifelse(site_level == "north", -0.1, 
                          ifelse(site_level == "south", 0.1, 0.3)),
         year = year(as_date(date)))  %>%
  ungroup()

ben_pars %>%
  ggplot(aes(temp,
             mi))+
  facet_wrap(~name, scales = "free_y", nrow = 3)+
  geom_point(size = 2)+
  scale_y_continuous(trans = "log")+
  geom_smooth(method = "lm",formula = 'y ~ x', se = F)

pel_pars %>%
  ggplot(aes(temp,
             mi))+
  facet_wrap(~name, scales = "free_y", nrow = 3)+
  geom_point(size = 2)+
  scale_y_continuous(trans = "log")+
  geom_smooth(method = "lm",formula = 'y ~ x', se = F)

#=========================================================================================
