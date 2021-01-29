#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(lemon)
library(cowplot)

# import data
pelagic_rds <- read_rds("analysis/output/pelagic.rds")
benthic_rds <- read_rds("analysis/output/benthic.rds")

# set theme
theme_set(theme_bw() %+replace%
            theme(panel.grid = element_blank(),
                  strip.background = element_blank(),
                  plot.margin = margin(1,1,1,1),
                  legend.margin = margin(0,0,0,-4),
                  legend.text = element_text(size = 8),
                  axis.text = element_text(size = 10, color="black", family = "sans"),
                  axis.title = element_text(size =10),
                  axis.title.y = element_text(angle = 90, margin=margin(0,10,0,0)),
                  axis.title.x = element_text(margin = margin(10,0,0,0)),
                  panel.spacing = unit(0.1, "lines"),
                  axis.ticks = element_line(size = 0.25)))








#==========
#========== PI curves
#==========

pelagic_sum <- pelagic_rds$fit_summary %>%
  filter(str_detect(.$var, "nep_sum")) %>%
  mutate(name = strsplit(var, "\\[|\\]|,") %>% map_chr(~.x[1]),
         id = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2]))) %>%
  filter(name %in% c("nep_sum")) %>%
  rename(lo = `16%`,
         mi = `50%`,
         hi = `84%`) %>%
  select(id, lo, mi, hi) %>%
  unique() %>%
  left_join(pelagic_rds$data_list$dd_sum %>%
              mutate(id = row_number())) %>%
  select(-id) %>%
  unique()









dd %>%
  ggplot(aes(par, do_flux, color = date))+
  facet_wrap(~site, nrow = 2)+
  geom_ribbon(data = nep_sum,
              aes(x = par,
                  ymin = lo,
                  ymax = hi,
                  fill = date),
              linetype = 0,
              alpha = 0.2,
              inherit.aes = F)+
  geom_point()+
  geom_line(data = nep_sum,
            aes(y = mi))+
  theme_bw()+
  scale_color_viridis_d()+
  scale_fill_viridis_d()


pars <-  fit_summary %>%
  filter(var %in% c(paste0("b[",1:11,"]"),
                    paste0("a[",1:11,"]"),
                    # paste0("o[",1:11,"]"),
                    paste0("r[",1:11,"]"))) %>%
  mutate(name = strsplit(var, "\\[|\\]|,") %>% map_chr(~as.character(.x[1])),
         id = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2]))) %>%
  rename(lo = `16%`,
         mi = `50%`,
         hi = `84%`) %>%
  select(name, id, lo, mi, hi) %>%
  left_join(dd %>%
              select(date, site, event) %>%
              unique() %>%
              mutate(id = row_number())) %>%
  select(-id) %>%
  unique() 


pars %>%
  mutate(dummya = as.numeric(as.factor(site))) %>%
  group_by(site) %>%
  mutate(dummy = as.numeric(date) - mean(as.numeric(date))) %>%
  ungroup() %>%
  mutate(dummy = dummy / abs(max(dummy)) / 3,
         dummyb = dummya + dummy) %>%
  ggplot(aes(dummyb, mi, color = site))+
  facet_wrap(~name, scales = "free_y", nrow = 2)+
  geom_point(size = 2)+
  geom_errorbar(aes(ymin = lo,
                    ymax = hi),
                width = 0)+
  theme_bw()+
  scale_color_manual(values = c("dodgerblue",
                                "firebrick",
                                "black"))


