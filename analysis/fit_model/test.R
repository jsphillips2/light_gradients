# set theme
theme_set(theme_bw() %+replace% 
            theme(panel.grid = element_blank(),
                  strip.background = element_blank(),
                  legend.margin = margin(0,0,0,0),
                  strip.text = element_text(size=10),
                  legend.text = element_text(size=10),
                  axis.text=element_text(size=10, color="black"),
                  axis.title.y=element_text(angle = 90 ,margin=margin(0,15,0,0)),
                  axis.title.x=element_text(margin=margin(15,0,0,0))))

pi_fun <- function( b_f,  a_f,  d_f,  l_f) {
  theta1 = b_f*l_f
  theta2 = (b_f*l_f*l_f)/(a_f*d_f*d_f)
  theta3 = (1 - 2*b_f/(a_f*d_f))*l_f
  theta4 = b_f/a_f
  pred = theta1/(theta2 + theta3 + theta4)
  return(pred)
}


test3 <- fit_clean %>%
  filter(name %in% c("beta0", "alpha", "rho", "delta","omega")) %>%
  select(site, sampledate, event_name, name, middle) %>%
  spread(name, middle) %>%
  mutate(g_b = {fit_clean %>% filter(name=="g_b")}$middle,
         g_r = {fit_clean %>% filter(name=="g_r")}$middle) %>%
  inner_join(prep_data %>%
               group_by(event_name) %>%
               expand(par = seq(min(par), max(par), length.out = 100),
                      temp = mean(temp))) %>%
  mutate(omega = 0,
         beta = beta0*g_b^(temp - 12),
         gpp = pi_fun(beta, alpha, delta, par),
         er = (rho + omega*par)*g_r^(temp - 12),
         nep = gpp - er)

# test3 <- fit_clean %>%
#   filter(name %in% c("beta0", "alpha", "rho", "omega")) %>%
#   select(site, sampledate, event_name, name, middle) %>%
#   spread(name, middle) %>%
#   mutate(g_b = {fit_clean %>% filter(name=="g_b")}$middle,
#          g_r = {fit_clean %>% filter(name=="g_r")}$middle) %>%
#   inner_join(prep_data %>%
#                group_by(event_name) %>%
#                expand(par = seq(min(par), max(par), length.out = 100),
#                       temp = mean(temp))) %>%
#   mutate(omega = 0,
#          beta = beta0*g_b^(temp - 12),
#          gpp = beta*tanh((alpha/beta)*par),
#          er = (rho + omega*par)*g_r^(temp - 12),
#          nep = gpp - er)



test3 %>%
  ggplot(aes(par, do_flux, color = event_name))+
  facet_wrap(~site, nrow = 2)+
  geom_line(aes(y = nep), size = 1)+
  geom_hline(yintercept = 0, size = 0.5, alpha = 0.5)+
  geom_point(data = prep_data %>%
               mutate(site = site_name),
             aes(y = do_flux), size=2.5, alpha = 0.6)+
  scale_y_continuous("Net Ecosystem Production")


fit_clean %>%
  filter(name %in% c("beta0", "alpha", "rho", "delta","omega")) %>%
  ggplot(aes(event_name, middle, color =site))+
  facet_wrap(~name, scales = "free_y", nrow = 3)+
  geom_point(size = 3.5)+
  geom_errorbar(aes(ymin = lower16, ymax = upper84), width = 0)+
  scale_color_manual(values=c("gray50","dodgerblue","firebrick"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))










test3 <- fit_clean %>%
  filter(name %in% c("beta0", "alpha", "rho")) %>%
  select(site, sampledate, event_name, name, middle) %>%
  spread(name, middle) %>%
  mutate(g_b = {fit_clean %>% filter(name=="g_b")}$middle,
         g_r = {fit_clean %>% filter(name=="g_r")}$middle) %>%
  inner_join(prep_data %>%
               group_by(event_name) %>%
               expand(par = seq(min(par), max(par), length.out = 100),
                      temp = mean(temp))) %>%
  mutate(beta = beta0*g_b^(temp - 12),
         gpp = beta*tanh((alpha/beta)*par),
         er = rho*g_r^(temp - 12),
         nep = gpp - er)

test3 %>%
  ggplot(aes(par, do_flux, color = event_name))+
  facet_wrap(~site, nrow = 2)+
  geom_line(aes(y = nep), size = 1)+
  geom_hline(yintercept = 0, size = 0.5, alpha = 0.5)+
  geom_point(data = prep_data %>%
               mutate(site = site_name),
             aes(y = do_flux), size=2.5, alpha = 0.6)+
  scale_y_continuous("Net Ecosystem Production")


fit_clean %>%
  filter(name %in% c("beta0", "alpha", "rho", "delta")) %>%
  ggplot(aes(event_name, middle, color =site))+
  facet_wrap(~name, scales = "free_y", nrow = 3)+
  geom_point(size = 3.5)+
  geom_errorbar(aes(ymin = lower16, ymax = upper84), width = 0)+
  scale_color_manual(values=c("gray50","dodgerblue","firebrick"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

