# Weekly projections of NZ electricity demand for 2020

# Data source: https://www.emi.ea.govt.nz/r/sepyr

# *****************************************************************************
# Setup ----

library(conflicted)
library(tidyverse)
library(lubridate)
library(my.r.functions)
library(scales)
library(here)
library(janitor)
library(sandwich)
library(lmtest)
library(broom)
library(forecast)
library(tsibble)
library(fable)

conflict_prefer("filter", "dplyr")
conflict_prefer("lag", "dplyr")

# *****************************************************************************


# *****************************************************************************
# Load data ----

# Daily electricity demand
dat_electricity <- read_csv(file = here("data/electricity-demand.csv"), 
                            skip = 10) %>%
  clean_names() %>%
  mutate(date = dmy(period_start)) %>%
  select(-period_start, -period_end, -trading_period, -region_id) 

# National public holidays
dat_public_hols <- read_csv(file = here("data/public-holidays.csv")) %>%
  arrange(holiday, date) %>%
  arrange(date)

# *****************************************************************************


# ***************************************************************************** 
# Weekly electricity demand model ----

# Generate week-of-year numbers from 1 to 52, with 8-day weeks in week 52 and
# in week 9 if the year is a leap year
generate_week_numbers <- function(d, y) {
  if (leap_year(y)) {
    weeks_in_year <- c(
      rep(1:8, each = 7), 
      rep(9, 8), 
      rep(10:51, each = 7), 
      rep(52, 8)
    )[1:nrow(d)]
  } else {
    weeks_in_year <- c(
      rep(1:51, each = 7), 
      rep(52, 8)
    )[1:nrow(d)]  
  }
  
  return(as.integer(weeks_in_year))
}

# Create weekly dataset from daily electricity data
dat_weekly <- dat_electricity %>%
  select(-region) %>%
  mutate(year = as.integer(year(date))) %>%
  nest(data = c(date, demand_g_wh)) %>%
  mutate(week_in_year = map2(.x = data, .y = year, .f = generate_week_numbers)) %>%
  unnest(c(data, week_in_year)) %>%
  left_join(y = dat_public_hols, by = "date") %>%
  mutate(holiday = ifelse(!is.na(holiday), 1L, 0L)) %>%
  group_by(year, week_in_year) %>%
  summarise(gwh = sum(demand_g_wh), 
            week_start_date = min(date), 
            holiday_days = sum(holiday), 
            n_days_in_week = n()) %>%
  ungroup() %>%
  filter(n_days_in_week >= 7) %>%
  mutate(t = row_number(), 
         month = as.integer(month(week_start_date)), 
         tb = ifelse(year >= 2016, 1L, 0L), 
         leap_week = ifelse((week_in_year == 9) & leap_year(year), 1, 0), 
         week_start_day = wday(x = week_start_date, week_start = 1)) %>%
  mutate(week_in_year = factor(week_in_year), 
         month = factor(month), 
         tb = factor(tb), 
         leap_week = factor(leap_week), 
         week_start_day = factor(week_start_day)) %>%
  as_tsibble(index = t)

# Plot weekly electricity demand pre March 2020
dat_weekly_chart <- dat_weekly %>%
  filter(week_start_date < ymd("2020-03-01")) %>%
  ggplot(mapping = aes(x = week_start_date, 
                       y = gwh)) + 
  geom_line(size = 0.25) + 
  geom_smooth(se = FALSE, method = "gam", size = 0.5) + 
  scale_x_date(breaks = date_breaks("1 year"), 
               labels = date_format("%Y"), 
               expand = expansion(0.02, 0)) + 
  scale_y_continuous(breaks = seq(0, 1000, 50))

output_chart(chart = dat_weekly_chart, 
             orientation = "wide", 
             path = here("outputs"), 
             xlab = "Week starting date", 
             ylab = "", 
             ggtitle = "Weekly electricity demand (GWh)", 
             plot.margin = margin(4, 4, 4, 4, "pt"), 
             axis.title.y = element_blank())


# Estimate weekly ARMA model using pre-March 2020 data
m_weekly <- dat_weekly %>%
  filter(week_start_date < ymd("2020-03-01")) %>%
  model(ar = ARIMA(
    formula = gwh ~ 
      tb * t + 
      week_in_year + 
      holiday_days + 
      leap_week
  ))

report(m_weekly)
accuracy(m_weekly)

# Plot fitted vs actual
m_weekly_fit_chart <- m_weekly %>%
  augment() %>%
  as_tibble() %>%
  left_join(y = dat_weekly %>% select(week_start_date, t), 
            by = "t") %>%
  select(week_start_date, gwh, .fitted) %>%
  mutate(pe = gwh / .fitted - 1) %>%
  mutate(sign = ifelse(pe > 0, "pos", "neg")) %>%
  ggplot(mapping = aes(x = week_start_date, 
                       y = pe, 
                       fill = sign)) + 
  geom_col(width = 5, 
           size = 0) + 
  geom_hline(yintercept = 0, size = 0.25, colour = grey(0.25)) + 
  scale_x_date(breaks = date_breaks("1 year"), 
               labels = date_format("%Y"), 
               expand = expansion(0.02, 0)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), 
                     breaks = seq(-0.07, 0.07, 0.01), 
                     limits = c(-0.07, 0.07)) + 
  scale_fill_brewer(palette = "Set2", guide = "none")

output_chart(chart = m_weekly_fit_chart, 
             path = here("outputs"), 
             orientation = "wide", 
             xlab = "Week starting date", 
             ylab = "", 
             ggtitle = "Actual vs predicted weekly electricity demand (%)", 
             plot.margin = margin(4, 4, 4, 4, "pt"), 
             axis.title.y = element_blank())

# Predictions for 2020 with confidence intervals
p2020 <- m_weekly %>%
  fabletools::forecast(new_data = dat_weekly %>% 
                         filter(year(week_start_date) == 2020) %>%
                         select(-gwh))

p2020_ci <- rbind(
  p2020 %>% mutate(interval = hilo(x = .distribution, level = 60)) %>% as_tibble(), 
  p2020 %>% mutate(interval = hilo(x = .distribution, level = 80)) %>% as_tibble(), 
  p2020 %>% mutate(interval = hilo(x = .distribution, level = 95)) %>% as_tibble()
) %>%
  unnest(interval) %>%
  select(t, week_start_date, predicted = gwh, .lower, .upper, .level) %>%
  left_join(y = dat_weekly %>% select(t, actual = gwh), 
            by = "t") %>%
  mutate(ap = actual/predicted - 1, 
         ap.lower = actual/.lower - 1, 
         ap.upper = actual/.upper - 1)

predictions_chart <- p2020_ci %>%
  ggplot(mapping = aes(x = week_start_date, 
                       y = ap, 
                       ymin = ap.lower, 
                       ymax = ap.upper, 
                       colour = fct_rev(factor(.level)))) + 
  geom_vline(xintercept = ymd("2020-03-26"), colour = "red") + 
  geom_vline(xintercept = ymd("2020-04-28"), colour = "darkorange") + 
  geom_vline(xintercept = ymd("2020-05-14"), colour = "darkgoldenrod1") + 
  geom_hline(yintercept = 0, colour = "black") + 
  geom_linerange(size = 6, data = p2020_ci %>% filter(.level == 95)) + 
  geom_linerange(size = 6, data = p2020_ci %>% filter(.level == 80)) + 
  geom_linerange(size = 6, data = p2020_ci %>% filter(.level == 60)) + 
  annotate(geom = "text", 
           x = ymd("2020-03-27"), 
           y = 0.11, 
           label = "Level 4", 
           colour = "red", 
           hjust = 0, 
           family = "Fira Sans", 
           fontface = "bold", 
           size = 3) +
  annotate(geom = "text", 
           x = ymd("2020-04-29"), 
           y = 0.11, 
           label = "Level 3", 
           colour = "darkorange", 
           hjust = 0, 
           family = "Fira Sans", 
           fontface = "bold", 
           size = 3) +
  annotate(geom = "text", 
           x = ymd("2020-05-15"), 
           y = 0.11, 
           label = "Level 2", 
           colour = "darkgoldenrod1", 
           hjust = 0, 
           family = "Fira Sans", 
           fontface = "bold", 
           size = 3) +
  scale_colour_manual(values = c("60" = "#3182bd", 
                                 "80" = "#9ecae1", 
                                 "95" = "#deebf7"), 
                      labels = c("60%", "80%", "95%"), 
                      name = "Confidence level") + 
  scale_y_continuous(labels = percent_format(accuracy = 1), 
                     breaks = seq(-0.20, 0.12, 0.02), 
                     limits = c(-0.20, 0.12), 
                     expand = expansion(0, 0)) + 
  scale_x_date(breaks = p2020_ci$week_start_date, 
               limits = c(ymd("2020-01-01", NA)), 
               labels = function(x) {
                 lx <- str_wrap(string = format(x = x, format = "%d %b"), 
                                width = 3)
                 return(lx)
               })

output_chart(chart = predictions_chart, 
             path = here("outputs"), 
             orientation = "wide", 
             xlab = "Week start date", 
             ylab = "", 
             ggtitle = "Actual vs predicted weekly electricity demand in 2020 (%)", 
             legend_position = "top", 
             plot.margin = margin(4, 4, 4, 4, "pt"), 
             axis.title.y = element_blank())

# *****************************************************************************
