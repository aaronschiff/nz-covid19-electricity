# Analysis of electricity demand

# Data source: Electricity Authority

# *****************************************************************************
# Setup ----

library(conflicted)
library(tidyverse)
library(lubridate)
library(here)
library(my.r.functions)
library(janitor)
library(scales)

conflict_prefer("filter", "dplyr")

# *****************************************************************************


# *****************************************************************************
# Load data ----

elec_dat <- read_csv(file = here("data/Grid_demand_trends_20200422194646.csv"), 
                     skip = 11) %>%
  clean_names() %>%
  mutate(period_start = dmy(period_start)) %>%
  select(date = period_start, 
         region, 
         demand = demand_g_wh) %>%
  mutate(weekday = wday(date, week_start = 1), 
         month = month(date), 
         year = year(date), 
         week = week(date))

akl_temp_dat <- read_csv(file = here("data/akl-temp.csv"), 
                         skip = 8) %>%
  clean_names() %>%
  separate(col = date_nzst, into = c("date", "junk"), sep = ":") %>%
  select(-junk) %>%
  mutate(date = ymd(date))

akl_rain_dat <- read_csv(file = here("data/akl-rain.csv"), 
                         skip = 8) %>%
  clean_names() %>%
  separate(col = date_nzst, into = c("date", "junk"), sep = ":") %>%
  select(-junk) %>%
  mutate(date = ymd(date))

akl_sunshine_dat <- read_csv(file = here("data/akl-sunshine.csv"), 
                             skip = 8) %>%
  clean_names() %>%
  separate(col = date_nzst, into = c("date", "junk"), sep = ":") %>%
  select(-junk) %>%
  mutate(date = ymd(date))

public_hols <- read_csv(file = here("data/public-holidays.csv")) %>%
  arrange(holiday, date)

# *****************************************************************************

test <- elec_dat %>% filter(region == "Upper North Island") %>%
  left_join(y = akl_temp_dat %>%
              select(date, tmax_c, tmin_c), 
            by = "date") %>%
  left_join(y = akl_rain_dat %>%
              select(date, amount_mm), 
            by = "date") %>%
  left_join(y = akl_sunshine_dat %>%
              select(date, amount_hrs), 
            by = "date") %>%
  left_join(y = public_hols, by = "date") %>%
  mutate(holiday = replace_na(holiday, "not holiday")) %>%
  mutate(xmas_break = case_when(
    (month(date) == 12) & (day(date) > 24) ~ 1, 
    (month(date) == 1) & (day(date) < 4) ~ 1, 
    TRUE ~ 0
  )) %>%
  mutate(t = row_number(), 
         e_tmax_c = exp(tmax_c / 100), 
         e_tmin_c = exp(tmin_c / 100))

test_m <- lm(demand ~ 
               factor(weekday) + 
               t + 
               tmax_c * factor(week) + 
               tmin_c * factor(week) + 
               amount_hrs * factor(week) +
               amount_mm + 
               relevel(x = as.factor(holiday), ref = "not holiday"), 
             data = test %>% filter(year < 2020))
summary(test_m)

test %>% 
  filter(year > 2017) %>%
  mutate(p = predict(test_m, newdata = test %>% filter(year > 2017))) %>%
  ggplot(mapping = aes(x = date)) + 
  geom_line(mapping = aes(y = demand)) + 
  geom_line(mapping = aes(y = p), colour = "red")

test %>%
  filter(year > 2019) %>%
  mutate(p = predict(test_m, newdata = test %>% filter(year > 2019))) %>%
  mutate(delta = demand / p - 1) %>%
  ggplot(mapping = aes(x = date, y = delta)) + 
  geom_col() + 
  scale_y_continuous(labels = percent_format(accuracy = 1.0))

weekly_test <- test %>%
  filter(date > ymd("2019-12-31")) %>%
  mutate(p = predict(test_m, newdata = test %>% filter(date > ymd("2019-12-31")))) %>%
  filter(!is.na(p)) %>%
  group_by(year, week) %>%
  summarise(demand = sum(demand), 
            p = sum(p)) %>%
  ungroup() %>%
  mutate(delta = demand / p - 1, 
         d = row_number()) %>%
  mutate(sign = ifelse(delta > 0, "pos", "neg")) %>%
  ggplot(mapping = aes(x = d, y = delta, fill = sign)) + 
  geom_col() + 
  scale_y_continuous(labels = percent_format(accuracy = 1.0), 
                     breaks = seq(-0.22, 0.1, 0.02)) + 
  scale_x_continuous(breaks = 1:20) + 
  scale_colour_manual(values = c("pos" = "cornflowerblue", 
                                 "neg" = "darkgoldenrod3"), 
                      aesthetics = c("colour", "fill"), 
                      guide = "none")

output_chart(chart = weekly_test, 
             path = here("outputs"), 
             orientation = "wide", 
             xlab = "Week of year", 
             ylab = "Percent")
