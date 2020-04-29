# Projections of NZ electricity demand for 2020

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

conflict_prefer("filter", "dplyr")
conflict_prefer("lag", "dplyr")

# *****************************************************************************


# *****************************************************************************
# Utility functions ----

# Parse horrible NIWA Cliflo data files
parse_niwa_data <- function(f) {
  # Load raw data
  d_raw <- read_csv(f, 
                    col_names = paste0("x", 1:12), 
                    col_types = "cccccccccccc", 
                    skip = 7)
  
  # Clean up raw data a bit
  d_clean <- d_raw %>%
    mutate(is_header_row = str_detect(x1, ": Daily") & is.na(x2)) %>%
    mutate(dataset = ifelse(is_header_row, x1, NA_character_)) %>%
    fill(dataset, .direction = "down") %>%
    select(-is_header_row)
  
  # Sunshine data
  d_sunshine <- d_clean %>% filter(dataset == "Sunshine: Daily") %>%
    filter(!is.na(x2)) %>%
    select(-x5:-dataset)
  names(d_sunshine) <- as.character(d_sunshine[1, ])
  d_sunshine <- d_sunshine %>%
    clean_names() %>%
    mutate(amount_hrs = as.numeric(amount_hrs), 
           period_hrs = as.numeric(period_hrs)) %>%
    filter(row_number() != 1) %>%
    separate(col = date_nzst, into = c("date", "junk"), sep = "\\:", convert = TRUE) %>%
    select(-junk) %>%
    mutate(date = ymd(date)) %>%
    mutate(dataset = "sunshine")
  
  # Rain data
  d_rain <- d_clean %>% filter(dataset == "Rain: Daily") %>%
    filter(!is.na(x2)) %>%
    select(-x9:-dataset)
  names(d_rain) <- as.character(d_rain[1, ])
  d_rain <- d_rain %>%
    clean_names() %>%
    mutate(amount_mm = as.numeric(amount_mm), 
           sof_g = as.numeric(sof_g), 
           deficit_mm = as.numeric(deficit_mm), 
           runoff_mm = as.numeric(runoff_mm), 
           period_hrs = as.integer(period_hrs)) %>%
    filter(row_number() != 1) %>%
    separate(col = date_nzst, into = c("date", "junk"), sep = "\\:", convert = TRUE) %>%
    select(-junk) %>%
    mutate(date = ymd(date)) %>%
    mutate(dataset = "rain")
  
  # Temperature data
  d_temp <- d_clean %>% filter(dataset == "Max_min: Daily") %>%
    filter(!is.na(x2)) %>%
    select(-dataset)
  names(d_temp) <- as.character(d_temp[1, ])
  d_temp <- d_temp %>%
    clean_names() %>%
    mutate(tmax_c = as.numeric(tmax_c), 
           period_hrs = as.integer(period_hrs), 
           tmin_c = as.numeric(tmin_c), 
           period_hrs_2 = as.integer(period_hrs_2), 
           tgmin_c = as.numeric(tgmin_c), 
           period_hrs_3 = as.integer(period_hrs_3), 
           tmean_c = as.numeric(tmean_c), 
           r_hmean_percent = as.numeric(r_hmean_percent), 
           period_hrs_4 = as.integer(period_hrs_4)) %>%
    filter(row_number() != 1) %>%
    separate(col = date_nzst, into = c("date", "junk"), sep = "\\:", convert = TRUE) %>%
    select(-junk) %>%
    mutate(date = ymd(date)) %>%
    mutate(dataset = "temp")
  
  # Combine
  w <- bind_rows(
    # Sunshine
    d_sunshine %>% 
      select(-period_hrs) %>%
      pivot_longer(cols = amount_hrs,
                   names_to = "measure", 
                   values_to = "value"), 
    
    # Rain
    d_rain %>% 
      select(station, date, amount_mm, dataset) %>%
      pivot_longer(cols = amount_mm, 
                   names_to = "measure", 
                   values_to = "value"),
    
    # Temperature
    d_temp %>%
      select(station, date, tmax_c, tmin_c, dataset) %>%
      pivot_longer(cols = c(tmax_c, tmin_c), 
                   names_to = "measure", 
                   values_to = "value")
  ) %>%
    mutate(measure = paste(dataset, measure, sep = ".")) %>%
    select(-dataset)
  
  return(w)
}

# *****************************************************************************


# *****************************************************************************
# Load data ----

# Electricity demand - base
dat_electricity_base <- read_csv(file = here("data/electricity-demand-base.csv"), 
                                 skip = 10) %>%
  clean_names() %>%
  mutate(date = dmy(period_start)) %>%
  select(-period_start, -period_end, -trading_period, -region_id) 

# Weather data - base
dat_weather_base <- bind_rows(
  parse_niwa_data(f = here("data/weather-base-auckland.csv")), 
  parse_niwa_data(f = here("data/weather-base-hamilton.csv")), 
  parse_niwa_data(f = here("data/weather-base-wellington.csv")), 
  parse_niwa_data(f = here("data/weather-base-christchurch.csv")), 
  parse_niwa_data(f = here("data/weather-base-dunedin.csv"))
) 

# Weather concordance with electricity regions
dat_weather_electricity_regions <- tribble(
  ~region, ~station, 
  "Upper North Island", "37852", 
  "Central North Island", "26117", 
  "Lower North Island", "3385", 
  "Lower North Island", "3445", 
  "Upper South Island", "4843", 
  "Lower South Island", "15752"
)

# Public holidays
# TODO: Get anniversary days for other provinces


# Combined daily data in wide format for modelling, nested by region
dat_daily <- bind_rows(
  # Electricity
  dat_electricity_base %>%
    mutate(measure = "demand_g_wh") %>%
    rename(value = demand_g_wh), 
  
  # Weather
  dat_weather_base %>%
    left_join(y = dat_weather_electricity_regions, by = "station") %>%
    select(-station)
) %>%
  mutate(weekday = wday(date, week_start = 1), 
         month = month(date), 
         year = year(date), 
         week = week(date)) %>%
  mutate(quarter = ceiling(month / 3)) %>%
  select(region, measure, date, year, month, quarter, week, weekday, value) %>%
  arrange(region, measure, date) %>%
  pivot_wider(names_from = measure, values_from = value) %>%
  group_by(region) %>%
  mutate(t = row_number()) %>%
  ungroup() %>%
  nest(data = -region)

# Quarterly real GDP
dat_real_gdp_actual <- read_csv(file = here("data/SNE445001_20200429_050839_57.csv"), 
                                skip = 1) %>%
  clean_names() %>%
  rename(yq = x1, gdp_actual = gross_domestic_product_production_measure) %>%
  filter(!is.na(gdp_actual)) %>%
  separate(col = yq, into = c("year", "quarter"), sep = "Q", convert = TRUE)

dat_real_gdp_seas_adj <- read_csv(file = here("data/SNE446901_20200430_084253_24.csv"), 
                                  skip = 1) %>%
  clean_names() %>%
  rename(yq = x1, gdp_seas_adj = gross_domestic_product_production_measure) %>%
  filter(!is.na(gdp_seas_adj)) %>%
  separate(col = yq, into = c("year", "quarter"), sep = "Q", convert = TRUE)

# *****************************************************************************


# *****************************************************************************
# Model daily electricity demand ----

daily_electricity_demand_model <- function(d) {
  m <- lm(
    demand_g_wh ~ 
      factor(weekday) + 
      t + 
      temp.tmax_c * factor(week) + 
      temp.tmin_c * factor(week) + 
      sunshine.amount_hrs * factor(week) +
      rain.amount_mm, 
      #relevel(x = as.factor(holiday), ref = "not holiday") 
    data = d %>% filter(year < 2020)
  )
  
  return(m)
}

models_daily <- dat_daily %>%
  mutate(m_daily = map(.x = data, .f = daily_electricity_demand_model))

# *****************************************************************************


# *****************************************************************************
# Model GDP vs electricity demand ----

dat_quarterly <- dat_electricity_base %>%
  mutate(year = year(date), quarter = ceiling(month(date) / 3)) %>%
  group_by(year, quarter) %>%
  summarise(gwh = sum(demand_g_wh)) %>%
  ungroup() %>%
  left_join(y = dat_real_gdp_actual, by = c("year", "quarter")) %>%
  filter(!is.na(gdp_actual)) %>%
  arrange(year, quarter) %>%
  mutate(t = row_number()) %>%
  mutate(log_gdp_actual = log(gdp_actual), 
         log_gwh = log(gwh)) %>%
  mutate(d_log_gdp_actual = log_gdp_actual - lag(log_gdp_actual, 4), 
         d_log_gwh = log_gwh - lag(log_gwh, 4)) %>%
  left_join(y = dat_real_gdp_seas_adj, by = c("year", "quarter")) %>%
  filter(!is.na(gdp_seas_adj)) %>%
  mutate(log_gdp_seas_adj = log(gdp_seas_adj)) %>%
  mutate(d_log_gdp_seas_adj = log_gdp_seas_adj - lag(log_gdp_seas_adj, 4))

m_actual <- lm(gdp_actual ~ gwh + t * factor(quarter) , 
               data = dat_quarterly)
summary(m_actual)
coeftest(x = m_actual, vcov = vcovHAC)

m_seas_adj <- lm(gdp_seas_adj ~ gwh + t * factor(quarter) , 
                 data = dat_quarterly)
summary(m_seas_adj)
coeftest(x = m_seas_adj, vcov = vcovHAC)



# ***************************************************************************** 