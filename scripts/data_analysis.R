# libs ----

pacman::p_load(
    tidyverse,
    ggplot2
)

# Set path to script location ----

# https://stackoverflow.com/questions/47044068/get-the-path-of-current-script
# get path of source file
getCurrentFileLocation <-  function()
{
    this_file <- commandArgs() %>% 
        tibble::enframe(name = NULL) %>%
        tidyr::separate(col=value, into=c("key", "value"), sep="=", fill='right') %>%
        dplyr::filter(key == "--file") %>%
        dplyr::pull(value)
    if (length(this_file)==0)
    {
        this_file <- rstudioapi::getSourceEditorContext()$path
    }
    return(dirname(this_file))
}
# sets path on source file location
script_path <- getCurrentFileLocation()
setwd(script_path)

# Animal metadata ----

metadata <-
    read_csv("../datasets/metadata/general_metadata.csv")

# Weight ----

weight_data <-
    read_csv("../datasets/weights/objective_1_weights.csv") %>% 
    left_join(
        .,
        metadata,
        by = c("animal"="ID")
    )

weight_data %>% 
    group_by(animal) %>% 
    summarise(
        start = min(date),
        end = max(date)
    )

# Intake ----

col_name <- c(
    "time",
    "animal",
    "pellets",
    "boops",
    "motorTurns",
    "battery",
    "delay",
    "protocol"
)

p0 <- readRDS("../../PHD_data/objective_1/data/raw/low_uncertainty.rds")
p1 <- readRDS("../../PHD_data/objective_1/data/raw/high_uncertainty.rds")
X <- readRDS("../../PHD_data/objective_1/analysis/intake_data.rds") %>% 
    mutate(
        epoch = as.integer(as.POSIXct(time)),
        animal = as.character(ID),
        pellets = as.numeric(pellets),
        boops = NA,
        motorTurns = NA,
        battery = NA,
        protocol = experimental_phase
    )
p2 <- X %>% filter(animal%in%c(413, 414, 415, 416, 417, 418, 419, 420, 421, 422))
p3 <- X %>% filter(animal%in%c(494, 495, 496, 497, 498, 499, 500, 501))

tmp <-
    p2 %>% 
    mutate(
        motorTurns = NA,
        boops = NA
    ) %>% 
    ungroup() %>% 
    select(
        epoch,
        animal,
        pellets,
        boops,
        motorTurns,
        battery,
        delay,
        protocol
    )

colnames(tmp) <- col_name

tmp %>% 
    group_by(animal) %>% 
    summarise(
        start = min(date),
        end = max(date)
    )

tmp <- tmp %>% 
    mutate(date = lubridate::date(as.POSIXct(time, origin="1970-01-01", tz="GMT")),
           time_human = as.POSIXct(time, origin="1970-01-01", tz="GMT")) %>% 
    filter(date > "2022-03-23") %>% 
    complete(animal, date) %>% 
    group_by(animal) %>% 
    arrange(date) %>% 
    ungroup() %>% 
    mutate(
        pp1 = zoo::na.locf(protocol))
tmp

tmp1 <-
    tmp %>% 
    mutate(date = lubridate::date(as.POSIXct(time, origin="1970-01-01", tz="GMT"))) %>% 
    group_by(animal, date) %>% 
    summarise(
        intake = n()
    ) %>% 
    ungroup() %>% 
    complete(animal, date) %>% 
    mutate(
        intake = replace_na(intake, 0),
        int_bin = as.factor(if_else(intake > 30, 1, 0))
    ) %>% 
    drop_na()

tmp1 %>% 
    ggplot(aes(
        date, animal, fill = as.factor(int_bin)
    )) +
    geom_tile()

fix <-
tmp1 %>% 
    filter() %>% 
    select(!intake) %>% 
    arrange((date))

tmp <- tmp %>% 
    left_join(
        ., fix, by = c("animal", "date")
    ) 


replacements <-
    fix %>% 
    filter(int_bin==0) %>% 
    mutate(r = row_number()) %>% 
    group_by(r) %>% 
    group_split() %>% 
    map_dfr(
        ., function(bad_data){
            d <- tmp %>% 
                filter(
                    animal == bad_data$animal,
                    date == bad_data$date
                )
            bad_data_protocol <-
                d %>% 
                pull(pp1) %>% 
                {.[1]}
            good_data <-
                tmp %>% 
                filter(
                    animal == bad_data$animal,
                    protocol == bad_data_protocol,
                    int_bin == 1
                )
            good_data_dates <-
                unique(good_data$date)
            sample_good_date <-
                sample(good_data_dates, size = 1, replace = TRUE)
            new_data <-
                good_data %>% 
                filter(date == sample_good_date)
            diff <- -1 * as.numeric(new_data$date[1] - d$date[1], unit="secs")
            replace_data <-
                new_data %>% 
                mutate(
                    time = time + diff,
                    date = as.POSIXct(time, origin="1970-01-01", tz="GMT")
                )
            return(replace_data)
        }
    )

fixed_data <-
    bind_rows(tmp, replacements) %>% 
    filter(int_bin == 1)

fixed_data %>% 
    mutate(date = lubridate::date(date)) %>% 
    group_by(animal, date, int_bin, protocol) %>% 
    summarise(
        intake = n()
    ) %>% 
    ggplot(aes(
        date, animal, fill = as.factor(protocol)
    )) +
    geom_tile(
        color = "white"
    ) +
    geom_text(aes(label = intake), color = "black", size = 3) +
    coord_fixed()




## Plots (A) ----

a1 <-
    weight_data %>% 
    filter(experiment_id=="uncertainty_A") %>% 
    ggplot(aes(
        date, weight, group=animal, color=group
    )) +
    geom_vline(xintercept = weight_data$uncertainty_start[1]) +
    geom_line() +
    geom_point()
a1

a2 <-
    weight_data %>% 
    filter(experiment_id=="uncertainty_A",
           date >= uncertainty_start) %>% 
    group_by(animal) %>% 
    mutate(
        start = head(weight,1),
        end = tail(weight, 1)
    ) %>% 
    pivot_longer(
        cols = c(start, end)
    ) %>% 
    ggplot(aes(
        interaction(group, name), value, group=animal
    )) +
    geom_point()
a2

# Lickometer
    
