####################################
####                            ####
####     Annotate KIWA Data     ####
####      Scott Yanco, PhD      ####
####    scott.yanco@yale.edu    ####
####                            ####
####################################


####---- Description ----####

# This script executes environmental  annotations of Kirtland's Warbler (KIWA)
# for inclusion as individual-specific covariates in multi-state full annual
# survival models  Annotating data with winter (March) EVI and landcover.  
# Annotations will be spatially linked to the geometric centroid of winter fix 
# locations.  We want march EVI data so will add a 'synthetic' timestamp of 
# March 31 (and the associated year) and get data with temporal buffer of 1 
# month (so will capture the preceding month). Irrelevant for the landcover data
# because that ahs temporal buffer of 1 year. Annotations completed using the
# STOAT annotation tool from mol.org/stoat


####---- Initialization ----####

library(tidyverse)
library(rstoat)
library(lubridate)
library(glue)
library(sf)


####---- LOAD & PROCESS DATA ----####


##-- Load data --##
message("Loading data...")
load("data/KW_Survival_With_Age July 27 2021.rda")

message(glue("Data loaded with {nrow(KW_Summary)} rows..."))
#Assign states based on detection tower (e.g., MI towers = breeding_)
message("Assigning seasonal states...")
kiwa_dat <- KW_Summary %>% 
  
  # TODO: still excluding daily duplicates in finding winter centroid, warranted?
  # group_by(motusTagID) %>% #group by tag ID
  # filter(!duplicated(ts_Yday)) %>%  #remove duplicates within a dAy
  # ungroup()%>% 
  
  mutate(year = year(ts)) %>% #create year var
  mutate(ch_ms = case_when( #create seasonal state cap history
    Period == "Tagged" ~ 1,
    Period == "Winter" ~ 1,
    Period == "Migration" ~ 2,
    Period == "Breeding" ~ 3,
    Period == "Died" ~ 4,
    #observations after breeding will get captured as breeding state for the 
    #purposes of feeding the model known states - these obs get censored below
    Period == "Breeding_Departure" ~ 3, 
    Period == "Fall_Migration" ~ 3)) %>% 
  
  #filter to only include winter fixes
  filter(Period == "Winter" | Period == "Tagged")


##-- Get Centroids --##
message("Getting winter centroids...")


#convert to sf, grab centroids, add datestamp
kiwa_cent <- kiwa_dat %>% 
  # convert to `sf` object
  st_as_sf(coords = c("recvLon","recvLat"),  crs = 4326) %>% 
  st_transform(26718) %>% # reproject to UTm for centroid calc 
  group_by(motusTagID) %>%  # group by individual
  summarize(geometry = st_combine(geometry), # combine individual geoms
            year = year[1], # carry year variable over
            cent = st_centroid(geometry)) %>%  # get individual centroid
  mutate(ts = ymd(glue("{year}-03-31")))  #create march date stamp for annos


# convert object to be readable by STOAT
message("Converting to STOAT readable format...")
kiwa_cent_stoat <- st_sf(kiwa_cent$cent, motusTagID = kiwa_cent$motusTagID, 
                         ts = kiwa_cent$ts) %>% 
  st_transform(4326) %>% 
  mutate(x = st_coordinates(.)[,1], #Add coordinates columns for stoat
         y = st_coordinates(.)[,2])



##-- Annotate --##
message("Annotating data...")

# EVI
#TODO: STOAT currently cannot find modis products... would be better with 
#those than current landsat due to revist freq - correct once STOAT is fixed
evi_lyr <- 'landsat8-evi-250-30'

kiwa_evi <- start_annotation_simple(kiwa_cent_stoat, layers = evi_lyr, 
                                    coords = c("x", "y"), date = "ts") %>% 
  rename_with(.fn = ~c(glue("variable_{evi_lyr}"),
                       glue("s_buff_{evi_lyr}"),
                       glue("product_{evi_lyr}"),
                       glue("valid_pixel_count_{evi_lyr}"),
                       glue("value_{evi_lyr}"),
                       glue("stdev_{evi_lyr}"),
                       glue("t_buff_{evi_lyr}")),
              .cols = c(variable, s_buff, product, valid_pixel_count, value, stdev, t_buff)) %>% 
  as.data.frame()

# # Landcover
#
# NOTE:  Depricating landcover for now b/c ESA CCI not picking up on enough details
# 
# lc_lyr <- 'esacci-landcover-300-365'
# 
# kiwa_lc <- start_annotation_simple(kiwa_cent_stoat, layers = lc_lyr, 
#                                    coords = c("x", "y"), date = "ts") %>% 
#   rename_with(.fn = ~c(glue("variable_{lc_lyr}"),
#                        glue("s_buff_{lc_lyr}"),
#                        glue("product_{lc_lyr}"),
#                        glue("valid_pixel_count_{lc_lyr}"),
#                        glue("value_{lc_lyr}"),
#                        glue("t_buff_{lc_lyr}")),
#               .cols = c(variable, s_buff, product, valid_pixel_count, value, t_buff)) %>% 
#   as.data.frame()


##-- Join Annotations --##

# kiwa_anno <- kiwa_evi %>%
#   full_join(kiwa_lc, by = "motusTagID")

# TODO: Kill this line and uncomment above if landcover product is added
kiwa_anno <- kiwa_evi


####---- FINALIZE SCRIPT ----####
message("Saving results...")
write_csv(kiwa_anno, file = "output/kiwa_annos.csv")
