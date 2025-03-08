catching <- read.csv("Data/catch_data.csv")
sighting<- read.csv("Data/sight_data.csv")

sighting$Date <- dmy(sighting$Date)
catching$Date <- dmy(catching$Date)

# sighting2 <- sighting %>%
#   dplyr::select(Name, Date, Field.Season.Number, Time, Lat, Long, Sex, Location) %>%
#   filter(Name != "", Field.Season.Number != "", Lat != "", Long != "") %>%
#   group_by(Name) %>% 
#   mutate(First.Contact = min(Date), Month.Year = paste0(month(Date), ".", year(Date))) %>% ungroup()
# 
# catching2 <- catching %>%
#   dplyr::select(Name, Date, Field.Season.Number, Diseased.Visual, PCR.Diseased) %>%
#   filter(Name != "", Diseased.Visual != "", !is.na(Date)) %>% mutate(Month.Year = paste0(month(Date), ".", year(Date))) 

sighting2 <- sighting %>%
  dplyr::select(Name, Date, Field.Season.Number, Time, Lat, Long, Sex, Location) %>%
  filter(Name != "", Field.Season.Number != "", Lat != "", Long != "", !is.na(Date)) %>%
  group_by(Name) %>%
  mutate(
    First.Contact.Date = min(Date),
    Month.Year = ifelse(
      month(Date) <= 6, # If month is Jan-Jun
      paste0("1.", year(Date)),
      paste0("2.", year(Date)) # If month is Jul-Dec
    )
  ) %>%
  ungroup()

catching2 <- catching %>%
  dplyr::select(Name, Date, Field.Season.Number, Diseased.Visual, PCR.Diseased, Age) %>%
  filter(Name != "", Diseased.Visual != "", !is.na(Date)) %>%
  mutate(
    Month.Year = ifelse(
      month(Date) <= 6, # If month is Jan-Jun
      paste0("1.", year(Date)),
      paste0("2.", year(Date)) # If month is Jul-Dec
    )
  )

catching2$Disease.Status <- NA

catching2$Disease.Status[catching2$Diseased.Visual == "Yes" & catching2$PCR.Diseased == ""] <- 1
catching2$Disease.Status[catching2$Diseased.Visual == "No" & catching2$PCR.Diseased == ""] <- 0
catching2$Disease.Status[catching2$Diseased.Visual == "Yes" | catching2$Diseased.Visual == "Unknown" & catching2$PCR.Diseased == "Yes"] <- 1
catching2$Disease.Status[catching2$Diseased.Visual == "No" | catching2$Diseased.Visual == "Unknown" & catching2$PCR.Diseased == "No"] <- 0
catching2$Disease.Status[catching2$Diseased.Visual == "No" & catching2$PCR.Diseased == "Yes"] <- 1

duplicated_rows <- catching2 %>%
  group_by(Name, Month.Year) %>%
  filter(n() > 1) %>% # Keep groups with more than one row
  ungroup()

catching2 <- catching2 %>%
  arrange(Name, Field.Season.Number, Date) %>%
  group_by(Name, Month.Year) %>%
  filter( #remove dupliate catches
    if (n() > 1) {
      row_number() == 1 
    } else {
      row_number() == 1 
    }
  ) %>%
  ungroup()

catching3 <- catching2 %>% filter(!is.na(Disease.Status)) %>%
  group_by(Name) %>% arrange(Date) %>%
  mutate(
    Date.First.Diseased = if (any(Disease.Status == 1)){
      min(Date[Disease.Status == 1])}
    else {
      NA
    }
  ) %>% ungroup() %>%
  dplyr::select(-Diseased.Visual, -PCR.Diseased) %>% distinct()

# 1. Sort catches2 by Date within each Name and Field.Season.Number
catching3 <- catching3 %>%
  arrange(Name, Field.Season.Number, Date) 

disease.date <- catching3 %>%
  dplyr::select(Name, Date.First.Diseased, Month.Year, Disease.Status) 

duplicated_rows <- disease.date %>%
  group_by(Name, Month.Year) %>%
  filter(n() > 1) %>% # Keep groups with more than one row
  ungroup()

merge <- sighting2 %>%
  left_join(disease.date, by = c("Name", "Month.Year"))

final_df <- data.frame()

for (name in unique(merge$Name)) {
  test_df <- merge %>% filter(Name == name)
  
  first_disease_date <- test_df$Date.First.Diseased[!is.na(test_df$Date.First.Diseased)][1]
  
  for (i in 1:nrow(test_df)) {
    if (is.na(test_df$Date.First.Diseased[i])) {
      test_df$Date.First.Diseased[i] <- first_disease_date
    }
  }
  final_df <- bind_rows(final_df, test_df)
}

final_df <- final_df %>%
  rename(Season = Field.Season.Number) %>%
  arrange(Name, Date) %>%
  group_by(Name) %>%
  mutate(
    Disease.Status = {
      status_values <- Disease.Status
      last_non_na <- NA
      
      for (i in 1:n()) {
        if (!is.na(status_values[i])) {
          last_non_na <- status_values[i]
        } else {
          if (!is.na(last_non_na)) {
            status_values[i] <- last_non_na
          }
        }
      }
      status_values
    }
  ) %>%
  ungroup()

final_df <- final_df %>%
  group_by(Name) %>%
  mutate(
    Disease.Status = {
      status_values <- Disease.Status
      
      for (i in 1:n()) {
        if (is.na(Disease.Status[i]) && !is.na(Date.First.Diseased[i])) {
          if (Date[i] < Date.First.Diseased[i]) {
            status_values <- 0
          } else if (Date[i] >= Date.First.Diseased[i]) {
            status_values <- 1
          }
        } else if (is.na(Disease.Status[i] && is.na(Date.First.Diseased[i]))) {
          status_values <- 0
        } else if (!is.na(Disease.Status[i])) {
          status_values <- Disease.Status[i]
        }
      }
      status_values
    }
  ) %>% ungroup()

for (i in 1:nrow(final_df)) {
  if (is.na(final_df$Disease.Status[i]) && !is.na(final_df$Date.First.Diseased[i])) {
    if (final_df$Date[i] >= final_df$Date.First.Diseased[i]) {
      final_df$Disease.Status[i] = 1
    } else if (final_df$Date[i] < final_df$Date.First.Diseased[i]) {
      final_df$Disease.Status[i] = 0
    }
  } else if (is.na(final_df$Disease.Status[i]) && is.na(final_df$Date.First.Diseased[i])) {
    final_df$Disease.Status[i] = 0
  }
}

write.csv(final_df, "Data/Dragon_Data.csv", row.names = F)

errors <- final_df %>%
  group_by(Name) %>%
  filter(Date >= Date.First.Diseased & is.na(Disease.Status)) %>% ungroup()

errors <- final_df %>%
  group_by(Name) %>%
  filter(is.na(Disease.Status)) %>% ungroup()

# final_df <- final_df %>%
#   arrange(Name, Date) %>%
#   group_by(Name) %>%
#   mutate(
#     Disease.Status = {
#       status_values <- Disease.Status
#       last_non_na <- NA
#       
#       for (i in 1:n()) {
#         if (!is.na(Date.First.Diseased[i])) {
#           if (Date[i] < Date.First.Diseased[i]) {
#             status_values[i] <- 0 # Set to 0 if Date < Date.First.Diseased
#           } else if (Date[i] > Date.First.Diseased[i] && is.na(status_values[i])) {
#             if (i > 1) {
#               previous_non_na <- status_values[1:(i - 1)][!is.na(status_values[1:(i - 1)])]
#               if (length(previous_non_na) > 0) {
#                 last_previous_non_na <- tail(previous_non_na, 1)
#                 if (last_previous_non_na == 0) {
#                   status_values[i] <- 0
#                 } else {
#                   status_values[i] <- 1
#                 }
#               } else {
#                 status_values[i] <- 1
#               }
#             } else {
#               status_values[i] <- 1
#             }
#           }
#         }
#         
#         if (!is.na(status_values[i])) {
#           last_non_na <- status_values[i]
#         } else {
#           if (!is.na(last_non_na)) {
#             status_values[i] <- last_non_na
#           }
#         }
#       }
#       status_values
#     }
#   ) %>%
#   ungroup()




