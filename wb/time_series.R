library(fpoiscpm)
library(data.table)
library(gridExtra)
library(readr)

#Function to get the sum of tmax nmet by isolate for each day in a ward between two dates
#Usefull only if you want to generate period data from database

time_series_for_ward <- function(ward, date1, date2){

  date1 <- as.Date(date1)
  date2 <- as.Date(date2)

  #initiate the table
  ward_time_series <- data.table(date = as.Date(character(0)), sum_max_value = numeric(0))

  #get the list of date between the two dates
  date_list <- seq(date1, date2, by = 'days')
  date_list <- as.character(date_list)

  for (datetoadd in date_list){

    #get the list of isolate for one day for one ward
    isolate_day_ward <- DB_get_isolate_day_ward(connec = con, date1 = datetoadd, uf_group = ward)

    #initiate the sum value
    sum_max_nmet <- 0

    for (isolate in isolate_day_ward$idisolate){

      #get the pairof nmet for the isolate
      pairnmet <- DB_get_pairofnmet_same_ward(connec = con , isolate = isolate)
      pairnmet <- pairnmet[!is.na(value)]

      if (length(pairnmet$value) > 0){

        #add the max of nmet to the sum value
        sum_max_nmet <- sum_max_nmet + max(pairnmet$value)

      }
    }

    #add the new row ro the table
    new_row <- data.table(date = as.Date(datetoadd), sum_max_value = sum_max_nmet)
    ward_time_series <- rbind(ward_time_series, new_row)
  }

  ward_time_series
}

#Loading data-------------------------------------------------------------------

nb_isolate_per_iduf_group <- read_csv("~/EpiTrack_GIT/epitrack-gui/data/nb_isolate_per_iduf_group.csv")
nb_nmet_per_iduf_group <- read_csv("~/EpiTrack_GIT/epitrack-gui/data/nb_nmet_per_iduf_group.csv")
cluster_ufgroup <- read_csv("~/EpiTrack_GIT/epitrack-gui/data/cluster_ufgroup.csv")
ward_list <- DB_get_ward_with_nmet(con) #get ward list from epitrack database using a connection to database

#Convert into datatable
cluster_ufgroup <- data.table(cluster_ufgroup)
nb_isolate_per_iduf_group <- data.table(nb_isolate_per_iduf_group)
nb_nmet_per_iduf_group <- data.table(nb_nmet_per_iduf_group)

#Constant-----------------------------------------------------------------------

SEUIL <- 0.05 #p-value threshold
SAVED_DIRECTORY <- "./" #directory where tables and plots will be saved

#Start of test------------------------------------------------------------------

#summary table for all ward
ward_msd_result <- data.table(iduf_group = character(0), nb_sum_nmet = numeric(0), nb_date = numeric(0), nb_cluster = numeric(0))

for (ward in ward_list$iduf_group){
  print("-------")
  print(ward)
  
  #summary table for ward in process
  ward_msd_result_continue <- data.table(iduf_group = character(0), p.value = numeric(0), p.lo = numeric(0), p.center = numeric(0), p.hi = numeric(0), date_tau = as.Date(character(0)), tau = numeric(0), rate0 = numeric(0), rate1 = numeric(0))
  
  #Line to create time series
  #ward_period_list[[ward]] <- time_series_for_ward(ward = ward, date1 = '2023-01-01', date2 = '2023-10-22')
  #ward_period_list[[ward]] <- ward_period_list[[ward]][order(-date)]

  #Period list for ward in process
  current_data_table <- ward_period_list[[ward]]
  
  #Take the oldest date of data to init date_tau
  date_tau <- current_data_table[length(current_data_table$date),date]
  
  #Init alternative option for msd test
  alt <- "greater"
  
  #Revers the period list
  reversed_dates <- rev(as.character(ward_period_list[[ward]]$date))
  
  for (current_date in reversed_dates){
    
    current_data_table_cut <- current_data_table[date <= current_date & date >= date_tau]
    tmax <- nrow(current_data_table_cut)
    
    if (tmax != 0){
      
      msdresult <- msd(current_data_table_cut$sum_max_value, tmax, alternative = alt)

      if (msdresult["p.value"] < SEUIL){
        
        print(paste0(alt, " point detected"))
        
        #We change the alternative option
        if (alt == "greater"){
          alt <- "less"
        }
        else{
          alt <- "greater"
        }
        
        #We change the date_tau with the new date break point
        date_tau = as.Date(current_data_table_cut$date[1] - msdresult["tau"])
        
        #Complete summary table for ward in process
        new_row <- data.table(iduf_group = ward,
                              p.value = msdresult["p.value"],
                              p.lo = msdresult["p.lo"],
                              p.center = msdresult["p.center"],
                              p.hi = msdresult["p.hi"],
                              date_tau = date_tau,
                              tau = msdresult["tau"],
                              rate0 = msdresult["rate0"],
                              rate1 = msdresult["rate1"])
        ward_msd_result_continue <- rbind(ward_msd_result_continue, new_row)

      }
    }
  }
  
  #If break point find we make the plot
  if (length(ward_msd_result_continue$iduf_group) > 0){

    pdf(paste0(SAVED_DIRECTORY, ward, "_graphique.pdf"))

    p <- plot(ward_period_list[[ward]], xlim = rev(range(ward_period_list[[ward]]$date)),
         main = ward,
         xlab = "Date",
         ylab = "Somme des Nmet",
         pch = 1)
    for (date_tau_to_display in as.character(ward_msd_result_continue$date_tau)){
      abline(v = as.Date(date_tau_to_display), lwd = 2)

      cluster_list <- unique(cluster_ufgroup[iduf_group == ward]$idcluster)

      if(length(cluster_list > 0)){

        for (cluster in cluster_list){

          cluster_date <- setorder(cluster_ufgroup[idcluster == cluster], sample_date)
          abline(v = cluster_date[2, sample_date], col = "red", lty = 2)

        }
      }
    }

    dev.off()

    new_row <- data.table(iduf_group = ward,
                          nb_sum_nmet = sum(ward_period_list[[ward]]$sum_max_value > 0),
                          nb_date = length(ward_msd_result_continue$iduf_group),
                          nb_cluster = length(cluster_list)
                          )
    ward_msd_result <- rbind(ward_msd_result, new_row)

    write.csv(ward_msd_result_continue, file=paste0(SAVED_DIRECTORY, ward, "_ward_mst_result.csv"), quote = TRUE, na = "", row.names = FALSE, fileEncoding = "UTF-8")
  }
}

#We complete the summary table for all wards
ward_msd_result <- nb_isolate_per_iduf_group[ward_msd_result, on = "iduf_group"]
ward_msd_result <- nb_nmet_per_iduf_group[ward_msd_result, on = "iduf_group"]
ward_msd_result <- ward_msd_result[,c("iduf_group", "hospital_group", "label", "nb_isolate", "nb_nmet", "nb_sum_nmet", "nb_date", "nb_cluster")]
setnames(ward_msd_result, 'nb_sum_nmet', 'nb_sum_nmet>0')
write.csv(ward_msd_result, file=paste0(SAVED_DIRECTORY, "all_ward_mst_result.csv"), quote = TRUE, na = "", row.names = FALSE, fileEncoding = "UTF-8")
