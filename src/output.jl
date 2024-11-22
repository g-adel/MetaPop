using PrettyTables

function generatePrettyTable(data)
    # Extract metrics
    path_lengths = data["pathLengths"]
    spread_inf_days = data["spreadInfInd"]
    nPopulations = length(spread_inf_days)
    peak_inf_days = data["peakInfInd"]
    cumu_inf_mob =[0; [data["cumuInfMob"][i, i+1] for i in 1:size(data["cumuInfMob"], 1)-1] ]
    isolation_day = data["isolationInd"]
    second_day_infection = data["infectedHistory"][2,:]
        
    second_day_restrictions = data["second_day_restrictions"]
    # Create table data
    table_data = hcat(1:length(path_lengths), path_lengths, spread_inf_days, 
                    peak_inf_days,cumu_inf_mob,isolation_day, second_day_restrictions,second_day_infection)
    header = ["Population", "Path Length", "Spread Day","Peak Inf. Day", 
    "Cumu Incoming Infected","isolation_day", "second day restriction",
      "second day infection"
     ]

    # Print the table
    println()
    pretty_table(table_data, header=header, crop = :none)
    # println( [data["metaHist"][t].populations[2].ρs[1] for t in 1:length(data["metaHist"])]')
end
