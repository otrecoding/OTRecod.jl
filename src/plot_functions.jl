export plot_benchmark, plot_scenario

"""
    plot_benchmark(datapath,outputpath,ymin=0.0, ymax=0.8)

Boxplot of each simulation with all methods
"""
function plot_benchmark(datapath, outputpath, ymin=0.0, ymax=0.8)

    datadirlist = readdir(datapath);
    outfilelist = readdir(outputpath);
    for dataset in datadirlist
        datasetpath = string(datapath,"/",dataset);

        if !isdir(datasetpath) continue end
        if (dataset == "Sn") continue end
        if (dataset == "Spi-4") continue end
        println(datasetpath);

        nbmethods = 0;
        dataplot = [];
        method = [];
        for outfile in outfilelist
            if occursin(string(dataset,"-"), outfile) & occursin(".out", outfile)
                nbmethods += 1;
                outfilepath = string(outputpath,"/",outfile);
                data = readdlm(outfilepath, ',');

                colindex = findall([strip(data[1,j]) for j in 1:size(data)[2]] .== "eDistravg");
                colmethod = findall([strip(data[1,j]) for j in 1:size(data)[2]] .== "method");
                colrelax = findall([strip(data[1,j]) for j in 1:size(data)[2]] .== "relax");
                if (isempty(colindex) | isempty(colmethod) | isempty(colrelax))
                    println("A required field is missing from ", outfile);
                    continue
                end

                if isempty(dataplot)
                    dataplot = Array{Float64}(data[2:end,colindex[1]]);
                    method = Array{String,1}([string(strip(data[2,colmethod[1]]), data[2,colrelax[1]])]);
                else
                    dataplot = [dataplot Array{Float64}(data[2:end,colindex[1]])];
                    method = [method  Array{String,1}([string(strip(data[2,colmethod[1]]), data[2,colrelax[1]])])];
                end
            end
        end
        method_tex = method;
        for i in 1:size(method,2)
            if method[1,i] == "group-0.0" method_tex[1,i] = "\\modgroup"
            elseif method[1,i] == "group-0.4" method_tex[1,i] = "\\modgrouprich{}"
            elseif method[1,i] == "joint-0.0" method_tex[1,i] ="\\modjoint{}"
            elseif method[1,i] == "joint-0.4" method_tex[1,i] = "\\modjointrich{}"
            end
        end
        indlist = sortperm(method[1,:])
        method[1,:] = method[1,indlist];
        method_tex[1,:] = method_tex[1,indlist];
        dataplot = dataplot[:,indlist];
        data_avg = mean(dataplot,dims=1);
        println(method_tex)
        println(data_avg)
        titleplot = string("Compare every method on scenario ", dataset);
        plt=StatsPlots.boxplot(method_tex,dataplot,title=titleplot,notch=false,xlab="Method", ylab="Error in the distribution",legend=false,color_palette=:grays, ylims=(ymin,ymax))
        Plots.savefig(plt, string(outputpath,"/",dataset,".tex"));
        # Plots.savefig(plt, string(outputpath,"/",dataset,".pdf"));
    end
end

function plot_scenario(outputpath, scenario, method, ymin=0.0, ymax=0.4)
    outfilelist = readdir(outputpath);
    dataplots = [];
    datasets = [];
    datanum = [];
    for outfile in outfilelist
        if occursin(scenario, outfile) & occursin(".out", outfile)
            if occursin(method, outfile)
                outfilepath = string(outputpath,"/",outfile);
                data = readdlm(outfilepath, ',');
                colindex = findall([strip(data[1,j]) for j in 1:size(data)[2]] .== "eDistravg");
                colmethod = findall([strip(data[1,j]) for j in 1:size(data)[2]] .== "method");
                colrelax = findall([strip(data[1,j]) for j in 1:size(data)[2]] .== "relax");
                if (isempty(colindex) | isempty(colmethod) | isempty(colrelax))
                    println("A required field is missing from ", outfile);
                    continue
                end
                if isempty(dataplots)
                    dataplots = Array{Float64}(data[2:101,colindex[1]]);
                    datasets = Array{String,1}([replace(outfile, '-' * method * ".out" => "")]);
                    splitstr = split(outfile,"-");
                    datanum = Array{Float64}([parse(Float64,splitstr[2])]);
                else
                    dataplots = hcat(dataplots, Array{Float64}(data[2:101,colindex[1]]));
                    datasets = hcat(datasets, [replace(outfile, '-' * method * ".out" => "")]);
                    splitstr = split(outfile,"-");
                    datanum = hcat(datanum, parse(Float64,splitstr[2]));
                end
            end
        end
    end
    indlist = sortperm(datanum[1,:])
    datasets[1,:] = datasets[1,indlist];
    dataplots = dataplots[:,indlist];

    for i in 1:size(datasets,2)
        if (datasets[1,i] == "SNL-0") datasets[1,i]= "\\Sref{}"
        elseif (datasets[1,i] == "SNL-3") datasets[1,i]= "\\text{SHG}"
        elseif (datasets[1,i] == "SR-0.5") datasets[1,i]= "\\Sref{}"
        elseif (datasets[1,i] == "Sa-0") datasets[1,i]= "\\Sref{}"
        elseif (datasets[1,i] == "SX-2.5") datasets[1,i]= "\\Sref{}"
        elseif (datasets[1,i] == "Sn-1000") datasets[1,i]= "\\Sref{}"
        else datasets[1,i] = "\\text{" * datasets[1,i] * "}"
        end
    end

    println(datasets)
    method_tex = method;
    if method == "group-basic" method_tex="\\modgroup"
    elseif method == "group-0.4-0.0" method_tex = "\\modgrouprich{}"
    elseif method == "joint-basic" method_tex="\\modjoint{}"
    elseif method == "joint-0.4-0.1" method_tex = "\\modjointrich{}"
    else  end

    titleplot = string("Results of ", method_tex);
    plt=StatsPlots.boxplot(datasets,dataplots,size=(400,500),whisker_width=0.25,title=titleplot,notch=true,legend=false,color_palette=:grays, ylims=(ymin,ymax), tickfontsize=14, titlefontsize=18)
    if scenario == "Sn"
        plt=StatsPlots.boxplot(datasets,dataplots,whisker_width=0.25,title=titleplot,notch=true,legend=false,color_palette=:grays, ylims=(ymin,ymax)) #xlab="Scenario", ylab="Error in the distribution",
    end
    Plots.savefig(plt, outputpath * "/" * scenario * '-' * method * ".tex");
    # Plots.savefig(plt, outputpath * "/" * scenario * '-' * method * ".pdf");

end

function plot_all_scenarios()

    mkpath("OutfilesJO")
    plot_scenario("OutfilesJO", "SX", "group-basic", 0.0, 0.7)
    plot_scenario("OutfilesJO", "SX", "group-0.4-0.0", 0.0, 0.7)
    plot_scenario("OutfilesJO", "SX", "joint-basic", 0.0, 0.7)
    plot_scenario("OutfilesJO", "SX", "joint-0.4-0.1", 0.0, 0.7)

    plot_scenario("OutfilesJO", "SR", "group-basic", 0.0, 1.0)
    plot_scenario("OutfilesJO", "SR", "group-0.4-0.0", 0.0, 1.0)
    plot_scenario("OutfilesJO", "SR", "joint-basic", 0.0, 1.0)
    plot_scenario("OutfilesJO", "SR", "joint-0.4-0.1", 0.0, 1.0)

    plot_scenario("OutfilesJO", "Sa", "group-basic", 0.0, 0.6)
    plot_scenario("OutfilesJO", "Sa", "group-0.4-0.0", 0.0, 0.6)
    plot_scenario("OutfilesJO", "Sa", "joint-basic", 0.0, 0.6)
    plot_scenario("OutfilesJO", "Sa", "joint-0.4-0.1", 0.0, 0.6)

    plot_scenario("OutfilesJO", "SNL", "group-basic", 0.0, 0.7)
    plot_scenario("OutfilesJO", "SNL", "group-0.4-0.0", 0.0, 0.7)
    plot_scenario("OutfilesJO", "SNL", "joint-basic", 0.0, 0.7)
    plot_scenario("OutfilesJO", "SNL", "joint-0.4-0.1", 0.0, 0.7)

    plot_scenario("OutfilesJO/Sn", "Sn", "group-basic", 0.0, 0.8)
    plot_scenario("OutfilesJO/Sn", "Sn", "group-0.4-0.0", 0.0, 0.8)
    plot_scenario("OutfilesJO/Sn", "Sn", "joint-basic", 0.0, 0.8)
    plot_scenario("OutfilesJO/Sn", "Sn", "joint-0.4-0.1", 0.0, 0.8)

end
