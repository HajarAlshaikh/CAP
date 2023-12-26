import Plots as plt
using PGFPlots
using XLSX
using DataFrames

input_file = "Correlation  approximation data result.xlsx"

data1 = XLSX.readtable(input_file, "PT_COR") |> DataFrame
describe(data1)
methods = repeat(["Pro", "SQV", "SQC"], inner=10)
sizes = [repeat([5, 9, 15, 25, 50, 100, 200, 500, 1000, 2000], outer=3)]
sizes = vcat(sizes...)
values = data1[:, 3]
df = DataFrame(Method=methods, Size=sizes, Value=values)

# Create label array
labels = permutedims(["Pro", "SQV", "SQC"])

# Get unique sizes
unique_sizes = unique(df.Size)
unique_Method = unique(df.Method)

function getDataByMethod(df, Method)
    filter(r -> r[:Method] == Method, df)
end

function getPlotingData(df)
    select(df, :Value).Value
end

function getPlot(df)
    data2 = getDataByMethod(df, "Pro") |> getPlotingData
    data3 = getDataByMethod(df, "SQV") |> getPlotingData
    data4 = getDataByMethod(df, "SQC") |> getPlotingData

    plt.plot(data2, marker=:circle, label="Pro")
    plt.plot!(data3, marker=:circle, label="SQV")
    plt.plot!(data4, xticks=(1:length(unique_sizes), unique_sizes), title="TIME", legend=:top, legend_column=:3,
        xtickfont=(10), titlefont=("Arial", 12, :bold), xrotation=45, marker=:circle, label="SQC")
    plt.xlabel!("Matrix Size")

end

plt.pgfplotsx()

pt = getPlot(df)

plt.savefig(pt, "PTCor.tex")

data2 = XLSX.readtable(input_file, "PI_COR") |> DataFrame
describe(data2)
methods = repeat(["Pro", "SQV", "SQC"], inner=10)
sizes = [repeat([5, 9, 15, 25, 50, 100, 200, 500, 1000, 2000], outer=3)]
sizes = vcat(sizes...)
values = data2[:, 3]
df = DataFrame(Method=methods, Size=sizes, Value=values)

# Create label array
labels = permutedims(["Pro", "SQV", "SQC"])

# Get unique sizes
unique_sizes = unique(df.Size)
unique_Method = unique(df.Method)

function getDataByMethod(df, Method)
    filter(r -> r[:Method] == Method, df)
end

function getPlotingData(df)
    select(df, :Value).Value
end

function getPlot(df)
    data5 = getDataByMethod(df, "Pro") |> getPlotingData
    data6 = getDataByMethod(df, "SQV") |> getPlotingData
    data7 = getDataByMethod(df, "SQC") |> getPlotingData

    plt.plot(data5, marker=:circle, label="Pro")
    plt.plot!(data6, marker=:circle, label="SQV")
    plt.plot!(data7, xticks=(1:length(unique_sizes), unique_sizes), title="ITER", legend=:top, legend_column=:3,
        xtickfont=(10), titlefont=("Arial", 12, :bold), xrotation=45, marker=:circle, label="SQC")
    plt.xlabel!("Matrix Size")

end

plt.pgfplotsx()

pi = getPlot(df)

plt.savefig(pi, "PICor.tex")

data3 = XLSX.readtable(input_file, "PE_COR") |> DataFrame
describe(data3)
methods = repeat(["Pro", "SQV", "SQC"], inner=10)
sizes = [repeat([5, 9, 15, 25, 50, 100, 200, 500, 1000, 2000], outer=3)]
sizes = vcat(sizes...)
values = data3[:, 3]
df = DataFrame(Method=methods, Size=sizes, Value=values)

# Create label array
labels = permutedims(["Pro", "SQV", "SQC"])

# Get unique sizes

unique_sizes = unique(df.Size)
unique_Method = unique(df.Method)

function getDataByMethod(df, Method)
    filter(r -> r[:Method] == Method, df)
end

function getPlotingData(df)
    select(df, :Value).Value
end

function getPlot(df)
    data5 = getDataByMethod(df, "Pro") |> getPlotingData
    data6 = getDataByMethod(df, "SQV") |> getPlotingData
    data7 = getDataByMethod(df, "SQC") |> getPlotingData

    plt.plot(data5, marker=:circle, label="Pro")
    plt.plot!(data6, marker=:circle, label="SQV")
    plt.plot!(data7, xticks=(1:length(unique_sizes), unique_sizes), title="Error", legend=:top, legend_column=:3,
        xtickfont=(10), titlefont=("Arial", 12, :bold), xrotation=45, marker=:circle, label="SQC")
    plt.xlabel!("Matrix Size")

end
plt.pgfplotsx()

pe = getPlot(df)

plt.savefig(pe, "PECor.tex")