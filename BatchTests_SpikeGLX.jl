includet("Batch/Batch.jl")


## Param and Metadata
param = Dict{Any,Any}(
    :dataroot => "X:/",
    :dataexportroot => "Y:/",
    :resultroot => "Z:/",
    :stimuliroot => "S:/",
    :spikesorter => "kilosort3")
meta = readmeta(joinpath(param[:dataexportroot],"metadata.mat"))


## Query Tests
tests = select!(filter(meta) do r
                    r.Subject_ID in ["AG5"] &&
                    # r.RecordSession == "V1" &&
                    # r.RecordSite == "ODR1" &&
                    r.ID in ["Flash2Color"] &&
                    r.sourceformat == "SpikeGLX"
                end,
                [:files,:ID,:UUID,:sourceformat])

## Batch Tests
batchtests(tests,param,plot=true)
