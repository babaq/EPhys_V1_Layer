using NeuroAnalysis,FileIO,JLD2,Mmap,DataFrames,StatsBase,StatsPlots,ProgressMeter,Logging
import NaNMath
import Base: close

includet("Batch_SpikeGLX.jl")

function batchtests(tests::DataFrame,param::Dict{Any,Any}=Dict{Any,Any}();log::Dict{Any,AbstractLogger}=Dict{Any,AbstractLogger}(),plot::Bool=true)
    p = ProgressMeter.Progress(size(tests,1),desc="Batch ... ",start=-1)
    for t in eachrow(tests)
        next!(p,showvalues = [(:Test, t.files)])
        try
            if t.ID in ["Flash","Flash2Color"] && t.sourceformat=="SpikeGLX"
                process_flash_spikeglx(t.files,param;uuid=t.UUID,log,plot)
            else
                @warn "No corresonding process function, skip `$(t.ID)` from `$(t.sourceformat)` ..."
            end
        catch exc
            display("============================================")
            display("Error In Processing: $(t.files)")
            display("============================================")
            display.(stacktrace(catch_backtrace()))
        end
    end
    finish!(p)
    close(log)
end

close(log::Dict{Any,AbstractLogger}) = foreach(l->close(l.stream),values(log))

function (log::Dict{Any,AbstractLogger})(key,f,logfile)
    if !haskey(log,key)
        log[key] = SimpleLogger(open(logfile,"w"))
    end
    logger = log[key]
    if f isa String
        println(logger.stream,f)
    else
        with_logger(f,logger)
    end
    flush(logger.stream)
end
