using NeuroAnalysis,FileIO,JLD2,Statistics,StatsPlots,StatsBase,Images,ProgressMeter,DataFrames,XLSX,Dierckx

function mergespikeinfo!(indir;unit=Dict(),datafile="spike.jld2")
    for (root,dirs,files) in walkdir(indir)
        if datafile in files
            spike,siteid = load(joinpath(root,datafile),"spike","siteid")
            foreach(k->delete!(spike,k),["unitspike","isspikesorted","unitsync","t0"])
            if haskey(unit,"siteid")
                if siteid == unit["siteid"]
                    newids = setdiff(spike["unitid"],unit["unitid"])
                    if !isempty(newids)
                        newididx = indexin(newids,spike["unitid"])
                        append!(unit["unitid"],newids)
                        append!(unit["unitgood"],spike["unitgood"][newididx])
                        unit["unitposition"] = [unit["unitposition"];spike["unitposition"][newididx,:]]
                        unit["unitwaveforms"] = [unit["unitwaveforms"];spike["unitwaveforms"][newididx,:,:]]
                        unit["unitwaveform"] = [unit["unitwaveform"];spike["unitwaveform"][newididx,:]]
                        foreach(k->append!(unit["unitfeature"][k],spike["unitfeature"][k][newididx]),keys(spike["unitfeature"]))
                        append!(unit["unittemplatemeanamplitude"],spike["unittemplatemeanamplitude"][newididx])
                        unit["unittemplateposition"] = [unit["unittemplateposition"];spike["unittemplateposition"][newididx,:]]
                        unit["unittemplatewaveform"] = [unit["unittemplatewaveform"];spike["unittemplatewaveform"][newididx,:]]
                        foreach(k->append!(unit["unittemplatefeature"][k],spike["unittemplatefeature"][k][newididx]),keys(spike["unittemplatefeature"]))
                        foreach(k->append!(unit["qm"][k],spike["qm"][k][newididx]),keys(spike["qm"]))
                    end
                end
            else
                unit = spike
                unit["siteid"] = siteid
            end
        end
    end
    return unit
end

resultroot = "Z:/"


## Merge ALL Spike Info of a RecordSite
penetration = DataFrame(XLSX.readtable(joinpath(resultroot,"penetration.xlsx"),"Sheet1"))
@showprogress "Batch Merging Spike Info ... " for r in eachrow(penetration)
    indir = joinpath(resultroot,r.Subject_ID,r.siteid)
    unit = mergespikeinfo!(indir)
    save(joinpath(indir,"unit.jld2"),"unit",unit)
end


## Spike Info of a RecordSite
figfmt = [".svg"]
indir = joinpath(resultroot,"AG2","AG2_V1_ODL3")
layer = load(joinpath(indir,"layer.jld2"),"layer")

plotdepthfeature = (feature,position,ks;good=trues(size(position,1)),kw=ones(1,length(ks)),size=(350,700),xlabel="",df=nothing,layer=nothing,grid=true,xticks=:auto,vl=[]) -> begin
    kn = length(ks)
    vs = hcat(map(k->feature[k][good],ks)...).*kw
    cs = permutedims(palette(:default).colors.colors[1:kn])
    p=plot(;size,grid,tickdir=:out)
    isempty(vl) || vline!(p,vl;linecolor=:gray10,legend=false,lw=1)
    scatter!(p,vs,position[good,2];label=ks,markerstrokewidth=0,alpha=0.6,markersize=3,xticks,
        xlabel,ylabel="Depth (μm)",color=cs,left_margin=4Plots.mm)

    yms = [unitdensity(position[good,2],w=vs[:,i],wfun=mean,bw=40,step=20,s=1.4) for i in 1:kn]
    y = yms[1].y
    yms = hcat(map(i->i.n,yms)...)
    if df isa Dict
        if haskey(df,"feature")
            df["feature"] = [df["feature"] ks]
            df["depth"] = y
            df["depthfeature"] = [df["depthfeature"] yms]
        else
            df["feature"] = ks
            df["depth"] = y
            df["depthfeature"] = yms
        end
    end

    plot!(p,[zeros(1,kn);yms;zeros(1,kn)],[minimum(y);y;maximum(y)], st=:shape,lw=0,label=false,alpha=0.2,color=cs)
    plot!(p,yms,y;label=false,lw=2,color=cs)
    if !isnothing(layer)
        xmin,xmax=extrema(vs)
        xm = 0.02(xmax-xmin)
        xmin-=7xm; xmax+=xm
        ann = [(xmin+0.02(xmax-xmin),mean(layer[k]),text(k,6,:gray10,:left,:vcenter)) for k in keys(layer)]
        hline!(p,[l[2] for l in values(layer)];linecolor=:gray25,leg=false,lw=0.5,ann,xlims=(xmin,xmax),
        yticks=0:200:3820)
    end
    p
end

spikeinfo = (indir;layer=nothing,figfmt=[".png"]) -> begin

    unit = load(joinpath(indir,"unit.jld2"),"unit")

    unitid = unit["unitid"];unitgood=unit["unitgood"]
    unittempposition=unit["unittemplateposition"];unittempamp = unit["unittemplatemeanamplitude"]
    unittempwave=unit["unittemplatewaveform"];unittempfeature=unit["unittemplatefeature"]
    unitposition=unit["unitposition"];unitwave = unit["unitwaveform"];unitfeature = unit["unitfeature"]
    unitqm=unit["qm"];siteid = unit["siteid"]

    if layer == :batch
        layerpath = joinpath(indir, "layer.jld2")
        layer = ispath(layerpath) ? load(layerpath, "layer") : nothing
    end

    # Unit Position
    plotunitposition(unitposition;unitgood,chposition=unit["chposition"],size=(400,700),layer)
    foreach(ext -> savefig(joinpath(indir, "$(isnothing(layer) ? "" : "Layer_")UnitPosition$ext")), figfmt)

    # Unit Density
    muc=1.2
    r=100
    w = replace(unitgood,0=>muc)

    # n,y = unitdensity(unitposition[unitgood,2];bw=40,step=20,r,s=1.4)
    n,y = unitdensity(unitposition[:,2];w,bw=40,step=20,r,s=1.4)
    plot(1e9n,y;xlabel="Density (unit/mm³)",ylabel="Depth (μm)",leg=false,size=(350,700),grid=false,tickdir=:out,
        lw=2,left_margin=4Plots.mm,title="Count(MU)=$muc, r=$(r)μm")
    if !isnothing(layer)
        xmin,xmax = extrema(1e9n)
        xm = 0.02(xmax-xmin)
        xmin-=xm;xmax+=xm
        ann = [(xmin+0.02(xmax-xmin),mean(layer[k]),text(k,6,:gray10,:left,:vcenter)) for k in keys(layer)]
        hline!([l[2] for l in values(layer)];linecolor=:gray25,leg=false,lw=0.5,ann,xlims=(xmin,xmax),
        yticks=0:200:3820)
    end
    foreach(ext -> savefig(joinpath(indir, "$(isnothing(layer) ? "" : "Layer_")UnitDensity$ext")), figfmt)

    df=Dict("feature"=>["Density";;],"depth"=>y,"depthfeature"=>[1e9n;;],"siteid"=>siteid)

    # Unit Feature
    un,wn = size(unitwave)
    wx = range(-wn/2,step=1,length=wn)

    uwy = [1e-2*unitwave[i,j]+unitposition[i,2] for i in 1:un,j in 1:wn]
    uwx = [0.05*wx[j]+unitposition[i,1] for i in 1:un,j in 1:wn]

    plot(uwx',uwy',leg=false,size=(350,700),xlabel="X (μm)",ylabel="Y (μm)",grid=false,lw=1,left_margin=4Plots.mm,alpha=0.6,tickdir=:out,
        color=permutedims(map(i->i ? RGB(0,0.55,0) : RGB(0),unitgood)),title="Unit Waveform")
    if !isnothing(layer)
        xmin=1.5;xmax=43
        ann = [(xmin+0.02(xmax-xmin),mean(layer[k]),text(k,6,:gray10,:left,:vcenter)) for k in keys(layer)]
        hline!([l[2] for l in values(layer)];linecolor=:gray25,leg=false,lw=0.5,ann,xlims=(xmin,xmax),
        yticks=0:200:3820)
    end
    foreach(ext -> savefig(joinpath(indir, "$(isnothing(layer) ? "" : "Layer_")UnitWaveform$ext")), figfmt)


    # n=findfirst(unitid.==432)
    # plot!(uwx[n,:],uwy[n,:],leg=false,lw=1,color=:red)
    # plot(0.1*(-40:40) .+ unit["chposition"][:,1]',  0.03*unit["unitwaveforms"][n,:,:] .+ unit["chposition"][:,2]';
    #     leg=false,ratio=0.1,size=(600,3840),color=:dodgerblue,lw=2,frame=:none)

    plotdepthfeature(unittempfeature,unitposition,["upspread" "downspread"];xlabel="Spread (μm)",layer,grid=false,vl=[0])
    foreach(ext -> savefig(joinpath(indir, "$(isnothing(layer) ? "" : "Layer_")Spread$ext")), figfmt)
    plotdepthfeature(unitfeature,unitposition,["duration"];kw=1000,xlabel="Time (ms)",layer,grid=false,vl=[0])
    foreach(ext -> savefig(joinpath(indir, "$(isnothing(layer) ? "" : "Layer_")Duration$ext")), figfmt)
    plotdepthfeature(unitfeature,unitposition,["peaktroughratio"];xlabel="Peak/Trough",layer,grid=false,xticks=-1:-2:-15,vl=[-1])
    foreach(ext -> savefig(joinpath(indir, "$(isnothing(layer) ? "" : "Layer_")PTRatio$ext")), figfmt)

    # plotdepthfeature(unittempfeature,unitposition,["upspread" "downspread" "leftspread" "rightspread"];xlabel="Spread (μm)",df,layer,good=unitgood)
    plotdepthfeature(unittempfeature,unitposition,["upspread" "downspread" "leftspread" "rightspread"];xlabel="Spread (μm)",df,layer)
    foreach(ext -> savefig(joinpath(indir, "$(isnothing(layer) ? "" : "Layer_")UnitFeature_spread$ext")), figfmt)

    # plotdepthfeature(unittempfeature,unitposition,["uppvinv" "downpvinv"],kw=1000;xlabel="1/Propagation Speed (ms/μm)",df,layer,good=unitgood)
    plotdepthfeature(unittempfeature,unitposition,["uppvinv" "downpvinv"],kw=1000;xlabel="1/Propagation Speed (ms/μm)",df,layer)
    foreach(ext -> savefig(joinpath(indir, "$(isnothing(layer) ? "" : "Layer_")UnitFeature_invspeed$ext")), figfmt)


    # plotdepthfeature(unitfeature,unitposition,["amplitude" "peaktroughratio"];kw=[1e-3 1],xlabel="A.U.",df,layer,good=unitgood)
    plotdepthfeature(unitfeature,unitposition,["amplitude" "peaktroughratio"];kw=[1e-3 1],xlabel="A.U.",df,layer)
    foreach(ext -> savefig(joinpath(indir, "$(isnothing(layer) ? "" : "Layer_")UnitFeature_amp$ext")), figfmt)

    # plotdepthfeature(unitfeature,unitposition,["halftroughwidth" "halfpeakwidth" "duration"];kw=1000,xlabel="Time (ms)",df,layer,good=unitgood)
    plotdepthfeature(unitfeature,unitposition,["halftroughwidth" "halfpeakwidth" "duration"];kw=1000,xlabel="Time (ms)",df,layer)
    foreach(ext -> savefig(joinpath(indir, "$(isnothing(layer) ? "" : "Layer_")UnitFeature_width$ext")), figfmt)

    # plotdepthfeature(unitfeature,unitposition,["repolarrate" "recoverrate"];kw=1e-3,xlabel="A.U.",df,layer,good=unitgood)
    plotdepthfeature(unitfeature,unitposition,["repolarrate" "recoverrate"];kw=1e-3,xlabel="A.U.",df,layer)
    foreach(ext -> savefig(joinpath(indir, "$(isnothing(layer) ? "" : "Layer_")UnitFeature_rate$ext")), figfmt)

    # plotdepthfeature(unitqm,unitposition,["fr" "fp" "pisi"];xlabel="A.U.",df,layer,good=unitgood)
    plotdepthfeature(unitqm,unitposition,["fr" "fp" "pisi"];xlabel="A.U.",df,layer)
    foreach(ext -> savefig(joinpath(indir, "$(isnothing(layer) ? "" : "Layer_")UnitFeature_qm$ext")), figfmt)

    save(joinpath(indir,"unitdepthfeature.jld2"),"df",df)
end

## Batch Spike Info
@showprogress "Batch Spike Info ... " for r in eachrow(penetration)
    # spikeinfo(joinpath(resultroot,r.Subject_ID,r.siteid))
    spikeinfo(joinpath(resultroot,r.Subject_ID,r.siteid),layer=:batch)
end
