using NeuroAnalysis,Statistics,StatsBase,FileIO,JLD2,StatsPlots,Images,DSP,ProgressMeter,DataFrames,XLSX,Dierckx,HypothesisTests

dataroot = "X:/"
dataexportroot = "Y:/"
resultroot = "Z:/"

cohcm = cgrad([HSL(170,1,0.94),HSL(190,1,0.56),HSL(233,1,0.38)])
powcm = cgrad([HSL(60,1,0.92),HSL(40,1,0.62),HSL(0,1,0.32)])
csdcm = cgrad([HSL(358,0.7,0.3),HSL(358,0.8,0.5),HSL(58,0.9,0.6),HSL(125,0.4,0.7),HSL(192,0.9,0.6),HSL(238,0.8,0.5),HSL(238,0.7,0.3)],
              [0,0.12,0.37,0.5,0.63,0.88,1])
huecm = Dict("DKL_HueL0"=>ColorMaps["lidkl_mcchue_l0"],"HSL_HueYm"=>ColorMaps["hsl_mshue_l0.4"])
absmax(x) = mapreduce(i->maximum(abs.(i)),max,x)
absmin(x) = mapreduce(i->minimum(abs.(i)),min,x)
abspct(x,p=99.9) = percentile(vec(abs.(x)),p)
abspct2(x,p=99.9) = percentile.((vec(abs.(x)),),(100-p,p))

plotlayerunitfeature=(udf,depths;w=140,h=600,color=:tab10,layer=nothing)->begin
    n = 4;cs = palette(color).colors.colors[1:n]
    p=plot(layout=(1,n),link=:y,legend=false,grid=false,size=(n*w,h))
    for i in 1:n
        yticks = i==1 ? (0:200:depths[end]) : false
        leftmargin = i==1 ? 4Plots.mm : -4Plots.mm
        bottommargin = 3Plots.mm

        if i == 1
            x = udf.Density
            xmin,xmax = extrema(x)
            plot!(p[i],x,depths;xlabel="Density",ylabel="Depth (μm)",leg=false,color=cs[i],
            lw=1.5,yticks,xticks=[ceil(xmin,sigdigits=1),floor(xmax,sigdigits=1)],
            leftmargin,bottommargin,tickor=:out)
        elseif i == 2
            x = [udf.upspread;reverse(udf.downspread)]
            y = [depths;reverse(depths)]
            plot!(p[i],x,y;st=:shape,xlabel="Spread",lw=0,yticks,alpha=0.8,xticks=[0],color=cs[i],
            leftmargin,bottommargin,tickor=:out)
        elseif i == 3
            x = [0;udf.peaktroughratio;0]
            y = [depths[begin];depths;depths[end]]
            plot!(p[i],x,y;st=:shape,xlabel="PTRatio",lw=0,yticks,alpha=0.8,xticks=[0],color=cs[i],
            leftmargin,bottommargin,tickor=:out)
        elseif i == 4
            x = [0;udf.duration;0]
            y = [depths[begin];depths;depths[end]]
            plot!(p[i],x,y;st=:shape,xlabel="Duration",lw=0,yticks,alpha=0.8,xticks=[0],color=cs[i],
            leftmargin,bottommargin,tickor=:out)
        end

        xmin,xmax = extrema(x)
        if !isnothing(layer)
            ann = [(xmin+0.02(xmax-xmin),mean(layer[k]),text(k,6,:gray10,:left,:vcenter)) for k in keys(layer)]
            hline!(p[i],[l[end] for l in values(layer)];linecolor=:gray25,legend=false,lw=0.5,ann)
        end
    end
    p
end

plotlayercondresponse=(resps,times,depths;colors=[],minmaxcolors=[],colormaps=[],xl="Time (ms)",xt=[1000,2000],rl="",rlims=(0,1),rscale=1,
                       rcolor=HSL(120,0.5,0.25),titles="",titlefontsize=7,w=140,h=600,color=:coolwarm,layer=nothing,layerext=nothing)->begin
    n = length(resps)
    p=plot(layout=(1,n),link=:y,legend=false,grid=false,size=(n*w,h))
    for i in 1:n
        yticks = i==1 ? (0:200:depths[i][end]) : false
        xticks = xl == "Time (ms)" ? (0:50:times[i][end]) : xt
        leftmargin = i==1 ? 6Plots.mm : -4Plots.mm
        xlabel = i==1 ? xl : ""
        rlabel = i==1 ? rl : ""
        ylabel = i==1 ? "Depth (μm)" : ""
        title = isempty(titles) ? titles : titles[i]

        if isnothing(layerext)
            dr = depths[i]
            rp = resps[i]
        else
            di = (layer["WM"][end]-layerext).<=depths[i].<=(layer["1"][end]+layerext)
            dr = depths[i][di]
            @views rp = resps[i][di,:]
        end
        if xl == "Time (ms)"
            @views lim = isnothing(layer) ? absmax(resps) : abspct(resps[i][layer["WM"][end].<=depths[i].<=layer["1"][end],:])
            clims = (-lim,lim)
        else
            @views clims = isnothing(layer) ? (absmin(resps),absmax(resps)) : abspct2(resps[i][layer["WM"][end].<=depths[i].<=layer["1"][end],:])
        end
        xmin,xmax = extrema(times[i])
        ymin,ymax = extrema(dr)
        annx = xmin+0.02(xmax-xmin)
        anny = ymin+0.01(ymax-ymin)
        adx = 30
        ann = isempty(colors) ? [] : [(annx,anny,Plots.text("■",15,colors[i],:left,:bottom))]
        ann = isempty(minmaxcolors) ? ann : [(annx,anny,Plots.text("▮",15,minmaxcolors[i][2],:left,:bottom)),(annx+adx,anny,Plots.text("▮",15,minmaxcolors[i][1],:left,:bottom))]
        ann = isempty(colormaps) ? ann : [(annx+adx*(j-1),anny,Plots.text("▮",8,colormaps[i][j],:left,:bottom)) for j in eachindex(colormaps[i])]

        heatmap!(p[i],times[i],dr,rp;color,clims,ylims=extrema(dr),
        title,titlefontsize,yticks,xticks,tickor=:out,xlabel,ylabel,
        ann,leftmargin,bottommargin=3Plots.mm,topmargin=12Plots.mm)

        if xl != "Time (ms)"
            rm = dropdims(mean(rp*rscale,dims=2),dims=2)
            rmse = dropdims(std(rp*rscale,dims=2),dims=2)/sqrt(size(rp,2))
            rticks=[rlims[1],round(0.5(rlims[2]-rlims[1]),sigdigits=1)]
            pp=mytwiny(p[i])
            plot!(pp,[rm.-rmse;reverse(rm.+rmse)],[dr;reverse(dr)]; st=:shape,lw=0,alpha=0.25,color=rcolor,leftmargin,xlabel=rlabel)
            plot!(pp,rm,dr;color=rcolor,leg=false,tickor=:out,xlims=rlims,xticks=rticks,ylims=extrema(dr),lw=1,leftmargin)
        end

        if !isnothing(layer)
            ann = [(annx,mean(layer[k]),text(k,6,:gray10,:left,:vcenter)) for k in keys(layer)]
            hline!(p[i],[l[end] for l in values(layer)];linecolor=:gray25,legend=false,lw=0.5,ann)
        end
    end
    p
end

plotlayertrialresponse=(resps,trials,depths;colors=[],minmaxcolors=[],colormaps=[],xl="Trial",rl="",rlims=(0,1),rscale=1,
                        rcolor=HSL(120,0.5,0.25),titles="",titlefontsize=7,w=140,h=600,color=:heat,layer=nothing,layerext=nothing)->begin
    n = length(resps)
    clims = (absmin(resps),absmax(resps))
    p=plot(layout=(1,n),link=:y,legend=false,grid=false,size=(n*w,h))
    for i in 1:n
        @views isnothing(layer) || (clims=abspct2(resps[i][layer["WM"][end].<=depths[i].<=layer["1"][end],:]))
        yticks = i==1 ? (0:200:depths[i][end]) : false
        xticks = range(start=1,stop=round(Int,0.75*trials[i][end]),length=2)
        leftmargin = i==1 ? 6Plots.mm : -4Plots.mm
        xlabel = i==1 ? xl : ""
        rlabel = i==1 ? rl : ""
        ylabel = i==1 ? "Depth (μm)" : ""
        title = isempty(titles) ? titles : titles[i]

        if isnothing(layerext)
            dr = depths[i]
            rp = resps[i]
        else
            di = (layer["WM"][end]-layerext).<=depths[i].<=(layer["1"][end]+layerext)
            dr = depths[i][di]
            @views rp = resps[i][di,:]
        end
        xmin,xmax = extrema(trials[i])
        ymin,ymax = extrema(dr)
        annx = xmin+0.02(xmax-xmin)
        anny = ymin+0.01(ymax-ymin)
        adx = 30
        ann = isempty(colors) ? [] : [(annx,anny,Plots.text("▮",15,colors[2+2(i-1)],:left,:bottom)),(annx+adx,anny,Plots.text("▮",15,colors[1+2(i-1)],:left,:bottom))]
        ann = isempty(minmaxcolors) ? ann : [(annx,anny,Plots.text("▮",15,minmaxcolors[i][2],:left,:bottom)),(annx+adx,anny,Plots.text("▮",15,minmaxcolors[i][1],:left,:bottom))]
        ann = isempty(colormaps) ? ann : [(annx+adx*(j-1),anny,Plots.text("▮",8,colormaps[i][j],:left,:bottom)) for j in eachindex(colormaps[i])]

        heatmap!(p[i],trials[i],dr,rp;color,clims,ylims=extrema(dr),
        title,titlefontsize,yticks,xticks,tickor=:out,xlabel,ylabel,
        ann,leftmargin,bottommargin=3Plots.mm)

        rm = dropdims(mean(rp*rscale,dims=2),dims=2)
        rmse = dropdims(std(rp*rscale,dims=2),dims=2)/sqrt(size(rp,2))
        rticks=[rlims[1],round(0.5(rlims[2]-rlims[1]),sigdigits=1)]
        pp=mytwiny(p[i])
        plot!(pp,[rm.-rmse;reverse(rm.+rmse)],[dr;reverse(dr)]; st=:shape,lw=0,alpha=0.25,color=rcolor,leftmargin,xlabel=rlabel)
        plot!(pp,rm,dr;color=rcolor,leg=false,tickor=:out,xlims=rlims,xticks=rticks,ylims=extrema(dr),lw=1,leftmargin)

        if !isnothing(layer)
            ann = [(annx,mean(layer[k]),text(k,6,:gray10,:left,:vcenter)) for k in keys(layer)]
            hline!(p[i],[l[end] for l in values(layer)];linecolor=:gray25,legend=false,lw=0.5,ann)
        end
    end
    p
end

function mytwiny(sp, letter=:y)
    plt = sp.plt
    # orig_sp = first(plt.subplots)
    orig_sp = sp
    for letter in filter(!=(letter), Plots.axes_letters(orig_sp, letter))
        ax = orig_sp[Plots.get_attr_symbol(letter, :axis)]
        ax[:grid] = false  # disable the grid (overlaps with twin axis)
    end
    if orig_sp[:framestyle] === :box
        # incompatible with shared axes (see github.com/JuliaPlots/Plots.jl/issues/2894)
        orig_sp[:framestyle] = :axes
    end
    plot!(
        # plt;
        sp;
        inset = (sp[:subplot_index], bbox(0, 0, 1, 1)),
        left_margin = orig_sp[:left_margin],
        top_margin = orig_sp[:top_margin],
        right_margin = orig_sp[:right_margin],
        bottom_margin = orig_sp[:bottom_margin],
    )
    twin_sp = last(plt.subplots)
    letters = Plots.axes_letters(twin_sp, letter)
    tax, oax = map(l -> twin_sp[Plots.get_attr_symbol(l, :axis)], letters)
    tax[:grid] = false
    tax[:showaxis] = false
    tax[:ticks] = :none
    oax[:grid] = false
    oax[:mirror] = true
    twin_sp[:background_color_inside] = RGBA{Float64}(0, 0, 0, 0)
    Plots.link_axes!(sp[Plots.get_attr_symbol(letter, :axis)], tax)
    twin_sp
end



## Manual Layer Assignment
subject = "AG5";recordsession = "V1";recordsite = "11R"
siteid = join(filter!(!isempty,[subject,recordsession,recordsite]),"_")
siteresultdir = joinpath(resultroot,subject,siteid)

layer = load(joinpath(siteresultdir,"layer.jld2"),"layer")
# layer = Dict()

layer["1"] = [3100,3665]
layer["2/3A"] = [3000,3515]
layer["3B"] = [1350,3100]
layer["4A"] = [2700,2425]
layer["4B"] = [2350,2258]
layer["4Cα"] = [2250,2770]
layer["4Cβ"] = [2050,1945]
layer["5"] = [1800,1895]
layer["6A"] = [1600,1735]
layer["6B"] = [550,2070]
layer["WM"] = [0,1910]
layer["GM"] = [0,1310]

w23 = layer["2/3A"][2]-layer["4A"][2]
w3b = round(Int,w23/3)
layer["3B"] = [1350,layer["4A"][2]+w3b]

# Finalize Layer
layer = checklayer!(layer)
jldsave(joinpath(siteresultdir,"layer.jld2");layer,siteid)

# Plot layer and responses
flashlayer(siteid,siteresultdir;layer)
flashlayer(siteid,siteresultdir;layer,layerext=200,figfmt=[".png"])



function flashlayer(siteid, siteresultdir; ii='0', layer=nothing, layerext=nothing, scdepth=500, gfreq=(30, 100), figfmt=[".png"])

    test = "Flash2Color"
    testids = ["$(siteid)_$(test)_$i" for i in 0:3]

    if layer == :batch
        layerpath = joinpath(siteresultdir, "layer.jld2")
        layer = ispath(layerpath) ? load(layerpath, "layer") : nothing
    end

    ## Unit
    if !isnothing(layer)
        df = load(joinpath(siteresultdir, "unitdepthfeature.jld2"), "df")
        udf = DataFrame(df["depthfeature"],df["feature"][:])
        plotlayerunitfeature(udf, df["depth"]; layer)
        foreach(ext -> savefig(joinpath(siteresultdir, "Layer_UnitFeature$ext")), figfmt)
    end


    ## AP
    aps = load.(joinpath.(siteresultdir, testids, "ap$ii.jld2"))
    titles = repeat(["$(i["exenv"]["eye"])_$(i["exenv"]["color"])" for i in aps], inner=2)
    colors = mapreduce(i -> [RGBA(parse.(Float64, split(match(r".*Color=\[(.*)\]", k).captures[1], ", "))...) for k in keys(i["crms"])], append!, aps)
    minmaxcolors=[colors[i:i+1] for i in 1:2:length(colors)]
    depths = mapreduce(i -> [i["depths"], i["depths"]], append!, aps)

    crms = mapreduce(i -> [v for v in values(i["crms"])], append!, aps)
    times = mapreduce(i -> [i["times"], i["times"]], append!, aps)
    plotlayercondresponse(crms, times, depths; colors, titles, layer, layerext)
    foreach(ext -> savefig(joinpath(siteresultdir, "$(isnothing(layer) ? "" : "Layer_")$(test)_dRMS$ext")), figfmt)

    tps = map(i -> i["tps"], aps)
    trials = map(i -> 1:size(i, 2), tps)
    plotlayertrialresponse(tps, trials, depths; minmaxcolors, titles=titles[1:2:end], layer, layerext, color=powcm,rl="Power (μV²)",rlims=(0,4000),rscale=1e12)
    foreach(ext -> savefig(joinpath(siteresultdir, "$(isnothing(layer) ? "" : "Layer_")$(test)_tPower$ext")), figfmt)
    # σ=1(20μm): 5(80μm) diameter gaussian kernal to filter depth line power
    ftps = map(i -> imfilter(i,Kernel.gaussian((1,0)),Fill(0)), tps)
    plotlayertrialresponse(ftps, trials, depths; minmaxcolors, titles=titles[1:2:end], layer, layerext, color=powcm,rl="Power (μV²)",rlims=(0,4000),rscale=1e12)
    foreach(ext -> savefig(joinpath(siteresultdir, "$(isnothing(layer) ? "" : "Layer_")$(test)_tPowerf1$ext")), figfmt)

    tlc = map(i -> i["tlc"], aps)
    plotlayertrialresponse(tlc, trials, depths; minmaxcolors, titles=titles[1:2:end], layer, layerext, color=cohcm,rl="Coherence")
    foreach(ext -> savefig(joinpath(siteresultdir, "$(isnothing(layer) ? "" : "Layer_")$(test)_tCoherence$ext")), figfmt)

    cfps = mapreduce(i -> [v for v in values(i["cfps"])], append!, aps)
    psfreqs = mapreduce(i->[i["psfreqs"],i["psfreqs"]],append!,aps)
    plotlayercondresponse(cfps, psfreqs, depths;xl="Frequency (Hz)", colors, titles, layer, layerext,color=powcm,rl="Power (μV²)",rlims=(0,4000),rscale=1e12)
    foreach(ext -> savefig(joinpath(siteresultdir, "$(isnothing(layer) ? "" : "Layer_")$(test)_fPower$ext")), figfmt)
    # σ=1(20μm): 5(80μm) diameter gaussian kernal to filter depth line power
    fcfps = map(i -> imfilter(i,Kernel.gaussian((1,0)),Fill(0)), cfps)
    # w=21(147Hz) frequency window to filter line noises
    w = 31
    fcfps = map(i -> mapwindow(minimum,i,(1,w)), fcfps)
    plotlayercondresponse(fcfps, psfreqs, depths;xl="Frequency (Hz)", colors, titles, layer, layerext,color=powcm,rl="Power (μV²)",rlims=(0,4000),rscale=1e12)
    foreach(ext -> savefig(joinpath(siteresultdir, "$(isnothing(layer) ? "" : "Layer_")$(test)_fPowerf1w$w$ext")), figfmt)

    cflc = mapreduce(i -> [v for v in values(i["cflc"])], append!, aps)
    lcfreqs = mapreduce(i->[i["lcfreqs"],i["lcfreqs"]],append!,aps)
    plotlayercondresponse(cflc, lcfreqs, depths;xl="Frequency (Hz)", colors, titles, layer, layerext,color=cohcm,rl="Coherence")
    foreach(ext -> savefig(joinpath(siteresultdir, "$(isnothing(layer) ? "" : "Layer_")$(test)_fCoherence$ext")), figfmt)
    fcflc = map(i -> mapwindow(minimum,i,(1,w)), cflc)
    plotlayercondresponse(fcflc, lcfreqs, depths;xl="Frequency (Hz)", colors, titles, layer, layerext,color=cohcm,rl="Coherence")
    foreach(ext -> savefig(joinpath(siteresultdir, "$(isnothing(layer) ? "" : "Layer_")$(test)_fCoherencew$w$ext")), figfmt)

    jldsave(joinpath(siteresultdir,"$(test)_ap.jld2");crms,tps,ftps,cfps,fcfps,tlc,cflc,fcflc,depths,times,trials,psfreqs,lcfreqs,colors,titles,siteid)


    # downsampled Coherence and interpolated
    # for d in 2:5
    # dpc = map(i->i["pcd"][d],aps)
    # ddepths = map(i->i[1:d:end],depths)
    # plottrialresponse(dpc,trials,ddepths;colors,titles=titles[1:2:end],layer,layerext,color=cohcm)
    # foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_Coherenced$d$ext")),figfmt)
    #
    # idpc = map(i->evalgrid(Spline2D(ddepths[1],trials[1],i),depths[1],trials[1]),dpc)
    # plottrialresponse(idpc,trials,depths;colors,titles=titles[1:2:end],layer,layerext,color=cohcm)
    # foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_Coherenced$(d)i$ext")),figfmt)
    # end

    # pdc = map(i->abs.(i["pc"].-i["pbc"]),aps)
    # wpdc = map(i->abs.(i["wpc"].-i["wpbc"]),aps)
    # plottrialresponse(pdc,trials,depths;colors,titles=titles[1:2:end],layer)
    # foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_dCoherence$ext")),figfmt)
    # plottrialresponse(wpdc,trials,depths;colors,titles=titles[1:2:end],layer)
    # foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_wdCoherence$ext")),figfmt)


    ## AP+
    # aps = load.(joinpath.(siteresultdir, testids, "ap$ii+.jld2"))

    # tps = map(i -> i["tps"], aps)
    # trials = map(i -> 1:size(i, 2), tps)
    # plotlayertrialresponse(tps, trials, depths; minmaxcolors, titles=titles[1:2:end], layer, layerext, color=powcm,rl="Power (μV²)",rlims=(0,4000),rscale=1e12)
    # foreach(ext -> savefig(joinpath(siteresultdir, "$(isnothing(layer) ? "" : "Layer_")$(test)_btPower$ext")), figfmt)
    # # σ=1(20μm): 5(80μm) diameter gaussian kernal to filter depth line power
    # ftps = map(i -> imfilter(i,Kernel.gaussian((1,0)),Fill(0)), tps)
    # plotlayertrialresponse(ftps, trials, depths; minmaxcolors, titles=titles[1:2:end], layer, layerext, color=powcm,rl="Power (μV²)",rlims=(0,4000),rscale=1e12)
    # foreach(ext -> savefig(joinpath(siteresultdir, "$(isnothing(layer) ? "" : "Layer_")$(test)_btPowerf1$ext")), figfmt)

    # tlc = map(i -> i["tlc"], aps)
    # plotlayertrialresponse(tlc, trials, depths; minmaxcolors, titles=titles[1:2:end], layer, layerext, color=cohcm,rl="Coherence")
    # foreach(ext -> savefig(joinpath(siteresultdir, "$(isnothing(layer) ? "" : "Layer_")$(test)_btCoherence$ext")), figfmt)

    # cfps = mapreduce(i -> [v for v in values(i["cfps"])], append!, aps)
    # psfreqs = mapreduce(i->[i["psfreqs"],i["psfreqs"]],append!,aps)
    # plotlayercondresponse(cfps, psfreqs, depths;xl="Frequency (Hz)", colors, titles, layer, layerext,color=powcm,rl="Power (μV²)",rlims=(0,4000),rscale=1e12)
    # foreach(ext -> savefig(joinpath(siteresultdir, "$(isnothing(layer) ? "" : "Layer_")$(test)_bfPower$ext")), figfmt)
    # # σ=1(20μm): 5(80μm) diameter gaussian kernal to filter depth line power
    # fcfps = map(i -> imfilter(i,Kernel.gaussian((1,0)),Fill(0)), cfps)
    # # w=21(147Hz) frequency window to filter line noises
    # fcfps = map(i -> mapwindow(minimum,i,(1,21)), fcfps)
    # plotlayercondresponse(fcfps, psfreqs, depths;xl="Frequency (Hz)", colors, titles, layer, layerext,color=powcm,rl="Power (μV²)",rlims=(0,4000),rscale=1e12)
    # foreach(ext -> savefig(joinpath(siteresultdir, "$(isnothing(layer) ? "" : "Layer_")$(test)_bfPowerf1w21$ext")), figfmt)

    # cflc = mapreduce(i -> [v for v in values(i["cflc"])], append!, aps)
    # lcfreqs = mapreduce(i->[i["lcfreqs"],i["lcfreqs"]],append!,aps)
    # plotlayercondresponse(cflc, lcfreqs, depths;xl="Frequency (Hz)", colors, titles, layer, layerext,color=cohcm,rl="Coherence")
    # foreach(ext -> savefig(joinpath(siteresultdir, "$(isnothing(layer) ? "" : "Layer_")$(test)_bfCoherence$ext")), figfmt)
    # # w=21(147Hz) frequency window to filter line noises
    # fcflc = map(i -> mapwindow(minimum,i,(1,21)), cflc)
    # plotlayercondresponse(fcflc, lcfreqs, depths;xl="Frequency (Hz)", colors, titles, layer, layerext,color=cohcm,rl="Coherence")
    # foreach(ext -> savefig(joinpath(siteresultdir, "$(isnothing(layer) ? "" : "Layer_")$(test)_bfCoherencew21$ext")), figfmt)

    # jldsave(joinpath(siteresultdir,"$(test)_ap+.jld2");tps,ftps,cfps,fcfps,tlc,cflc,fcflc,depths,trials,psfreqs,lcfreqs,colors,titles,siteid)


    # ## Unit PSTH
    # psths = load.(joinpath.(siteresultdir, testids, "psth$ii.jld2"))
    # cpsth = mapreduce(i -> [v for v in values(i["cpsth"])], append!, psths)
    # cdpsth = mapreduce(i -> [v for v in values(i["cdpsth"])], append!, psths)
    # xs = mapreduce(i -> [i["x"], i["x"]], append!, psths)
    # ys = mapreduce(i -> [i["y"], i["y"]], append!, psths)

    # plotlayercondresponse(cpsth,xs,ys;colors,titles,layer,layerext)
    # foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_PSTH$ext")),figfmt)
    # plotlayercondresponse(cdpsth,xs,ys;colors,titles,layer,layerext)
    # foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_dPSTH$ext")),figfmt)

    # # σ=1(20μm): 5(80μm) diameter gaussian kernal to filter depth line psth
    # fcdpsth = map(i -> imfilter(i,Kernel.gaussian((1,0)),Fill(0)), cdpsth)
    # plotlayercondresponse(fcdpsth,xs,ys;colors,titles,layer,layerext)
    # foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_dPSTHf1$ext")),figfmt)

    # jldsave(joinpath(siteresultdir,"$(test)_psth.jld2");cpsth,cdpsth,fcdpsth,xs,ys,colors,titles,siteid)


    ## LF
    lfs = load.(joinpath.(siteresultdir, testids, "lf$ii.jld2"))
    clfp = mapreduce(i -> [1e6v for v in values(i["clfp"])], append!, lfs)
    cdcsd = mapreduce(i -> [v for v in values(i["ccsd"])], append!, lfs)
    times = mapreduce(i -> [i["times"], i["times"]], append!, lfs)
    depths = mapreduce(i -> [i["depths"], i["depths"]], append!, lfs)
    hy = lfs[1]["exenv"]["hy"]
    baseindex = lfs[1]["baseindex"]

    plotlayercondresponse(clfp,times,depths;colors,titles,layer,layerext)
    foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_LFP$ext")),figfmt)
    plotlayercondresponse(cdcsd,times,depths;colors,titles,layer,layerext,color=csdcm)
    foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_dCSD$ext")),figfmt)

    # σ=1(20μm): 5(80μm) diameter gaussian kernal
    f1cdcsd = map(i->imfilter(i,Kernel.gaussian((1,0)),Fill(0)),cdcsd)
    plotlayercondresponse(f1cdcsd,times,depths;colors,titles,layer,layerext,color=csdcm)
    foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_dCSDf1$ext")),figfmt)

    # σ=1.5(30μm): 9(160μm) diameter gaussian kernal
    fcdcsd = map(i -> imfilter(i, Kernel.gaussian((1.5, 0)), Fill(0)), cdcsd)
    plotlayercondresponse(fcdcsd, times, depths; colors, titles, layer, color=csdcm, layerext)
    foreach(ext -> savefig(joinpath(siteresultdir, "$(isnothing(layer) ? "" : "Layer_")$(test)_dCSDf1.5$ext")), figfmt)

    # σ=2(40μm): 9(160μm) diameter gaussian kernal
    f2cdcsd = map(i->imfilter(i,Kernel.gaussian((2,0)),Fill(0)),cdcsd)
    plotlayercondresponse(f2cdcsd,times,depths;colors,titles,layer,color=csdcm,layerext)
    foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_dCSDf2$ext")),figfmt)

    # CSD
    ccsd = map(v->csd(v,h=hy),clfp)
    plotlayercondresponse(ccsd,times,depths;colors,titles,layer,layerext,color=csdcm)
    foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_CSD$ext")),figfmt)
    # CSD filtered
    fccsds=[]
    for f in 1:0.5:2
        fccsd = map(i->imfilter(i,Kernel.gaussian((f,0)),Fill(0)),ccsd)
        push!(fccsds,fccsd)
        plotlayercondresponse(fccsd,times,depths;colors,titles,layer,layerext,color=csdcm)
        foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_CSDf$f$ext")),figfmt)
    end

    jldsave(joinpath(siteresultdir,"$(test)_lf.jld2");clfp,cdcsd,fcdcsd,ccsd,fccsd=fccsds[2],times,depths,colors,titles,siteid)

    # # downsampled and interpolated CSD
    # for d in 2:5
    #     dccsd = map(v->csd(v[1:d:end,:],h=d*hy),clfp)
    #     ddepths = map(i->i[1:d:end],depths)
    #     plotlayercondresponse(dccsd,times,ddepths;colors,titles,layer,layerext,color=csdcm)
    #     foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_CSDd$d$ext")),figfmt)

    #     idccsd = map(i->evalgrid(Spline2D(ddepths[1],times[1],i),depths[1],times[1]),dccsd)
    #     plotlayercondresponse(idccsd,times,depths;colors,titles,layer,layerext,color=csdcm)
    #     foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_CSDd$(d)i$ext")),figfmt)
    # end

    # # downsampled and interpolated dCSD
    # for d in 2:5
    #     dcdcsd = map(v->stfilter(csd(v[1:d:end,:],h=d*hy),temporaltype=:sub,ti=baseindex),clfp)
    #     ddepths = map(i->i[1:d:end],depths)
    #     plotlayercondresponse(dcdcsd,times,ddepths;colors,titles,layer,layerext,color=csdcm)
    #     foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_dCSDd$d$ext")),figfmt)

    #     idcdcsd = map(i->evalgrid(Spline2D(ddepths[1],times[1],i),depths[1],times[1]),dcdcsd)
    #     plotlayercondresponse(idcdcsd,times,depths;colors,titles,layer,layerext,color=csdcm)
    #     foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_dCSDd$(d)i$ext")),figfmt)
    # end


    # ## LF+
    # lfs = load.(joinpath.(siteresultdir, testids, "lf$ii+.jld2"))
    # depths = mapreduce(i -> [i["depths"], i["depths"]], append!, lfs)

    # tps = map(i -> i["tps"], lfs)
    # trials = map(i -> 1:size(i, 2), tps)
    # plotlayertrialresponse(tps, trials, depths; minmaxcolors, titles=titles[1:2:end], layer, layerext, color=powcm,rl="Power (μV²)",rlims=(0,100),rscale=1e12)
    # foreach(ext -> savefig(joinpath(siteresultdir, "$(isnothing(layer) ? "" : "Layer_")$(test)_tPower.lf$ext")), figfmt)

    # tlc = map(i -> i["tlc"], lfs)
    # plotlayertrialresponse(tlc, trials, depths; minmaxcolors, titles=titles[1:2:end], layer, layerext, color=cohcm,rl="Coherence")
    # foreach(ext -> savefig(joinpath(siteresultdir, "$(isnothing(layer) ? "" : "Layer_")$(test)_tCoherence.lf$ext")), figfmt)

    # cfps = mapreduce(i -> [v for v in values(i["cfps"])], append!, lfs)
    # psfreqs = mapreduce(i->[i["psfreqs"],i["psfreqs"]],append!,lfs)
    # plotlayercondresponse(cfps, psfreqs, depths;xl="Frequency (Hz)",xt=[30,60], colors, titles, layer, layerext,color=powcm,rl="Power (μV²)",rlims=(0,100),rscale=1e12)
    # foreach(ext -> savefig(joinpath(siteresultdir, "$(isnothing(layer) ? "" : "Layer_")$(test)_fPower.lf$ext")), figfmt)

    # gi = gfreq[1] .<= psfreqs[1] .<= gfreq[2]
    # psgfreqs = map(i->i[gi],psfreqs)
    # cgfps = map(i->i[:,gi],cfps)
    # plotlayercondresponse(cgfps, psgfreqs, depths;xl="Frequency (Hz)",xt=[60,90], colors, titles, layer, layerext,color=powcm,rl="Power (μV²)",rlims=(0,20),rscale=1e12)
    # foreach(ext -> savefig(joinpath(siteresultdir, "$(isnothing(layer) ? "" : "Layer_")$(test)_gfPower.lf$ext")), figfmt)

    # cflc = mapreduce(i -> [v for v in values(i["cflc"])], append!, lfs)
    # lcfreqs = mapreduce(i->[i["lcfreqs"],i["lcfreqs"]],append!,lfs)
    # plotlayercondresponse(cflc, lcfreqs, depths;xl="Frequency (Hz)",xt=[30,60], colors, titles, layer, layerext,color=cohcm,rl="Coherence")
    # foreach(ext -> savefig(joinpath(siteresultdir, "$(isnothing(layer) ? "" : "Layer_")$(test)_fCoherence.lf$ext")), figfmt)

    # gi = gfreq[1] .<= lcfreqs[1] .<= gfreq[2]
    # lcgfreqs = map(i->i[gi],lcfreqs)
    # cgflc = map(i->i[:,gi],cflc)
    # plotlayercondresponse(cgflc, lcgfreqs, depths;xl="Frequency (Hz)",xt=[60,90], colors, titles, layer, layerext,color=cohcm,rl="Coherence")
    # foreach(ext -> savefig(joinpath(siteresultdir, "$(isnothing(layer) ? "" : "Layer_")$(test)_gfCoherence.lf$ext")), figfmt)

    # jldsave(joinpath(siteresultdir,"$(test)_lf+.jld2");tps,cfps,cgfps,tlc,cflc,cgflc,depths,trials,psfreqs,psgfreqs,lcfreqs,lcgfreqs,colors,titles,siteid)

    # Superficial Cortical Gamma Power Contrast
    # if !isnothing(layer)
    #     dgp(xs,di,fi) = @views dropdims(sum(xs[di,fi,:],dims=(1,2)),dims=(1,2))
    #     di = (layer["1"][2]-scdepth).<depths[1].<layer["1"][2]
    #     fi = first(gfreq).<freqs[1].<last(gfreq)
    #     sgp = mapreduce(i->[dgp(v,di,fi) for v in values(i["cps"])],append!,lfs)
    #     sgbp = mapreduce(i->[dgp(v,di,fi) for v in values(i["cbps"])],append!,lfs)
    #     sgpc = map((p,bp)->log2.(p./bp),sgp,sgbp)
    #
    #     bar(1:length(sgpc),mean.(sgpc),yerror=sem.(sgpc),color=colors,xformatter=i->titles[Int(i)],xrot=20,
    #         xticks=1:length(sgpc),tickdir=:out,leg=false,ylabel="Superficial γ PowerContrast")
    #     foreach(ext->savefig(joinpath(siteresultdir,"Layer_$(test)_sgPowerContrast$ext")),figfmt)
    # end

end



## Batch Penetration Sites
penetration = DataFrame(XLSX.readtable(joinpath(resultroot,"penetration.xlsx"),"Sheet1"))
@showprogress "Batch Layer ... " for r in eachrow(penetration)
    flashlayer(r.siteid,joinpath(resultroot,r.Subject_ID,r.siteid),layer=nothing)

    # flashlayer(r.siteid,joinpath(resultroot,r.Subject_ID,r.siteid),layer=:batch)

    # flashlayer(r.siteid,joinpath(resultroot,r.Subject_ID,r.siteid),layer=:batch,layerext=200,figfmt=[".svg"])
end




## Collect Layers of All RecordSites
function collectlayer!(indir;layer=DataFrame(),datafile="layer.jld2")
    for (root,dirs,files) in walkdir(indir)
        if datafile in files
            try
                siteid,l = load(joinpath(root,datafile),"siteid","layer")
                df = DataFrame("siteid"=>siteid,(k=>[l[k]] for k in sort(collect(keys(l))))...)
                append!(layer,df,cols=:union)
            catch
            end
        end
    end
    return layer
end

alllayer = collectlayer!(resultroot)
jldsave(joinpath(resultroot,"alllayer.jld2");alllayer)




## Layer Normalization and Statistics
using KernelDensity
alllayer = load(joinpath(resultroot,"alllayer.jld2"),"alllayer")
excludesites = ["AG1_V1_ODR7","AG1_V1_ODR8","AG1_V1_ODL17","AG2_V1_ODL18"]
# excludesites = ["AG1_V1_ODL17"]
palllayer = filter(r->r.siteid ∉ excludesites,alllayer)


# normalized layer template
transform!(palllayer,"1"=>ByRow(i->(x->last(i).-x))=>"e2cfun") # electrode coordinates to cortical surface coordinates
transform!(palllayer,vcat.("e2cfun",names(palllayer,Not([:siteid,:e2cfun]))) .=> ByRow((f,x)->f(x)) => last)
layerwidth = combine(palllayer,names(palllayer,Not([:siteid,:e2cfun,:WM])) .=> ByRow(i->first(i)-last(i)) => identity)
layerboundary = combine(palllayer,names(palllayer,Not([:siteid,:e2cfun])) .=> ByRow(last)=>identity)

layerwidthtemplate = combine(layerwidth,names(layerwidth) .=> round ∘ mean => identity)
layerboundarytemplate = combine(layerwidthtemplate,names(layerwidthtemplate)=>ByRow((i...)->cumsum([0,i...]))=>[names(layerwidthtemplate);"WM"])
ln = names(layerwidthtemplate)
lbt = collect(first(layerboundarytemplate))
nlbt = lbt./last(lbt)
layertemplate = Dict("Out"=>[-Inf,0],(ln[i]=>[nlbt[i],nlbt[i+1]] for i in eachindex(ln))...,"WM"=>[1,Inf])
jldsave(joinpath(resultroot,"layertemplate.jld2");layertemplate,ln,lbt,nlbt)


"linearly map each layer boundary to each layertemplate boundary"
function gettcfun(x,y,isnorm=false)
    i->begin
        ii = Spline1D([first(x)-100;x;last(x)+100],[first(y)-100;y;last(y)+100],k=1,bc="extrapolate")(i)
        isnorm ? ii./last(y) : ii
    end
end
palllayer.tcfun = [gettcfun(collect(r),lbt,true) for r in eachrow(layerboundary)]


plotallsitelayer = (lwdf;dir=nothing,p="",figfmt=[".png"],sm=:width,xr=0)-> begin
    isempty(p) || (p='_'*p)
    if isnothing(sm)
        lw = lwdf
    elseif sm == :width
        lw = select!(sort!(transform(lwdf,names(lwdf,Not(1))=>ByRow(+)=>:width),:width),Not(:width))
    end

    pl = @df lw groupedbar(cols(ncol(lw):-1:2);size=(850,650),bar_position=:stack,yflip=true,palette=palette(:tab10,rev=true),lw=0,bar_width=1,
        xticks=1:nrow(lw),xlim=(0,nrow(lw)+1),xtickfontsize=6,xformatter=i->lw[Int(i),1],xr,tickor=:out,leftmargin=2Plots.mm,grid=:x,
        leg=(0.11,0.20),legendfontsize=6,ylim=(0,2250),ylabel="Cortical Depth (μm)",xlabel="Penetration",yticks=0:100:2500)
    isnothing(dir) ? pl : foreach(ext->savefig(joinpath(dir,"AllSiteLayer$p$ext")),figfmt)
end

plotlayertemplate = (;dir=nothing,figfmt=[".png"])-> begin
    p=plot(size=(250,650),palette=:tab10,yflip=true,xlim=(-0.009,0.035),ylim=(0,2250),yticks=0:100:2500,
        xticks=0:0.02:0.5,ylabel="Cortical Depth (μm)",xlabel="Probability",grid=false,tickor=:out,leftmargin=5Plots.mm)
    for l in ln
        ld = kde(layerwidth[!,l])
        plot!(p,ld.density,ld.x .+ layerboundarytemplate[1,l],label=l,lw=2)
    end
    ann = [(-0.008,mean(lbt[i:i+1]),text(ln[i],7,:gray10,:left,:vcenter)) for i in eachindex(ln)]
    hline!(p,lbt;linecolor=:gray25,legend=false,lw=1,ann)
    isnothing(dir) ? p : foreach(ext->savefig(joinpath(dir,"LayerTemplate$ext")),figfmt)
end

palllwdf = [palllayer[:,[:siteid]] layerwidth]
pallcwdf = select(palllwdf,:siteid,Not(:siteid)=>ByRow(+)=>:cw)
palllwdf = rightjoin(penetration[:,[:siteid,:pid,:od,:cofd,:d2b]],palllwdf,on=:siteid)
pallcwdf = rightjoin(penetration[:,[:siteid,:pid,:od,:cofd,:d2b]],pallcwdf,on=:siteid)

plotallsitelayer(select(palllwdf,Not([:siteid,:od,:cofd,:d2b])))
plotallsitelayer(select(palllwdf,Not([:siteid,:od,:cofd,:d2b])),dir=resultroot,figfmt=[".svg",".png"])
plotlayertemplate()
plotlayertemplate(dir=resultroot,figfmt=[".svg",".png"])


mai = map(r->startswith(r.pid,"A"),eachrow(palllwdf))
mbi = map(r->startswith(r.pid,"B"),eachrow(palllwdf))
odbi = map(r->r.od == "Both",eachrow(palllwdf))
odli = map(r->r.od == "Left",eachrow(palllwdf))
odri = map(r->r.od == "Right",eachrow(palllwdf))
odi = odli .| odri
cofdai = map(r->r.cofd ∈ ["A+","A-"],eachrow(palllwdf))
cofdsi = map(r->r.cofd ∈ ["S+","S-"],eachrow(palllwdf))
cofdlmi = map(r->!any(contains.(r.cofd,["A","S","None"])),eachrow(palllwdf))

plotallsitelayer(select(palllwdf[mai,:],Not([:siteid,:od,:cofd,:d2b])))
plotallsitelayer(select(palllwdf[odbi,:],Not([:siteid,:od,:cofd,:d2b])))
plotallsitelayer(select(palllwdf[odi,:],Not([:siteid,:od,:cofd,:d2b])))
plotallsitelayer(select(palllwdf[cofdai,:],Not([:siteid,:od,:cofd,:d2b])))
plotallsitelayer(select(palllwdf[cofdsi,:],Not([:siteid,:od,:cofd,:d2b])))
plotallsitelayer(select(palllwdf[cofdlmi,:],Not([:siteid,:od,:cofd,:d2b])))


# transform!(pallcwdf,:pid=>ByRow(i->split(i," ")[1])=>:pid,:pid=>ByRow(first)=>:mid)
# transform!(pallcwdf,:mid=>ByRow(i->i=='A' ? :royalblue : :tomato)=>:mcolor)

# @df pallcwdf scatter(:d2b,:cw,size=(750,600),leg=:inline,color=:mcolor,group=:pid,ms=5,xlabel="Distance to V1/V2 Border (mm)",ylabel="Cortical Thickness (μm)")
# foreach(ext->savefig(joinpath(resultroot,"Thickness_Border$ext")),[".png",".svg"])


# p=plot(layout=grid(1,3,widths=[0.25, 0.375, 0.375]),size=(800,550),leg=false,grid=false,tickdir=:out,ylims=(1350,2150))
# @df pallcwdf boxplot!(p[1],:mid,:cw,fillalpha=0.75, linewidth=2,ylabel="Cortical Thickness (μm)",xlabel="Monkey")
# ht = KruskalWallisTest(pallcwdf.cw[mai],pallcwdf.cw[mbi])

# @df pallcwdf[odbi .| odi,:] boxplot!(p[2],:od,:cw,fillalpha=0.75, linewidth=2,yticks=[],xlabel="Ocular Dominance")
# ht = KruskalWallisTest(pallcwdf.cw[odbi],pallcwdf.cw[odli],pallcwdf.cw[odri])

# insertcols!(pallcwdf,:cofdg=>[i ? "COFD-A" : "" for i in cofdai])
# pallcwdf.cofdg[cofdsi] .= "COFD-S"
# pallcwdf.cofdg[cofdlmi] .= "COFD-LM"
# @df pallcwdf[cofdai .| cofdsi .| cofdlmi,:] boxplot!(p[3],:cofdg,:cw,fillalpha=0.75, linewidth=2,yticks=[],xlabel="COFD Group")
# ht = KruskalWallisTest(pallcwdf.cw[cofdai],pallcwdf.cw[cofdsi],pallcwdf.cw[cofdlmi])

# foreach(ext->savefig(joinpath(resultroot,"Thickness_Group$ext")),[".png",".svg"])



## Aligned Layer Unit Features
function collectunitdepthfeature!(indir;udf=DataFrame(),datafile="unitdepthfeature.jld2")
    for (root,dirs,files) in walkdir(indir)
        if datafile in files
            d = load(joinpath(root,datafile),"df")
            df = DataFrame("siteid"=>d["siteid"],"depth"=>[d["depth"]],(d["feature"][i]=>[d["depthfeature"][:,i]] for i in 1:length(d["feature"]))...)
            append!(udf,df,cols=:union)
        end
    end
    return udf
end

alludf = collectunitdepthfeature!(resultroot)
palludf = select(alludf,[:siteid,:depth,:Density,:duration,:peaktroughratio,:uppvinv,:downpvinv],[:upspread,:downspread]=>ByRow(.-)=>:spread,
    :uppvinv=>ByRow(i->inv.(i))=>:uppv,:downpvinv=>ByRow(i->inv.(i))=>:downpv)
palludf = leftjoin(palllayer[:,[:siteid,:e2cfun,:tcfun]],palludf,on=:siteid)
transform!(palludf,[:tcfun,:e2cfun,:depth]=>ByRow((g,f,x)->g(f(x)))=>:depth)
palludf = leftjoin!(palludf,penetration[:,[:siteid,:od,:cofd,:d2b]],on=:siteid)

"Cubic Spline interpolation of depth 1D feature to layer depth template"
function getdffun(x,y)
    i->Spline1D(reverse(x),reverse(y),k=3,bc="zero")(i)
end

plotlayeralignedfeature = (ltf,f;dir=nothing,p="",color=HSL(1,0.95,0.46),xlabel="",xlims=(),figfmt=[".png"],norm=nothing,vl=[])->begin
    isempty(p) || (p='_'*p)
    isnothing(dir) || isdir(dir) || mkpath(dir)
    y = range(0,1.2,length=2000)
    ci = 0 .<= y .<= 1
    af = hcat(map((i,j)->getdffun(i,j)(y),ltf.depth,ltf[!,f])...)
    if norm==:minmax
        @views foreach(i->af[:,i] = clampscale(af[:,i],extrema(af[ci,i])...),1:size(af,2))
    elseif norm==:absmax
        @views foreach(i->af[:,i] = af[:,i]/maximum(abs.(af[:,i])),1:size(af,2))
    end
    xmin,xmax = isempty(xlims) ? extrema(af) : xlims
    afm = mean(af,dims=2)
    afse = std(af,dims=2)/sqrt(size(af,2))

    pl = plot(size=(300,650),leftmargin=4Plots.mm)
    plot!(pl,af,y;grid=false,yflip=true,color=:gray80,lw=1,ylims=(-0.01,1.21),xlims,ylabel="Normalized Cortical Depth",tickor=:out)
    isempty(vl) || vline!(pl,vl;linecolor=:gray10,legend=false,lw=1)
    plot!(pl,[afm.-afse;reverse(afm.+afse)],[y;reverse(y)]; st=:shape,lw=0,alpha=0.3,color,xlabel)
    plot!(pl,afm,y;label=false,lw=2,color)
    ann = [(xmin+0.02(xmax-xmin),mean(nlbt[i:i+1]),text(ln[i],7,:gray10,:left,:vcenter)) for i in 1:length(ln)]
    hline!(pl,nlbt;linecolor=:gray25,legend=false,lw=1,ann)

    isnothing(dir) ? pl : foreach(ext->savefig(joinpath(dir,"AlignedLayer_$f$p$ext")),figfmt)
end

unitaligneddir = joinpath(resultroot,"UnitAligned")
plotlayeralignedfeature(palludf,"Density";xlabel="Density (unit/mm³)")
plotlayeralignedfeature(palludf,"Density";xlabel="Density (unit/mm³)",dir=unitaligneddir,figfmt=[".svg",".png"])
plotlayeralignedfeature(palludf,"Density";xlabel="Normalized Density",dir=nothing,norm=:minmax)
plotlayeralignedfeature(palludf,"Density";xlabel="Normalized Density",dir=unitaligneddir,norm=:minmax,figfmt=[".svg",".png"])
plotlayeralignedfeature(palludf,"spread";xlabel="Spike Spread (μm)")
plotlayeralignedfeature(palludf,"spread";xlabel="Spike Spread (μm)",dir=unitaligneddir,figfmt=[".svg",".png"])
plotlayeralignedfeature(palludf,"spread";xlabel="Normalized Spike Spread",dir=nothing,norm=:minmax)
plotlayeralignedfeature(palludf,"spread";xlabel="Normalized Spike Spread",dir=unitaligneddir,norm=:minmax,figfmt=[".svg",".png"])
plotlayeralignedfeature(palludf,"duration";xlabel="Spike Duration (ms)",vl=[0])
plotlayeralignedfeature(palludf,"duration";xlabel="Spike Duration (ms)",dir=unitaligneddir,vl=[0],figfmt=[".svg",".png"])
plotlayeralignedfeature(palludf,"duration";xlabel="Normalized Spike Duration",dir=nothing,norm=:absmax,vl=[0])
plotlayeralignedfeature(palludf,"duration";xlabel="Normalized Spike Duration",dir=unitaligneddir,norm=:absmax,vl=[0],figfmt=[".svg",".png"])
plotlayeralignedfeature(palludf,"peaktroughratio";xlabel="Spike Peak/Trough",vl=[-1])
plotlayeralignedfeature(palludf,"peaktroughratio";xlabel="Spike Peak/Trough",dir=unitaligneddir,vl=[-1],figfmt=[".svg",".png"])

plotlayeralignedfeature(palludf,"uppvinv",xlabel="1/Up Propagation Speed (ms/μm)",xlims=(-2e-3,2e-3),dir=nothing,vl=[0])
plotlayeralignedfeature(palludf,"downpvinv",xlabel="1/Down Propagation Speed (ms/μm)",xlims=(-2e-3,2e-3),dir=nothing,vl=[0])
plotlayeralignedfeature(palludf,"uppv",xlabel="Up Propagation Speed (μm/ms)",xlims=(-3e5,3e5),dir=nothing,vl=[0])
plotlayeralignedfeature(palludf,"downpv",xlabel="Down Propagation Speed (μm/ms)",xlims=(-3e5,3e5),dir=nothing,vl=[0])


plotlayeralignedfeature(palludf[mbi,:],"Density";xlabel="Density (unit/mm³)",dir=unitaligneddir,figfmt=[".svg",".png"],p="MB")
plotlayeralignedfeature(palludf[odi,:],"Density";xlabel="Density (unit/mm³)",dir=unitaligneddir,figfmt=[".svg",".png"],p="od")
plotlayeralignedfeature(palludf[cofdlmi,:],"Density";xlabel="Density (unit/mm³)",dir=unitaligneddir,figfmt=[".svg",".png"],p="cofd-LM")
plotlayeralignedfeature(palludf[mbi,:],"spread";xlabel="Spike Spread (μm)",dir=unitaligneddir,figfmt=[".svg",".png"],p="MB")
plotlayeralignedfeature(palludf[odi,:],"spread";xlabel="Spike Spread (μm)",dir=unitaligneddir,figfmt=[".svg",".png"],p="od")
plotlayeralignedfeature(palludf[cofdlmi,:],"spread";xlabel="Spike Spread (μm)",dir=unitaligneddir,figfmt=[".svg",".png"],p="cofd-LM")
plotlayeralignedfeature(palludf[mbi,:],"duration";xlabel="Spike Duration (ms)",dir=unitaligneddir,vl=[0],figfmt=[".svg",".png"],p="MB")
plotlayeralignedfeature(palludf[odi,:],"duration";xlabel="Spike Duration (ms)",dir=unitaligneddir,vl=[0],figfmt=[".svg",".png"],p="od")
plotlayeralignedfeature(palludf[cofdlmi,:],"duration";xlabel="Spike Duration (ms)",dir=unitaligneddir,vl=[0],figfmt=[".svg",".png"],p="cofd-LM")
plotlayeralignedfeature(palludf[mbi,:],"peaktroughratio";xlabel="Spike Peak/Trough",dir=unitaligneddir,vl=[-1],figfmt=[".svg",".png"],p="MB")
plotlayeralignedfeature(palludf[odi,:],"peaktroughratio";xlabel="Spike Peak/Trough",dir=unitaligneddir,vl=[-1],figfmt=[".svg",".png"],p="od")
plotlayeralignedfeature(palludf[cofdlmi,:],"peaktroughratio";xlabel="Spike Peak/Trough",dir=unitaligneddir,vl=[-1],figfmt=[".svg",".png"],p="cofd-LM")



## Aligned Layer Flash2Color Features
function collectflashfeature!(indir;aplf=DataFrame(),test="Flash2Color")
    apfile = "$(test)_ap.jld2"
    lffile = "$(test)_lf.jld2"
    lfpfile = "$(test)_lf+.jld2"
    psthfile = "$(test)_psth.jld2"
    appfile = "$(test)_ap+.jld2"
    for (root,dirs,files) in walkdir(indir)
        if apfile in files
            ap = load(joinpath(root,apfile))
            lf = load(joinpath(root,lffile))
            lfp = load(joinpath(root,lfpfile))
            psth = load(joinpath(root,psthfile))
            app = load(joinpath(root,appfile))
            df = DataFrame("siteid"=>ap["siteid"],"colors"=>[ap["colors"]],"titles"=>[ap["titles"]],
                "aptimes"=>[ap["times"]],"appsfreqs"=>[ap["psfreqs"]],"aplcfreqs"=>[ap["lcfreqs"]],"apdepths"=>[ap["depths"]],"aptrials"=>[ap["trials"]],
                "crms"=>[ap["crms"]],"ftps"=>[ap["ftps"]],"fcfps"=>[ap["fcfps"]],"tlc"=>[ap["tlc"]],"fcflc"=>[ap["fcflc"]],

                "clfp"=>[lf["clfp"]],"ccsd"=>[lf["ccsd"]],"cdcsd"=>[lf["cdcsd"]],"fccsd"=>[lf["fccsd"]],"fcdcsd"=>[lf["fcdcsd"]],
                "lftimes"=>[lf["times"]],"lfdepths"=>[lf["depths"]],

                "lfpsgfreqs"=>[lfp["psgfreqs"]],"lflcgfreqs"=>[lfp["lcgfreqs"]],"cgfps"=>[lfp["cgfps"]],"cgflc"=>[lfp["cgflc"]],

                "fcdpsth"=>[psth["fcdpsth"]],"xs"=>[psth["xs"]],"ys"=>[psth["ys"]],

                "apbpsfreqs"=>[app["psfreqs"]],"apblcfreqs"=>[app["lcfreqs"]],
                "bftps"=>[app["ftps"]],"bfcfps"=>[app["fcfps"]],"btlc"=>[app["tlc"]],"bfcflc"=>[app["fcflc"]])
            append!(aplf,df,cols=:union)
        end
    end
    return aplf
end

allflashf = collectflashfeature!(resultroot)
pallflashf = leftjoin(palllayer[:,[:siteid,:e2cfun,:tcfun]],allflashf,on=:siteid)
transform!(pallflashf,[:tcfun,:e2cfun,:apdepths]=>ByRow((g,f,x)->g.(f.(x)))=>:apdepths,
            [:tcfun,:e2cfun,:ys]=>ByRow((g,f,x)->g.(f.(x)))=>:ys,
            [:tcfun,:e2cfun,:lfdepths]=>ByRow((g,f,x)->g.(f.(x)))=>:lfdepths)
# DKL_Y and DKL_Z colors order are different in AG1 and AG2, aligned to AG2
transform!(pallflashf,vcat.([:siteid],[:colors,:crms,:ftps,:fcfps,:tlc,:fcflc,:clfp,:ccsd,:cdcsd,:fccsd,:fcdcsd,:cgfps,:cgflc,:fcdpsth,:bftps,:bfcfps,:btlc,:bfcflc])
                            .=>ByRow((s,v)->startswith(s,"AG1") ? permute!(v,[1,2,3,4,6,5,8,7]) : v).=>last)
pallflashf = leftjoin!(pallflashf,penetration[:,[:siteid,:od,:cofd,:d2b]],on=:siteid)
minmaxcolors = [(pallflashf.colors[1][i], pallflashf.colors[1][i+1]) for i in 1:2:length(pallflashf.colors[1])]

absex(x;dims=1) = x[argmax(abs.(x);dims)]
absmaxnorm(c,d) = c/maximum(abs.(c[0 .<= d .<= 1,:]))/2
minmaxnorm(c,d) = clampscale(c,extrema(c[0 .<= d .<= 1,:])...)

"Cubic Spline interpolation of depth 2D feature to layer depth template"
function getdffun2(x,y,z)
    (i,j)->evalgrid(Spline2D(reverse(x),y,reverse(z,dims=1)),i,j)
end

plotlayeralignedcondresponse=(resps,times,depths;colors=[],xl="Time (ms)",xt=[1000,2000],rl="",rlims=(0,1),rscale=1,
                        rcolor=HSL(120,0.5,0.25),titles="",titlefontsize=7,w=140,h=600,color=:coolwarm,eachcm=false,clims=())->begin
    n = length(resps)
    if xl == "Time (ms)"
        if isempty(clims)
            lim = absmax(resps)
            clims = (-lim,lim)
        end
    else
        isempty(clims) && (clims = (absmin(resps),absmax(resps)))
    end
    p=plot(layout=(1,n),link=:y,legend=false,grid=false,size=(n*w,h))
    for i in 1:n
        yticks = i==1 ? (0:0.1:1) : false
        xticks = xl == "Time (ms)" ? (0:50:times[i][end]) : xt
        leftmargin = i==1 ? 6Plots.mm : -4Plots.mm
        xlabel = i==1 ? xl : ""
        rlabel = i==1 ? rl : ""
        ylabel = i==1 ? "Normalized Cortical Depth" : ""
        title = isempty(titles) ? titles : titles[i]

        if xl == "Time (ms)"
            if eachcm
                @views lim = abspct(resps[i][0 .<= depths[i] .<= 1,:])
                clims = (-lim,lim)
            end
        else
            @views eachcm && (clims = abspct2(resps[i][0 .<= depths[i] .<= 1,:]))
        end

        xmin,xmax = extrema(times[i])
        annx = xmin+0.02(xmax-xmin)
        ann = isempty(colors) ? [] : [(annx,1.1,Plots.text("■",15,colors[i],:left,:bottom))]

        heatmap!(p[i],times[i],depths[i],resps[i];yflip=true,color,clims,ylims=extrema(depths[i]),
        title,titlefontsize,yticks,xticks,tickor=:out,xlabel,ylabel,ann,
        leftmargin,bottommargin=3Plots.mm,topmargin=10Plots.mm)

        if xl != "Time (ms)" && rl != ""
            rm = dropdims(mean(resps[i]*rscale,dims=2),dims=2)
            rmse = dropdims(std(resps[i]*rscale,dims=2),dims=2)/sqrt(size(resps[i],2))
            rticks=[rlims[1],round(0.5(rlims[2]-rlims[1]),sigdigits=1)]
            pp=mytwiny(p[i])
            plot!(pp,[rm.-rmse;reverse(rm.+rmse)],[depths[i];reverse(depths[i])];yflip=true, st=:shape,lw=0,alpha=0.25,color=rcolor,leftmargin,xlabel=rlabel)
            plot!(pp,rm,depths[i];yflip=true,color=rcolor,leg=false,tickor=:out,xlims=rlims,xticks=rticks,ylims=extrema(depths[i]),lw=1,leftmargin)
        end

        ann = [(annx,mean(nlbt[i:i+1]),text(ln[i],7,:gray10,:left,:vcenter)) for i in 1:length(ln)]
        hline!(p[i],nlbt;linecolor=:gray25,legend=false,lw=1,ann)
    end
    p
end

plotlayeralignedtrialresponse=(resps,trials,depths;minmaxcolors=[],xl="Trial",rl="",rlims=(0,1),rscale=1,
                         rcolor=HSL(120,0.5,0.25),titles="",titlefontsize=7,w=140,h=600,color=:heat,eachcm=false,clims=())->begin
    n = length(resps)
    isempty(clims) && (clims = (absmin(resps),absmax(resps)))
    p=plot(layout=(1,n),link=:y,legend=false,grid=false,size=(n*w,h))
    for i in 1:n
        eachcm && (clims=abspct2(resps[i][0 .<= depths[i] .<= 1,:]))
        yticks = i==1 ? (0:0.1:1) : false
        xticks = range(start=1,stop=round(Int,0.75*trials[i][end]),length=2)
        leftmargin = i==1 ? 6Plots.mm : -4Plots.mm
        xlabel = i==1 ? xl : ""
        rlabel = i==1 ? rl : ""
        ylabel = i==1 ? "Normalized Cortical Depth" : ""
        title = isempty(titles) ? titles : titles[i]

        xmin,xmax = extrema(trials[i])
        annx = xmin+0.02(xmax-xmin)
        adx = 15
        ann = isempty(minmaxcolors) ? [] : [(annx,1.1,Plots.text("▮",15,minmaxcolors[i][2],:left,:bottom)),(annx+adx,1.1,Plots.text("▮",15,minmaxcolors[i][1],:left,:bottom))]

        heatmap!(p[i],trials[i],depths[i],resps[i];yflip=true,color,clims,ylims=extrema(depths[i]),
        title,titlefontsize,yticks,xticks,tickor=:out,xlabel,ylabel,ann,
        leftmargin,bottommargin=3Plots.mm)

        rm = dropdims(mean(resps[i]*rscale,dims=2),dims=2)
        rmse = dropdims(std(resps[i]*rscale,dims=2),dims=2)/sqrt(size(resps[i],2))
        rticks=[rlims[1],round(0.5(rlims[2]-rlims[1]),sigdigits=1)]
        pp=mytwiny(p[i])
        plot!(pp,[rm.-rmse;reverse(rm.+rmse)],[depths[i];reverse(depths[i])];yflip=true, st=:shape,lw=0,alpha=0.25,color=rcolor,leftmargin,xlabel=rlabel)
        plot!(pp,rm,depths[i];yflip=true,color=rcolor,leg=false,tickor=:out,xlims=rlims,xticks=rticks,ylims=extrema(depths[i]),lw=1,leftmargin)

        ann = [(annx,mean(nlbt[i:i+1]),text(ln[i],7,:gray10,:left,:vcenter)) for i in 1:length(ln)]
        hline!(p[i],nlbt;linecolor=:gray25,legend=false,lw=1,ann)
    end
    p
end

plotlayeralignedresponsefeature=(arfs,depths;colors=[],xl="Power",color=:red,titles="",titlefontsize=7,w=140,h=600,xticks = [0,0.5],xlims=(),vl=[],cm=nothing,clims=())->begin
    n = length(arfs)
    p=plot(layout=(1,n),link=:y,legend=false,grid=false,size=(n*w,h))
    for i in 1:n
        yticks = i==1 ? (0:0.1:1) : false
        leftmargin = i==1 ? 6Plots.mm : -4Plots.mm
        xlabel = i==1 ? xl : ""
        ylabel = i==1 ? "Normalized Cortical Depth" : ""
        title = isempty(titles) ? titles : titles[i]

        xmin,xmax = isempty(xlims) ? extrema(arfs[i]) : xlims
        annx = xmin+0.02(xmax-xmin)
        ann = isempty(colors) ? [] : [(annx,1.1,Plots.text("■",15,colors[i],:left,:bottom))]
        rm = mean(arfs[i],dims=2)
        rmse = std(arfs[i],dims=2)/sqrt(size(arfs[i],2))

        plot!(p[i],arfs[i],depths;grid=false,yflip=true,color=:gray80,lw=1,ylim=(-0.11,1.11),xlims,xlabel,ylabel,ann,
                title,titlefontsize,yticks,xticks,tickor=:out,leftmargin,bottommargin=3Plots.mm,topmargin=10Plots.mm)
        isempty(vl) || vline!(p[i],vl;linecolor=:gray10,legend=false,lw=1)
        plot!(p[i],[rm.-rmse;reverse(rm.+rmse)],[depths;reverse(depths)]; st=:shape,lw=0,alpha=0.3,color=isnothing(cm) ? color : :black)

        plot!(p[i],rm,depths;label=false,lw=2,color=isnothing(cm) ? color : get(cm,clampscale(rm,clims...)))

        ann = [(annx,mean(nlbt[i:i+1]),text(ln[i],7,:gray10,:left,:vcenter)) for i in 1:length(ln)]
        hline!(p[i],nlbt;linecolor=:gray25,legend=false,lw=1,ann)
    end
    p
end


function flashlayeralignedresponse(pallflashf,resultdir;p="",figfmt=[".png"])
    isempty(p) || (p='_'*p)
    isdir(resultdir) || mkpath(resultdir)
    ndepths = range(-0.1,1.1,length=2000)

    # LF clfp
    # aclfps = [cat(map(r->getdffun2(r.lfdepths[i],r.lftimes[i],r.clfp[i])(ndepths,r.lftimes[i]),eachrow(pallflashf))...,dims=3) for i in 1:8]
    # naclfps = [cat(map(r->absmaxnorm(getdffun2(r.lfdepths[i],r.lftimes[i],r.clfp[i])(ndepths,r.lftimes[i]),ndepths),eachrow(pallflashf))...,dims=3) for i in 1:8]
    # maclfp = map(i->dropdims(mean(i,dims=3),dims=3),aclfps)
    # mnaclfp = map(i->dropdims(mean(i,dims=3),dims=3),naclfps)

    # plotlayeralignedcondresponse(maclfp,pallflashf.lftimes[1],fill(ndepths,8);colors=pallflashf.colors[1])
    # foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_LFP$p$ext")),figfmt)
    # plotlayeralignedcondresponse(mnaclfp,pallflashf.lftimes[1],fill(ndepths,8);colors=pallflashf.colors[1],eachcm=true)
    # foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_nLFP$p$ext")),figfmt)

    # # LF ccsd
    # accsds = [cat(map(r->getdffun2(r.lfdepths[i],r.lftimes[i],r.ccsd[i])(ndepths,r.lftimes[i]),eachrow(pallflashf))...,dims=3) for i in 1:8]
    # naccsds = [cat(map(r->absmaxnorm(getdffun2(r.lfdepths[i],r.lftimes[i],r.ccsd[i])(ndepths,r.lftimes[i]),ndepths),eachrow(pallflashf))...,dims=3) for i in 1:8]
    # maccsd = map(i->dropdims(mean(i,dims=3),dims=3),accsds)
    # mnaccsd = map(i->dropdims(mean(i,dims=3),dims=3),naccsds)

    # plotlayeralignedcondresponse(maccsd,pallflashf.lftimes[1],fill(ndepths,8);colors=pallflashf.colors[1],color=csdcm)
    # foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_CSD$p$ext")),figfmt)
    # plotlayeralignedcondresponse(mnaccsd,pallflashf.lftimes[1],fill(ndepths,8);colors=pallflashf.colors[1],color=csdcm,eachcm=true)
    # foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_nCSD$p$ext")),figfmt)

    # # LF cdcsd
    # acdcsds = [cat(map(r->getdffun2(r.lfdepths[i],r.lftimes[i],r.cdcsd[i])(ndepths,r.lftimes[i]),eachrow(pallflashf))...,dims=3) for i in 1:8]
    # nacdcsds = [cat(map(r->absmaxnorm(getdffun2(r.lfdepths[i],r.lftimes[i],r.cdcsd[i])(ndepths,r.lftimes[i]),ndepths),eachrow(pallflashf))...,dims=3) for i in 1:8]
    # macdcsd = map(i->dropdims(mean(i,dims=3),dims=3),acdcsds)
    # mnacdcsd = map(i->dropdims(mean(i,dims=3),dims=3),nacdcsds)

    # plotlayeralignedcondresponse(macdcsd,pallflashf.lftimes[1],fill(ndepths,8);colors=pallflashf.colors[1],color=csdcm)
    # foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_dCSD$p$ext")),figfmt)
    # plotlayeralignedcondresponse(mnacdcsd,pallflashf.lftimes[1],fill(ndepths,8);colors=pallflashf.colors[1],color=csdcm,eachcm=true)
    # foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_ndCSD$p$ext")),figfmt)

    # # LF fccsd
    # afccsds = [cat(map(r->getdffun2(r.lfdepths[i],r.lftimes[i],r.fccsd[i])(ndepths,r.lftimes[i]),eachrow(pallflashf))...,dims=3) for i in 1:8]
    # nafccsds = [cat(map(r->absmaxnorm(getdffun2(r.lfdepths[i],r.lftimes[i],r.fccsd[i])(ndepths,r.lftimes[i]),ndepths),eachrow(pallflashf))...,dims=3) for i in 1:8]
    # mafccsd = map(i->dropdims(mean(i,dims=3),dims=3),afccsds)
    # mnafccsd = map(i->dropdims(mean(i,dims=3),dims=3),nafccsds)

    # plotlayeralignedcondresponse(mafccsd,pallflashf.lftimes[1],fill(ndepths,8);colors=pallflashf.colors[1],color=csdcm)
    # foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_CSDf$p$ext")),figfmt)
    # plotlayeralignedcondresponse(mnafccsd,pallflashf.lftimes[1],fill(ndepths,8);colors=pallflashf.colors[1],color=csdcm,eachcm=true)
    # foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_nCSDf$p$ext")),figfmt)

    # # LF fcdcsd
    # afcdcsds = [cat(map(r->getdffun2(r.lfdepths[i],r.lftimes[i],r.fcdcsd[i])(ndepths,r.lftimes[i]),eachrow(pallflashf))...,dims=3) for i in 1:8]
    # nafcdcsds = [cat(map(r->absmaxnorm(getdffun2(r.lfdepths[i],r.lftimes[i],r.fcdcsd[i])(ndepths,r.lftimes[i]),ndepths),eachrow(pallflashf))...,dims=3) for i in 1:8]
    # mafcdcsd = map(i->dropdims(mean(i,dims=3),dims=3),afcdcsds)
    # mnafcdcsd = map(i->dropdims(mean(i,dims=3),dims=3),nafcdcsds)

    # plotlayeralignedcondresponse(mafcdcsd,pallflashf.lftimes[1],fill(ndepths,8);colors=pallflashf.colors[1],color=csdcm)
    # foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_dCSDf$p$ext")),figfmt)
    # plotlayeralignedcondresponse(mnafcdcsd,pallflashf.lftimes[1],fill(ndepths,8);colors=pallflashf.colors[1],color=csdcm,eachcm=true)
    # foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_ndCSDf$p$ext")),figfmt)

    # # LF f3ccsd
    # afccsds = [cat(map(r->getdffun2(r.lfdepths[i],r.lftimes[i],imfilter(r.ccsd[i],Kernel.gaussian((3,0)),Fill(0)))(ndepths,r.lftimes[i]),eachrow(pallflashf))...,dims=3) for i in 1:8]
    # nafccsds = [cat(map(r->absmaxnorm(getdffun2(r.lfdepths[i],r.lftimes[i],imfilter(r.ccsd[i],Kernel.gaussian((3,0)),Fill(0)))(ndepths,r.lftimes[i]),ndepths),eachrow(pallflashf))...,dims=3) for i in 1:8]
    # mafccsd = map(i->dropdims(mean(i,dims=3),dims=3),afccsds)
    # mnafccsd = map(i->dropdims(mean(i,dims=3),dims=3),nafccsds)

    # plotlayeralignedcondresponse(mafccsd,pallflashf.lftimes[1],fill(ndepths,8);colors=pallflashf.colors[1],color=csdcm)
    # foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_CSDf3$p$ext")),figfmt)
    # plotlayeralignedcondresponse(mnafccsd,pallflashf.lftimes[1],fill(ndepths,8);colors=pallflashf.colors[1],color=csdcm,eachcm=true)
    # foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_nCSDf3$p$ext")),figfmt)

    # # LF f3cdcsd
    # afcdcsds = [cat(map(r->getdffun2(r.lfdepths[i],r.lftimes[i],imfilter(r.cdcsd[i],Kernel.gaussian((3,0)),Fill(0)))(ndepths,r.lftimes[i]),eachrow(pallflashf))...,dims=3) for i in 1:8]
    # nafcdcsds = [cat(map(r->absmaxnorm(getdffun2(r.lfdepths[i],r.lftimes[i],imfilter(r.cdcsd[i],Kernel.gaussian((3,0)),Fill(0)))(ndepths,r.lftimes[i]),ndepths),eachrow(pallflashf))...,dims=3) for i in 1:8]
    # mafcdcsd = map(i->dropdims(mean(i,dims=3),dims=3),afcdcsds)
    # mnafcdcsd = map(i->dropdims(mean(i,dims=3),dims=3),nafcdcsds)

    # plotlayeralignedcondresponse(mafcdcsd,pallflashf.lftimes[1],fill(ndepths,8);colors=pallflashf.colors[1],color=csdcm)
    # foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_dCSDf3$p$ext")),figfmt)
    # plotlayeralignedcondresponse(mnafcdcsd,pallflashf.lftimes[1],fill(ndepths,8);colors=pallflashf.colors[1],color=csdcm,eachcm=true)
    # foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_ndCSDf3$p$ext")),figfmt)


    # # LF gamma power
    # acgfpss = [cat(map(r->getdffun2(r.lfdepths[i],r.lfpsgfreqs[i],r.cgfps[i])(ndepths,r.lfpsgfreqs[i]),eachrow(pallflashf))...,dims=3) for i in 1:8]
    # nacgfpss = [cat(map(r->minmaxnorm(getdffun2(r.lfdepths[i],r.lfpsgfreqs[i],r.cgfps[i])(ndepths,r.lfpsgfreqs[i]),ndepths),eachrow(pallflashf))...,dims=3) for i in 1:8]
    # macgfps = map(i->dropdims(mean(i,dims=3),dims=3),acgfpss)
    # mnacgfps = map(i->dropdims(mean(i,dims=3),dims=3),nacgfpss)

    # plotlayeralignedcondresponse(macgfps,pallflashf.lfpsgfreqs[1],fill(ndepths,8);xl="Freqency (Hz)",xt=[60,90],colors=pallflashf.colors[1],color=powcm,rl="Power (μV²)",rlims=(0,40),rscale=1e12)
    # foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_gfPower.lf$p$ext")),figfmt)
    # plotlayeralignedcondresponse(mnacgfps,pallflashf.lfpsgfreqs[1],fill(ndepths,8);xl="Freqency (Hz)",xt=[60,90],colors=pallflashf.colors[1],color=powcm,eachcm=true)
    # foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_ngfPower.lf$p$ext")),figfmt)

    # # LF gamma coherence
    # acgflcs = [cat(map(r->getdffun2(r.lfdepths[i],r.lflcgfreqs[i],r.cgflc[i])(ndepths,r.lflcgfreqs[i]),eachrow(pallflashf))...,dims=3) for i in 1:8]
    # nacgflcs = [cat(map(r->minmaxnorm(getdffun2(r.lfdepths[i],r.lflcgfreqs[i],r.cgflc[i])(ndepths,r.lflcgfreqs[i]),ndepths),eachrow(pallflashf))...,dims=3) for i in 1:8]
    # macgflc = map(i->dropdims(mean(i,dims=3),dims=3),acgflcs)
    # mnacgflc = map(i->dropdims(mean(i,dims=3),dims=3),nacgflcs)

    # plotlayeralignedcondresponse(macgflc,pallflashf.lflcgfreqs[1],fill(ndepths,8);xl="Freqency (Hz)",xt=[60,90],colors=pallflashf.colors[1],color=cohcm,rl="Coherence")
    # foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_gfCoherence.lf$p$ext")),figfmt)
    # plotlayeralignedcondresponse(mnacgflc,pallflashf.lflcgfreqs[1],fill(ndepths,8);xl="Freqency (Hz)",xt=[60,90],colors=pallflashf.colors[1],color=cohcm,eachcm=true)
    # foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_ngfCoherence.lf$p$ext")),figfmt)


    # # Unit PSTH
    # afcdpsths = [cat(map(r->getdffun2(r.ys[i],r.xs[i],r.fcdpsth[i])(ndepths,r.xs[i]),eachrow(pallflashf))...,dims=3) for i in 1:8]
    # nafcdpsths = [cat(map(r->absmaxnorm(getdffun2(r.ys[i],r.xs[i],r.fcdpsth[i])(ndepths,r.xs[i]),ndepths),eachrow(pallflashf))...,dims=3) for i in 1:8]
    # mafcdpsth = map(i->dropdims(mean(i,dims=3),dims=3),afcdpsths)
    # mnafcdpsth = map(i->dropdims(mean(i,dims=3),dims=3),nafcdpsths)

    # plotlayeralignedcondresponse(mafcdpsth,pallflashf.xs[1],fill(ndepths,8);colors=pallflashf.colors[1])
    # foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_dPSTH$p$ext")),figfmt)
    # plotlayeralignedcondresponse(mnafcdpsth,pallflashf.xs[1],fill(ndepths,8);colors=pallflashf.colors[1],eachcm=true)
    # foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_ndPSTH$p$ext")),figfmt)


    # AP rms
    # acrmss = [cat(map(r->getdffun2(r.apdepths[i],r.aptimes[i],r.crms[i])(ndepths,r.aptimes[i]),eachrow(pallflashf))...,dims=3) for i in 1:8]
    # nacrmss = [cat(map(r->absmaxnorm(getdffun2(r.apdepths[i],r.aptimes[i],r.crms[i])(ndepths,r.aptimes[i]),ndepths),eachrow(pallflashf))...,dims=3) for i in 1:8]
    # macrms = map(i->dropdims(mean(i,dims=3),dims=3),acrmss)
    # mnacrms = map(i->dropdims(mean(i,dims=3),dims=3),nacrmss)

    # plotlayeralignedcondresponse(macrms,pallflashf.aptimes[1],fill(ndepths,8);colors=pallflashf.colors[1])
    # foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_dRMS$p$ext")),figfmt)
    # plotlayeralignedcondresponse(mnacrms,pallflashf.aptimes[1],fill(ndepths,8);colors=pallflashf.colors[1],eachcm=true)
    # foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_ndRMS$p$ext")),figfmt)

    # # AP power
    # aftpss = [cat(map(r->getdffun2(r.apdepths[i],r.aptrials[i],r.ftps[i])(ndepths,r.aptrials[i]),eachrow(pallflashf))...,dims=3) for i in 1:4]
    # naftpss = [cat(map(r->minmaxnorm(getdffun2(r.apdepths[i],r.aptrials[i],r.ftps[i])(ndepths,r.aptrials[i]),ndepths),eachrow(pallflashf))...,dims=3) for i in 1:4]
    # maftps = map(i->dropdims(mean(i,dims=3),dims=3),aftpss)
    # mnaftps = map(i->dropdims(mean(i,dims=3),dims=3),naftpss)

    # plotlayeralignedtrialresponse(maftps,pallflashf.aptrials[1],fill(ndepths,4);minmaxcolors,color=powcm,rl="Power (μV²)",rlims=(0,2000),rscale=1e12)
    # foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_tPower$p$ext")),figfmt)
    # plotlayeralignedtrialresponse(mnaftps,pallflashf.aptrials[1],fill(ndepths,4);minmaxcolors,color=powcm,eachcm=true,rl="Normalized Power")
    # foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_ntPower$p$ext")),figfmt)

    # afcfpss = [cat(map(r->getdffun2(r.apdepths[i],r.appsfreqs[i],r.fcfps[i])(ndepths,r.appsfreqs[i]),eachrow(pallflashf))...,dims=3) for i in 1:8]
    # nafcfpss = [cat(map(r->minmaxnorm(getdffun2(r.apdepths[i],r.appsfreqs[i],r.fcfps[i])(ndepths,r.appsfreqs[i]),ndepths),eachrow(pallflashf))...,dims=3) for i in 1:8]
    # mafcfps = map(i->dropdims(mean(i,dims=3),dims=3),afcfpss)
    # mnafcfps = map(i->dropdims(mean(i,dims=3),dims=3),nafcfpss)

    # plotlayeralignedcondresponse(mafcfps,pallflashf.appsfreqs[1],fill(ndepths,8);xl="Freqency (Hz)",colors=pallflashf.colors[1],color=powcm,rl="Power (μV²)",rlims=(0,2000),rscale=1e12)
    # foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_fPower$p$ext")),figfmt)
    # plotlayeralignedcondresponse(mnafcfps,pallflashf.appsfreqs[1],fill(ndepths,8);xl="Freqency (Hz)",colors=pallflashf.colors[1],color=powcm,eachcm=true)
    # foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_nfPower$p$ext")),figfmt)

    # # AP coherence
    # atlcs = [cat(map(r->getdffun2(r.apdepths[i],r.aptrials[i],r.tlc[i])(ndepths,r.aptrials[i]),eachrow(pallflashf))...,dims=3) for i in 1:4]
    # natlcs = [cat(map(r->minmaxnorm(getdffun2(r.apdepths[i],r.aptrials[i],r.tlc[i])(ndepths,r.aptrials[i]),ndepths),eachrow(pallflashf))...,dims=3) for i in 1:4]
    # matlc = map(i->dropdims(mean(i,dims=3),dims=3),atlcs)
    # mnatlc = map(i->dropdims(mean(i,dims=3),dims=3),natlcs)

    # plotlayeralignedtrialresponse(matlc,pallflashf.aptrials[1],fill(ndepths,4);minmaxcolors,color=cohcm,rl="Coherence")
    # foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_tCoherence$p$ext")),figfmt)
    # plotlayeralignedtrialresponse(mnatlc,pallflashf.aptrials[1],fill(ndepths,4);minmaxcolors,color=cohcm,eachcm=true,rl="Normalized Coherence")
    # foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_ntCoherence$p$ext")),figfmt)

    # afcflcs = [cat(map(r->getdffun2(r.apdepths[i],r.aplcfreqs[i],r.fcflc[i])(ndepths,r.aplcfreqs[i]),eachrow(pallflashf))...,dims=3) for i in 1:8]
    # nafcflcs = [cat(map(r->minmaxnorm(getdffun2(r.apdepths[i],r.aplcfreqs[i],r.fcflc[i])(ndepths,r.aplcfreqs[i]),ndepths),eachrow(pallflashf))...,dims=3) for i in 1:8]
    # mafcflc = map(i->dropdims(mean(i,dims=3),dims=3),afcflcs)
    # mnafcflc = map(i->dropdims(mean(i,dims=3),dims=3),nafcflcs)

    # plotlayeralignedcondresponse(mafcflc,pallflashf.aplcfreqs[1],fill(ndepths,8);xl="Freqency (Hz)",colors=pallflashf.colors[1],color=cohcm,rl="Coherence")
    # foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_fCoherence$p$ext")),figfmt)
    # plotlayeralignedcondresponse(mnafcflc,pallflashf.aplcfreqs[1],fill(ndepths,8);xl="Freqency (Hz)",colors=pallflashf.colors[1],color=cohcm,eachcm=true)
    # foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_nfCoherence$p$ext")),figfmt)


    # AP baseline power
    abftpss = [cat(map(r->getdffun2(r.apdepths[i],r.aptrials[i],r.bftps[i])(ndepths,r.aptrials[i]),eachrow(pallflashf))...,dims=3) for i in 1:4]
    nabftpss = [cat(map(r->minmaxnorm(getdffun2(r.apdepths[i],r.aptrials[i],r.bftps[i])(ndepths,r.aptrials[i]),ndepths),eachrow(pallflashf))...,dims=3) for i in 1:4]
    mabftps = map(i->dropdims(mean(i,dims=3),dims=3),abftpss)
    mnabftps = map(i->dropdims(mean(i,dims=3),dims=3),nabftpss)

    plotlayeralignedtrialresponse(mabftps,pallflashf.aptrials[1],fill(ndepths,4);minmaxcolors,color=powcm,rl="Power (μV²)",rlims=(0,2000),rscale=1e12)
    foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_btPower$p$ext")),figfmt)
    plotlayeralignedtrialresponse(mnabftps,pallflashf.aptrials[1],fill(ndepths,4);minmaxcolors,color=powcm,eachcm=true,rl="Normalized Power")
    foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_nbtPower$p$ext")),figfmt)

    abfcfpss = [cat(map(r->getdffun2(r.apdepths[i],r.apbpsfreqs[i],r.bfcfps[i])(ndepths,r.apbpsfreqs[i]),eachrow(pallflashf))...,dims=3) for i in 1:8]
    nabfcfpss = [cat(map(r->minmaxnorm(getdffun2(r.apdepths[i],r.apbpsfreqs[i],r.bfcfps[i])(ndepths,r.apbpsfreqs[i]),ndepths),eachrow(pallflashf))...,dims=3) for i in 1:8]
    mabfcfps = map(i->dropdims(mean(i,dims=3),dims=3),abfcfpss)
    mnabfcfps = map(i->dropdims(mean(i,dims=3),dims=3),nabfcfpss)

    plotlayeralignedcondresponse(mabfcfps,pallflashf.apbpsfreqs[1],fill(ndepths,8);xl="Freqency (Hz)",colors=pallflashf.colors[1],color=powcm,rl="Power (μV²)",rlims=(0,2000),rscale=1e12)
    foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_bfPower$p$ext")),figfmt)
    plotlayeralignedcondresponse(mnabfcfps,pallflashf.apbpsfreqs[1],fill(ndepths,8);xl="Freqency (Hz)",colors=pallflashf.colors[1],color=powcm,eachcm=true)
    foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_nbfPower$p$ext")),figfmt)

    # AP baseline coherence
    abtlcs = [cat(map(r->getdffun2(r.apdepths[i],r.aptrials[i],r.btlc[i])(ndepths,r.aptrials[i]),eachrow(pallflashf))...,dims=3) for i in 1:4]
    nabtlcs = [cat(map(r->minmaxnorm(getdffun2(r.apdepths[i],r.aptrials[i],r.btlc[i])(ndepths,r.aptrials[i]),ndepths),eachrow(pallflashf))...,dims=3) for i in 1:4]
    mabtlc = map(i->dropdims(mean(i,dims=3),dims=3),abtlcs)
    mnabtlc = map(i->dropdims(mean(i,dims=3),dims=3),nabtlcs)

    plotlayeralignedtrialresponse(mabtlc,pallflashf.aptrials[1],fill(ndepths,4);minmaxcolors,color=cohcm,rl="Coherence")
    foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_btCoherence$p$ext")),figfmt)
    plotlayeralignedtrialresponse(mnabtlc,pallflashf.aptrials[1],fill(ndepths,4);minmaxcolors,color=cohcm,eachcm=true,rl="Normalized Coherence")
    foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_nbtCoherence$p$ext")),figfmt)

    abfcflcs = [cat(map(r->getdffun2(r.apdepths[i],r.apblcfreqs[i],r.bfcflc[i])(ndepths,r.apblcfreqs[i]),eachrow(pallflashf))...,dims=3) for i in 1:8]
    nabfcflcs = [cat(map(r->minmaxnorm(getdffun2(r.apdepths[i],r.apblcfreqs[i],r.bfcflc[i])(ndepths,r.apblcfreqs[i]),ndepths),eachrow(pallflashf))...,dims=3) for i in 1:8]
    mabfcflc = map(i->dropdims(mean(i,dims=3),dims=3),abfcflcs)
    mnabfcflc = map(i->dropdims(mean(i,dims=3),dims=3),nabfcflcs)

    plotlayeralignedcondresponse(mabfcflc,pallflashf.apblcfreqs[1],fill(ndepths,8);xl="Freqency (Hz)",colors=pallflashf.colors[1],color=cohcm,rl="Coherence")
    foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_bfCoherence$p$ext")),figfmt)
    plotlayeralignedcondresponse(mnabfcflc,pallflashf.apblcfreqs[1],fill(ndepths,8);xl="Freqency (Hz)",colors=pallflashf.colors[1],color=cohcm,eachcm=true)
    foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_nbfCoherence$p$ext")),figfmt)
end

flashaligneddir = joinpath(resultroot,"FlashAligned")
flashlayeralignedresponse(pallflashf,flashaligneddir,figfmt=[".png",".svg"])
flashlayeralignedresponse(pallflashf[mai,:],flashaligneddir,p="MA",figfmt=[".png",".svg"])
flashlayeralignedresponse(pallflashf[mbi,:],flashaligneddir,p="MB",figfmt=[".png",".svg"])
flashlayeralignedresponse(pallflashf[odbi,:],flashaligneddir,p="od-B",figfmt=[".png",".svg"])
flashlayeralignedresponse(pallflashf[odi,:],flashaligneddir,p="od",figfmt=[".png",".svg"])
flashlayeralignedresponse(pallflashf[cofdai,:],flashaligneddir,p="cofd-A",figfmt=[".png",".svg"])
flashlayeralignedresponse(pallflashf[cofdsi,:],flashaligneddir,p="cofd-S",figfmt=[".png",".svg"])
flashlayeralignedresponse(pallflashf[cofdlmi,:],flashaligneddir,p="cofd-LM",figfmt=[".png",".svg"])


function flashlayeralignedresponsefeature(pallflashf,resultdir;p="",figfmt=[".png"],win=(30,100))
    isempty(p) || (p='_'*p)
    isdir(resultdir) || mkpath(resultdir)
    ndepths = range(-0.1,1.1,length=2000)

    # LF fcdcsd
    # lfwini = win[1] .<= pallflashf.lftimes[1][1] .<= win[2]
    # afcrmsdcsds = [cat(map(r->minmaxnorm(sqrt.(mean(getdffun2(r.lfdepths[i],r.lftimes[i],r.fcdcsd[i])(ndepths,r.lftimes[i])[:,lfwini].^2,dims=2)),ndepths),eachrow(pallflashf))...,dims=2) for i in 1:8]
    # plotlayeralignedresponsefeature(afcrmsdcsds,ndepths;xl="Normalized ΔCSD Amplitude",colors=pallflashf.colors[1],color=HSL(358,0.7,0.3))
    # foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_nrmsdCSD$p$ext")),figfmt)

    # afcmdcsds = [cat(map(r->absmaxnorm(mean(getdffun2(r.lfdepths[i],r.lftimes[i],r.fcdcsd[i])(ndepths,r.lftimes[i])[:,lfwini],dims=2),ndepths),eachrow(pallflashf))...,dims=2) for i in 1:8]
    # plotlayeralignedresponsefeature(afcmdcsds,ndepths;xl="Normalized ΔCSD",colors=pallflashf.colors[1],color=HSL(358,0.7,0.3),xlims=(-0.52,0.52),vl=[0],cm=csdcm,clims=(-0.4,0.4))
    # foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_nmdCSD$p$ext")),figfmt)

    # afcexdcsds = [cat(map(r->absmaxnorm(absex(getdffun2(r.lfdepths[i],r.lftimes[i],r.fcdcsd[i])(ndepths,r.lftimes[i])[:,lfwini],dims=2),ndepths),eachrow(pallflashf))...,dims=2) for i in 1:8]
    # plotlayeralignedresponsefeature(afcexdcsds,ndepths;xl="Normalized ΔCSD Extrema",colors=pallflashf.colors[1],color=HSL(358,0.7,0.3),xlims=(-0.52,0.52),vl=[0],cm=csdcm,clims=(-0.4,0.4))
    # foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_nexdCSD$p$ext")),figfmt)

    # LF fcdcsd delay (30ms window)
    afcddcsds = [cat(map(r->r.lftimes[i][getindex.(argmax(mapwindow(rms,getdffun2(r.lfdepths[i],r.lftimes[i],r.fcdcsd[i])(ndepths,r.lftimes[i]),(1,31),border="symmetric"),dims=2),2)],eachrow(pallflashf))...,dims=2) for i in 1:8]
    plotlayeralignedresponsefeature(afcddcsds,ndepths;xl="ΔCSD Delay",colors=pallflashf.colors[1],color=HSL(358,0.7,0.3),xlims=(-0.1,150.1),xticks=[0,50,100])
    foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_ddCSD$p$ext")),figfmt)

    # # LF fccsd
    # afcrmscsds = [cat(map(r->minmaxnorm(sqrt.(mean(getdffun2(r.lfdepths[i],r.lftimes[i],r.fccsd[i])(ndepths,r.lftimes[i])[:,lfwini].^2,dims=2)),ndepths),eachrow(pallflashf))...,dims=2) for i in 1:8]
    # plotlayeralignedresponsefeature(afcrmscsds,ndepths;xl="Normalized CSD Amplitude",colors=pallflashf.colors[1],color=HSL(358,0.7,0.3))
    # foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_nrmsCSD$p$ext")),figfmt)

    # afcmcsds = [cat(map(r->absmaxnorm(mean(getdffun2(r.lfdepths[i],r.lftimes[i],r.fccsd[i])(ndepths,r.lftimes[i])[:,lfwini],dims=2),ndepths),eachrow(pallflashf))...,dims=2) for i in 1:8]
    # plotlayeralignedresponsefeature(afcmcsds,ndepths;xl="Normalized CSD",colors=pallflashf.colors[1],color=HSL(358,0.7,0.3),xlims=(-0.52,0.52),vl=[0],cm=csdcm,clims=(-0.4,0.4))
    # foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_nmCSD$p$ext")),figfmt)

    # afcexcsds = [cat(map(r->absmaxnorm(absex(getdffun2(r.lfdepths[i],r.lftimes[i],r.fccsd[i])(ndepths,r.lftimes[i])[:,lfwini],dims=2),ndepths),eachrow(pallflashf))...,dims=2) for i in 1:8]
    # plotlayeralignedresponsefeature(afcexcsds,ndepths;xl="Normalized CSD Extrema",colors=pallflashf.colors[1],color=HSL(358,0.7,0.3),xlims=(-0.52,0.52),vl=[0],cm=csdcm,clims=(-0.4,0.4))
    # foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_nexCSD$p$ext")),figfmt)


    # LF gamma power
    # acgfpss = [cat(map(r->minmaxnorm(mean(getdffun2(r.lfdepths[i],r.lfpsgfreqs[i],r.cgfps[i])(ndepths,r.lfpsgfreqs[i]),dims=2),ndepths),eachrow(pallflashf))...,dims=2) for i in 1:8]
    # plotlayeralignedresponsefeature(acgfpss,ndepths;xl="Normalized Power",colors=pallflashf.colors[1],color=HSL(0,1,0.32))
    # foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_nmgPower.lf$p$ext")),figfmt)

    # # LF gamma coherence
    # acgflcs = [cat(map(r->minmaxnorm(mean(getdffun2(r.lfdepths[i],r.lflcgfreqs[i],r.cgflc[i])(ndepths,r.lflcgfreqs[i]),dims=2),ndepths),eachrow(pallflashf))...,dims=2) for i in 1:8]
    # plotlayeralignedresponsefeature(acgflcs,ndepths;xl="Normalized Coherence",colors=pallflashf.colors[1],color=HSL(233,1,0.38))
    # foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_nmgCoherence.lf$p$ext")),figfmt)


    # Unit PSTH
    # xwini = win[1] .<= pallflashf.xs[1][1] .<= win[2]
    # acrmsdpsth = [cat(map(r->minmaxnorm(sqrt.(mean((getdffun2(r.ys[i],r.xs[i],r.fcdpsth[i])(ndepths,r.xs[i])[:,xwini]).^2,dims=2)),ndepths),eachrow(pallflashf))...,dims=2) for i in 1:8]
    # plotlayeralignedresponsefeature(acrmsdpsth,ndepths;xl="Normalized ΔPSTH Amplitude",colors=pallflashf.colors[1],color=HSL(348,1,0.36))
    # foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_nrmsdPSTH$p$ext")),figfmt)

    # acmdpsth = [cat(map(r->absmaxnorm(mean(getdffun2(r.ys[i],r.xs[i],r.fcdpsth[i])(ndepths,r.xs[i])[:,xwini],dims=2),ndepths),eachrow(pallflashf))...,dims=2) for i in 1:8]
    # plotlayeralignedresponsefeature(acmdpsth,ndepths;xl="Normalized ΔPSTH",colors=pallflashf.colors[1],color=HSL(348,1,0.36),xlims=(-0.52,0.52),vl=[0],cm=cgrad(:coolwarm),clims=(-0.3,0.3))
    # foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_nmdPSTH$p$ext")),figfmt)

    # acexdpsth = [cat(map(r->absmaxnorm(absex(getdffun2(r.ys[i],r.xs[i],r.fcdpsth[i])(ndepths,r.xs[i])[:,xwini],dims=2),ndepths),eachrow(pallflashf))...,dims=2) for i in 1:8]
    # plotlayeralignedresponsefeature(acexdpsth,ndepths;xl="Normalized ΔPSTH Extrema",colors=pallflashf.colors[1],color=HSL(348,1,0.36),xlims=(-0.52,0.52),vl=[0],cm=cgrad(:coolwarm),clims=(-0.3,0.3))
    # foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_nexdPSTH$p$ext")),figfmt)


    # AP Prc
    # apwini = win[1] .<= pallflashf.aptimes[1][1] .<= win[2]
    # acrmsprcs = [cat(map(r->minmaxnorm(sqrt.(mean((getdffun2(r.apdepths[i],r.aptimes[i],r.crms[i])(ndepths,r.aptimes[i])[:,apwini]).^2,dims=2)),ndepths),eachrow(pallflashf))...,dims=2) for i in 1:8]
    # plotlayeralignedresponsefeature(acrmsprcs,ndepths;xl="Normalized Prc Amplitude",colors=pallflashf.colors[1],color=HSL(348,1,0.36))
    # foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_nrmsPrc$p$ext")),figfmt)

    # acmprcs = [cat(map(r->absmaxnorm(mean(getdffun2(r.apdepths[i],r.aptimes[i],r.crms[i])(ndepths,r.aptimes[i])[:,apwini],dims=2),ndepths),eachrow(pallflashf))...,dims=2) for i in 1:8]
    # plotlayeralignedresponsefeature(acmprcs,ndepths;xl="Normalized Prc",colors=pallflashf.colors[1],color=HSL(348,1,0.36),xlims=(-0.52,0.52),vl=[0],cm=cgrad(:coolwarm),clims=(-0.3,0.3))
    # foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_nmPrc$p$ext")),figfmt)

    # acexprcs = [cat(map(r->absmaxnorm(absex(getdffun2(r.apdepths[i],r.aptimes[i],r.crms[i])(ndepths,r.aptimes[i])[:,apwini],dims=2),ndepths),eachrow(pallflashf))...,dims=2) for i in 1:8]
    # plotlayeralignedresponsefeature(acexprcs,ndepths;xl="Normalized Prc Extrema",colors=pallflashf.colors[1],color=HSL(348,1,0.36),xlims=(-0.52,0.52),vl=[0],cm=cgrad(:coolwarm),clims=(-0.3,0.3))
    # foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_nexPrc$p$ext")),figfmt)

    # AP Prc delay (30ms window)
    acdprcs = [cat(map(r->r.aptimes[i][getindex.(argmax(mapwindow(rms,getdffun2(r.apdepths[i],r.aptimes[i],r.crms[i])(ndepths,r.aptimes[i]),(1,301),border="symmetric"),dims=2),2)],eachrow(pallflashf))...,dims=2) for i in 1:8]
    plotlayeralignedresponsefeature(acdprcs,ndepths;xl="Prc Delay",colors=pallflashf.colors[1],color=HSL(348,1,0.36),xlims=(-0.1,150.1),xticks=[0,50,100])
    foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_dPrc$p$ext")),figfmt)

    # # AP power
    # afcfpss = [cat(map(r->minmaxnorm(mean(getdffun2(r.apdepths[i],r.appsfreqs[i],r.fcfps[i])(ndepths,r.appsfreqs[i]),dims=2),ndepths),eachrow(pallflashf))...,dims=2) for i in 1:8]
    # plotlayeralignedresponsefeature(afcfpss,ndepths;xl="Normalized Power",colors=pallflashf.colors[1],color=HSL(0,1,0.32))
    # foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_nmPower$p$ext")),figfmt)

    # # AP coherence
    # afcflcs = [cat(map(r->minmaxnorm(mean(getdffun2(r.apdepths[i],r.aplcfreqs[i],r.fcflc[i])(ndepths,r.aplcfreqs[i]),dims=2),ndepths),eachrow(pallflashf))...,dims=2) for i in 1:8]
    # plotlayeralignedresponsefeature(afcflcs,ndepths;xl="Normalized Coherence",colors=pallflashf.colors[1],color=HSL(233,1,0.38))
    # foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_nmCoherence$p$ext")),figfmt)

    # # AP baseline power
    # abfcfpss = [cat(map(r->minmaxnorm(mean(getdffun2(r.apdepths[i],r.apbpsfreqs[i],r.bfcfps[i])(ndepths,r.apbpsfreqs[i]),dims=2),ndepths),eachrow(pallflashf))...,dims=2) for i in 1:8]
    # plotlayeralignedresponsefeature(abfcfpss,ndepths;xl="Normalized Power",colors=pallflashf.colors[1],color=HSL(0,1,0.32))
    # foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_nmbPower$p$ext")),figfmt)

    # # AP baseline coherence
    # abfcflcs = [cat(map(r->minmaxnorm(mean(getdffun2(r.apdepths[i],r.apblcfreqs[i],r.bfcflc[i])(ndepths,r.apblcfreqs[i]),dims=2),ndepths),eachrow(pallflashf))...,dims=2) for i in 1:8]
    # plotlayeralignedresponsefeature(abfcflcs,ndepths;xl="Normalized Coherence",colors=pallflashf.colors[1],color=HSL(233,1,0.38))
    # foreach(ext->savefig(joinpath(resultdir,"AlignedLayer_nmbCoherence$p$ext")),figfmt)

    # Dict(:rmsdcsd=>afcrmsdcsds,:mdcsd=>afcmdcsds,:exdcsd=>afcexdcsds,:ddcsd=>afcddcsds,:rmscsd=>afcrmscsds,:mcsd=>afcmcsds,:excsd=>afcexcsds,
    #     :rmsdpsth=>acrmsdpsth,:mdpsth=>acmdpsth,:exdpsth=>acexdpsth,
    #     :rmsprc=>acrmsprcs,:mprc=>acmprcs,:exprc=>acexprcs,:dprc=>acdprcs,:mpower=>afcfpss,:mlc=>afcflcs,:ndepths=>ndepths)
end

flarf=flashlayeralignedresponsefeature(pallflashf,flashaligneddir,figfmt=[".png",".svg"])
flashlayeralignedresponsefeature(pallflashf[mai,:],flashaligneddir,p="MA",figfmt=[".png",".svg"])
flashlayeralignedresponsefeature(pallflashf[mbi,:],flashaligneddir,p="MB",figfmt=[".png",".svg"])
flashlayeralignedresponsefeature(pallflashf[odbi,:],flashaligneddir,p="od-B",figfmt=[".png",".svg"])
flashlayeralignedresponsefeature(pallflashf[odi,:],flashaligneddir,p="od",figfmt=[".png",".svg"])
flashlayeralignedresponsefeature(pallflashf[cofdai,:],flashaligneddir,p="cofd-A",figfmt=[".png",".svg"])
flashlayeralignedresponsefeature(pallflashf[cofdsi,:],flashaligneddir,p="cofd-S",figfmt=[".png",".svg"])
flashlayeralignedresponsefeature(pallflashf[cofdlmi,:],flashaligneddir,p="cofd-LM",figfmt=[".png",".svg"])



function layerresponsedf(lrf,ndepths;sti=["White~","Black~","White","Black","Red","Green","Purple","Lime"])
    lrdf = DataFrame(lr=map(c->vcat(map(i->mean(c[nlbt[i] .<= ndepths .< nlbt[i+1],:],dims=1),eachindex(ln))...),lrf), sti=sti)
    transform!(lrdf,:lr=>ByRow(i->([i[:,p] for p in 1:size(i,2)],pallflashf.siteid))=>[:lr,:siteid])
    lrdf = flatten(lrdf,[:lr,:siteid])
    lrdf[!,:layer].=[ln]
    lrdf=flatten(lrdf,[:lr,:layer])
    leftjoin!(lrdf,penetration[:,[:siteid,:pid,:od,:cofd,:d2b]],on=:siteid)
    transform!(lrdf,:pid=>ByRow(first)=>:mid,:siteid=>ByRow(i->begin
    if i in pallflashf.siteid[cofdai]
        "COFD-A"
    elseif i in pallflashf.siteid[cofdsi]
        "COFD-S"
    elseif i in pallflashf.siteid[cofdlmi]
        "COFD-LM"
    else
        missing
    end
    end)=>:cofdg,:od=>ByRow(i->i in ["Left-","Right-"] ? missing : i)=>:odg,:sti=>ByRow(i->begin
        if i in ["White~","Black~"]
            "A~"
        elseif i in ["White","Black"]
            "A"
        elseif i in ["Green","Red"]
            "LM"
        else
            "S"
        end
    end)=>:stig)
end

lmdcsd = layerresponsedf(flarf[:mdcsd],flarf[:ndepths])
lmpower = layerresponsedf(flarf[:mpower],flarf[:ndepths])
lmprc = layerresponsedf(flarf[:mprc],flarf[:ndepths])


using VegaLite

# pl = lmdcsd |> @vlplot(
#     mark={:boxplot, extent="min-max"},
#     x={:mid,axis={title="Monkey"}},
#     y={:lr, axis={title="Average ΔCSD"}},
#     row="layer",
#     column={"sti",sort=["White~","Black~","White","Black","Red","Green","Purple","Lime"]}
# )
# foreach(ext->save(joinpath(flashaligneddir,"lmdCSD_mid$ext"),pl),[".png",".svg"])

# pl = dropmissing(lmdcsd,:odg) |> @vlplot(
#     mark={:boxplot, extent="min-max"},
#     x={:odg,axis={title="OD Group"}},
#     y={:lr, axis={title="Average ΔCSD"}},
#     row="layer",
#     column={"sti",sort=["White~","Black~","White","Black","Red","Green","Purple","Lime"]}
# )
# foreach(ext->save(joinpath(flashaligneddir,"lmdCSD_odg$ext"),pl),[".png",".svg"])

# pl = dropmissing(lmdcsd,:cofdg) |> @vlplot(
#     mark={:boxplot, extent="min-max"},
#     x={:cofdg,axis={title="COFD Group"}},
#     y={:lr, axis={title="Average ΔCSD"}},
#     row="layer",
#     column={"sti",sort=["White~","Black~","White","Black","Red","Green","Purple","Lime"]}
# )
# foreach(ext->save(joinpath(flashaligneddir,"lmdCSD_cofdg$ext"),pl),[".png",".svg"])


lcofddcsd=dropmissing(filter(r->r.stig in ["A","LM","S"],lmdcsd),:cofdg)

pl = lcofddcsd |> @vlplot(width=300,height=180,
    mark={:boxplot, extent="min-max"},
    x={:cofdg,axis={title="COFD Group",titleFontSize=24,titleFontWeight="normal",labelAngle=0,labelFontSize=18}},
    y={:lr, axis={title="Average ΔCSD",titleFontSize=20,titleFontWeight="normal",grid=false,labelFontSize=14}},
    row={"layer",header={labelAngle=0,labelFontSize=26,labelAlign="left",labelFontWeight="bold",title=""}},
    column={"stig",sort=["A~","A","LM","S"],header={labelFontSize=26,labelFontWeight="bold",title=""}}
)
foreach(ext->save(joinpath(flashaligneddir,"ldCSD_cofdg_stig$ext"),pl),[".png",".svg"])

t=filter(r->r.layer in ["1"] && r.stig in ["A"],lcofddcsd)
ht = MannWhitneyUTest(filter(r->r.cofdg=="COFD-A",t).lr,filter(r->r.cofdg=="COFD-LM",t).lr)
ht = MannWhitneyUTest(filter(r->r.cofdg=="COFD-A",t).lr,filter(r->r.cofdg=="COFD-S",t).lr)
ht = MannWhitneyUTest(filter(r->r.cofdg=="COFD-S",t).lr,filter(r->r.cofdg=="COFD-LM",t).lr)


lprc=dropmissing(filter(r->r.stig in ["A","LM","S"],lmpower),:cofdg)

pl = lprc |> @vlplot(width=300,height=180,
    mark={:boxplot, extent="min-max"},
    x={:stig,axis={title="Stimuli",titleFontSize=24,titleFontWeight="normal",labelAngle=0,labelFontSize=18}},
    y={:lr, axis={title="Average Prc",titleFontSize=20,titleFontWeight="normal",grid=false,labelFontSize=14}},
    row={"layer",header={labelAngle=0,labelFontSize=26,labelAlign="left",labelFontWeight="bold",title=""}}
    # column={"stig",sort=["A~","A","LM","S"],header={labelFontSize=26,labelFontWeight="bold",title=""}}
)
foreach(ext->save(joinpath(flashaligneddir,"lprc_stig$ext"),pl),[".png",".svg"])

t=filter(r->r.layer in ["4A"],lprc)
ht = MannWhitneyUTest(filter(r->r.stig=="A",t).lr,filter(r->r.stig=="LM",t).lr)
ht = MannWhitneyUTest(filter(r->r.stig=="A",t).lr,filter(r->r.stig=="S",t).lr)
ht = MannWhitneyUTest(filter(r->r.stig=="LM",t).lr,filter(r->r.stig=="S",t).lr)

using GLM

lmr = lm(@formula(lr~cofdg+stig+cofdg*stig),t)
lmr = lm(@formula(lr~cofdg),t)
