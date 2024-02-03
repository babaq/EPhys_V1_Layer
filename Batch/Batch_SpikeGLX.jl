using Images,DSP

function process_flash_spikeglx(files,param;uuid="",log=nothing,plot=true)
    dataset = prepare(joinpath(param[:dataexportroot],files),spikesorter=param[:spikesorter])

    ex = dataset["ex"];subject = ex["Subject_ID"];recordsession = ex["RecordSession"];recordsite = ex["RecordSite"]
    siteid = join(filter(!isempty,[subject,recordsession,recordsite]),"_")
    datadir = joinpath(param[:dataroot],subject,siteid)
    siteresultdir = joinpath(param[:resultroot],subject,siteid)
    # testid = join(filter(!isempty,[siteid,ex["TestID"]]),"_")
    testid = splitdir(splitext(ex["source"])[1])[2]
    resultdir = joinpath(siteresultdir,testid)
    isdir(resultdir) || mkpath(resultdir)

    spike = dataset["spike"];unitspike = spike["unitspike"];unitid = spike["unitid"]
    unitgood=spike["unitgood"];unitposition=spike["unitposition"];unitsync=spike["unitsync"]
    layer = haskey(param,:layer) ? param[:layer] : nothing
    figfmt = haskey(param,:figfmt) ? param[:figfmt] : [".png"]
    jldsave(joinpath(resultdir,"spike.jld2");spike,siteid)
    haskey(param,:onlyspike) && param[:onlyspike] && return

    # Condition Tests
    envparam = ex["EnvParam"];exparam = ex["Param"];preicidur = ex["PreICI"];conddur = ex["CondDur"];suficidur = ex["SufICI"]
    condon = ex["CondTest"]["CondOn"]
    condoff = ex["CondTest"]["CondOff"]
    condidx = ex["CondTest"]["CondIndex"]
    minconddur=minimum(condoff.-condon)
    exenv=Dict()
    exenv["ID"] = ex["ID"]
    exenv["conddur"] = conddur
    if haskey(ex,"Eye")
        exenv["eye"] = ex["Eye"]
    else
        exenv["eye"] = ex["Log"]
    end
    exenv["color"] = "$(exparam["ColorSpace"])_$(exparam["Color"])"

    # Prepare Conditions
    ctc = condtestcond(ex["CondTestCond"])
    cond = condin(ctc)
    factors = finalfactor(ctc)

    for ii in dataset["imecindex"]
        # Prepare AP
        # d = "ap$ii"
        # meta = dataset[d]["meta"]
        # binfile = meta["fileName"]
        # nsavedch = meta["nSavedChans"]
        # nsample = meta["nFileSamp"]
        # fs = meta["fs"]
        # nch = nsavedch-1 # exclude sync channel
        # hx,hy,hz = meta["probespacing"]
        # exchmask = exchmasknp(dataset,imecindex=ii,datatype=d)
        # pnrow,pncol = size(exchmask)
        # depths = hy*(0:(pnrow-1))
        # mmbinfile = mmap(binfile,Matrix{Int16},(nsavedch,nsample),0)

        # Prepare AP
        d = "ap$ii"
        meta = dataset[d]["meta"]
        binfile = spike["datapath"]
        chmapraw = spike["chmapraw"]
        nch = spike["nch"]
        nsample = spike["nsample"]
        hx,hy,hz = meta["probespacing"]
        pversion = meta["probeversion"]
        t0 = spike["t0"]
        winvraw = spike["winvraw"]
        exchmask = exchmasknp(dataset,imecindex=ii,datatype=d)
        pnrow,pncol = size(exchmask)
        depths = hy*(0:(pnrow-1))
        mmbinfile = mmap(binfile,Matrix{Int16},(nch,nsample),0)
        synccondon = ref2sync(condon,dataset,ii)
        exenv["hy"] = hy
        exenv["t0"] = t0
        exenv["synccondon"] = synccondon
        exenv["exchmask"] = exchmask
        exenv["chmapraw"] = chmapraw
        exenv["nch"] = nch
        exenv["nsample"] = nsample
        exenv["winvraw"] = winvraw
        exenv["fs"] = spike["fs"]
        exenv["binfile"] = binfile

        # epoch AP
        epoch = [-40 150]
        epochs = synccondon.+t0.+epoch
        # All AP epochs(mapped to concat file), unwhiten, gain corrected(voltage), bandpass filtered,
        # with all channels mapped in the shape of probe where excluded channels are replaced with rand local average
        ys = fill2mask(epochsamplenp(mmbinfile,spike["fs"],epochs,1:nch;meta,bandpass=[300,3000],whiten=winvraw),exchmask,chmap=chmapraw,randreplace=true)
        # 3 fold downsampling(10kHz)
        ys = resample(ys,1//3,dims=3)
        fs = spike["fs"]/3
        baseindex = epoch2sampleindex([0 50],fs) # [-40, 10] ms

        # power spectrum of same depth
        @views pys = ys[.!exchmask,:,:] # exclude and flat channels
        chpos = vcat(chpositionnp(pversion)[.!exchmask]...) # exclude and flat channel positions
        chgi = [findall(chpos[:,2].==up) for up in unique(chpos[:,2])] # group channels with same depth
        @views ppss,psfreqs = powerspectrum(pys[:,baseindex[end]+1:end,:],fs;freqrange=[300,3000])
        pss = Array{Float64}(undef,length(chgi),size(ppss)[2:end]...)
        @views foreach(i->pss[i,:,:] = mean(ppss[chgi[i],:,:],dims=1),eachindex(chgi)) # average of same depth
        tps = dropdims(mean(pss,dims=2),dims=2) # freq average
        @views cfps = Dict(condstring(r)=>
            dropdims(mean(pss[:,:,r.i],dims=3),dims=3) # trial average
        for r in eachrow(cond))

        # power contrast of same depth
        prmss = mapwindow(x->rms(x)^2,pys,(1,101,1),border="symmetric") # 10ms rms window
        rmss = Array{Float64}(undef,length(chgi),size(prmss)[2:end]...)
        @views foreach(i->rmss[i,:,:] = mean(prmss[chgi[i],:,:],dims=1),eachindex(chgi)) # average of same depth
        @views crms = Dict(condstring(r)=>
            stfilter(dropdims(mean(rmss[:,:,r.i],dims=3),dims=3),temporaltype=:rcb,ti=baseindex) # trial average
        for r in eachrow(cond))
        times = range(start=epoch[1],step=1/fs/SecondPerUnit,length=size(pys,2))

        # local coherence
        @views lcs,lcfreqs = localcoherence(pys[:,baseindex[end]+1:end,:],chpos,fs;freqrange=[300,3000],lr=55,sigma=25,chgroupdim=2)
        tlc = dropdims(mean(lcs,dims=2),dims=2) # freq average
        @views cflc = Dict(condstring(r)=>
            dropdims(mean(lcs[:,:,r.i],dims=3),dims=3) # trial average
        for r in eachrow(cond))

        # @views pc,freqs = coherencespectrum(ys[:,baseindex[end]+1:end,:],fs,freqrange=[300,3000],ismean=true) # freq average

        # dsn = 2:5 # downsample n, 40,60,80,100μm for hy=20μm
        # dsr = [3,2,1,1] # downsample r, 120,120,80,100μm for dsn
        # @views pcd = Dict(dsn[d]=>hcat(map(i->bandmean(pc[1:dsn[d]:end,1:dsn[d]:end,i],r=dsr[d],s=1),1:size(pc,3))...) for d in 1:4)
        # @views pc = hcat(map(i->bandmean(pc[:,:,i],r=5,s=1),1:size(pc,3))...) # local average with gaussian weights of -100:100μm, σ=20μm

        if plot
            plotanalog(tps;hy,color=:heat,n=mean(tps,dims=2),xlabel="Trial",xunit="",cunit=:fr,layer)
            foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_tPower$ext")),figfmt)
            plotanalog(tlc;hy,xlabel="Trial",xunit="",color=:heat,cunit=:fr,n=mean(tlc,dims=2),layer)
            foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_tCoherence$ext")),figfmt)
            fpslims = extrema(mapreduce(extrema,union,values(cfps)))
            rmslim = mapreduce(x->maximum(abs.(x)),max,values(crms))
            flclims = extrema(mapreduce(extrema,union,values(cflc)))
            for k in keys(crms)
                plotanalog(cfps[k];hy,x=psfreqs,xlabel="Frequency",xunit="Hz",clims=fpslims,color=:heat,cunit=:fr,n=mean(cfps[k],dims=2),layer)
                foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_$(k)_fPower$ext")),figfmt)
                plotanalog(crms[k];x=times,hy,clims=(-rmslim,rmslim),layer)
                foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_$(k)_dRMS$ext")),figfmt)
                plotanalog(cflc[k];hy,x=lcfreqs,xlabel="Frequency",xunit="Hz",clims=flclims,color=:heat,cunit=:fr,n=mean(cflc[k],dims=2),layer)
                foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_$(k)_fCoherence$ext")),figfmt)
            end
        end
        jldsave(joinpath(resultdir,"ap$(ii).jld2");tps,cfps,crms,tlc,cflc,times,psfreqs,lcfreqs,depths,baseindex,fs,siteid,exenv)


        # # baseline power spectrum of same depth
        # @views ppss,psfreqs = powerspectrum(pys[:,baseindex,:],fs;freqrange=[300,3000])
        # pss = Array{Float64}(undef,length(chgi),size(ppss)[2:end]...)
        # @views foreach(i->pss[i,:,:] = mean(ppss[chgi[i],:,:],dims=1),eachindex(chgi)) # average of same depth
        # tps = dropdims(mean(pss,dims=2),dims=2) # freq average
        # @views cfps = Dict(condstring(r)=>
        #     dropdims(mean(pss[:,:,r.i],dims=3),dims=3) # trial average
        # for r in eachrow(cond))

        # # baseline local coherence
        # @views lcs,lcfreqs = localcoherence(pys[:,baseindex,:],chpos,fs;freqrange=[300,3000],lr=55,sigma=25,chgroupdim=2)
        # tlc = dropdims(mean(lcs,dims=2),dims=2) # freq average
        # @views cflc = Dict(condstring(r)=>
        #     dropdims(mean(lcs[:,:,r.i],dims=3),dims=3) # trial average
        # for r in eachrow(cond))

        # if plot
        #     plotanalog(tps;hy,color=:heat,n=mean(tps,dims=2),xlabel="Trial",xunit="",cunit=:fr,layer)
        #     foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_tbPower$ext")),figfmt)
        #     plotanalog(tlc;hy,xlabel="Trial",xunit="",color=:heat,cunit=:fr,n=mean(tlc,dims=2),layer)
        #     foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_tbCoherence$ext")),figfmt)
        #     fpslims = extrema(mapreduce(extrema,union,values(cfps)))
        #     flclims = extrema(mapreduce(extrema,union,values(cflc)))
        #     for k in keys(cfps)
        #         plotanalog(cfps[k];hy,x=psfreqs,xlabel="Frequency",xunit="Hz",clims=fpslims,color=:heat,cunit=:fr,n=mean(cfps[k],dims=2),layer)
        #         foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_$(k)_fbPower$ext")),figfmt)
        #         plotanalog(cflc[k];hy,x=lcfreqs,xlabel="Frequency",xunit="Hz",clims=flclims,color=:heat,cunit=:fr,n=mean(cflc[k],dims=2),layer)
        #         foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_$(k)_fbCoherence$ext")),figfmt)
        #     end
        # end
        # jldsave(joinpath(resultdir,"ap$(ii)+.jld2");tps,cfps,tlc,cflc,psfreqs,lcfreqs,depths,baseindex,fs,siteid,exenv)

        # # Whiten epoch AP
        # epoch = [-40 150]
        # epochs = synccondon.+t0.+epoch
        # # All AP epochs(mapped to concat file), remain whiten, bandpass filtered,
        # # with all channels mapped in the shape of probe where excluded channels are replaced with rand local average
        # ys = fill2mask(epochsamplenp(mmbinfile,spike["fs"],epochs,1:nch;bandpass=[300,3000]),exchmask,chmap=chmapraw,randreplace=true)
        # # 3 fold downsampling(10kHz)
        # ys = resample(ys,1//3,dims=3)
        # fs = spike["fs"]/3
        # baseindex = epoch2sampleindex([0 50],fs) # [-40, 10] ms

        # # power spectrum of same depth
        # @views pys = ys[.!exchmask,:,:] # exclude and flat channels
        # chpos = vcat(chpositionnp(pversion)[.!exchmask]...) # exclude and flat channel positions
        # chgi = [findall(chpos[:,2].==up) for up in unique(chpos[:,2])] # group channels with same depth
        # @views ppss,psfreqs = powerspectrum(pys[:,baseindex[end]+1:end,:],fs;freqrange=[300,3000])
        # pss = Array{Float64}(undef,length(chgi),size(ppss)[2:end]...)
        # @views foreach(i->pss[i,:,:] = mean(ppss[chgi[i],:,:],dims=1),eachindex(chgi)) # average of same depth
        # wtps = dropdims(mean(pss,dims=2),dims=2) # freq average
        # @views wcfps = Dict(condstring(r)=>
        #     dropdims(mean(pss[:,:,r.i],dims=3),dims=3) # trial average
        # for r in eachrow(cond))

        # # power contrast of same depth
        # prmss = mapwindow(x->rms(x)^2,pys,(1,101,1),border="symmetric") # 10ms rms window
        # rmss = Array{Float64}(undef,length(chgi),size(prmss)[2:end]...)
        # @views foreach(i->rmss[i,:,:] = mean(prmss[chgi[i],:,:],dims=1),eachindex(chgi)) # average of same depth
        # @views wcrms = Dict(condstring(r)=>
        #     stfilter(dropdims(mean(rmss[:,:,r.i],dims=3),dims=3),temporaltype=:rcb,ti=baseindex) # trial average
        # for r in eachrow(cond))
        # times = range(start=epoch[1],step=1/fs/SecondPerUnit,length=size(pys,2))

        # # local coherence
        # @views lcs,lcfreqs = localcoherence(pys[:,baseindex[end]+1:end,:],chpos,fs;freqrange=[300,3000],lr=55,sigma=25,chgroupdim=2)
        # wtlc = dropdims(mean(lcs,dims=2),dims=2) # freq average
        # @views wcflc = Dict(condstring(r)=>
        #     dropdims(mean(lcs[:,:,r.i],dims=3),dims=3) # trial average
        # for r in eachrow(cond))

        # if plot
        #     plotanalog(wtps;hy,color=:heat,n=mean(wtps,dims=2),xlabel="Trial",xunit="",cunit=:fr,layer)
        #     foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_wtPower$ext")),figfmt)
        #     plotanalog(wtlc;hy,xlabel="Trial",xunit="",color=:heat,cunit=:fr,n=mean(wtlc,dims=2),layer)
        #     foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_wtCoherence$ext")),figfmt)
        #     fpslims = extrema(mapreduce(extrema,union,values(wcfps)))
        #     rmslim = mapreduce(x->maximum(abs.(x)),max,values(wcrms))
        #     flclims = extrema(mapreduce(extrema,union,values(wcflc)))
        #     for k in keys(wcrms)
        #         plotanalog(wcfps[k];hy,x=psfreqs,xlabel="Frequency",xunit="Hz",clims=fpslims,color=:heat,cunit=:fr,n=mean(wcfps[k],dims=2),layer)
        #         foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_$(k)_wfPower$ext")),figfmt)
        #         plotanalog(wcrms[k];x=times,hy,clims=(-rmslim,rmslim),layer)
        #         foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_$(k)_wdRMS$ext")),figfmt)
        #         plotanalog(wcflc[k];hy,x=lcfreqs,xlabel="Frequency",xunit="Hz",clims=flclims,color=:heat,cunit=:fr,n=mean(wcflc[k],dims=2),layer)
        #         foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_$(k)_wfCoherence$ext")),figfmt)
        #     end
        # end
        # jldsave(joinpath(resultdir,"wap$(ii).jld2");wtps,wcfps,wcrms,wtlc,wcflc,times,psfreqs,lcfreqs,depths,baseindex,fs,siteid,exenv)


        # # Depth Unit PSTH
        # epoch = [-40 150]
        # epochs = synccondon.+epoch
        # bw = 2 # ms
        # psthbins = epoch[1]:bw:epoch[2]
        # baseindex = epoch2sampleindex([0 50],1/(bw*SecondPerUnit)) # [-40, 10] ms

        # ui = unitsync.==ii
        # @views unitepochpsth = map(st->psthspiketrains(epochspiketrain(st,epochs,isminzero=true,shift=-epoch[1]).y,psthbins,israte=true,ismean=false),unitspike[ui])
        # x = unitepochpsth[1].x
        # @views cpsth = Dict(condstring(r)=>begin
        #                     ucp = map(p->meanse(p.mat[r.i,:]).m,unitepochpsth)
        #                     p = spacepsth(ucp,unitposition[ui,:],dims=2,spacerange=depths,bw=2hy,step=hy)
        #                     (;psth=mapwindow(mean,p.psth,(1,5),border="symmetric"),p.y) # 10ms mean window
        #                 end  for r in eachrow(cond))
        # y = first(values(cpsth)).y
        # cpsth = Dict(k=>cpsth[k].psth for k in keys(cpsth))
        # cdpsth = Dict(k=>stfilter(cpsth[k],temporaltype=:sub,ti=baseindex) for k in keys(cpsth))

        # if plot
        #     plim = mapreduce(i->maximum(abs.(i)),max,values(cdpsth))
        #     for k in keys(cdpsth)
        #         plotanalog(cdpsth[k];x,y,clims=(-plim,plim))
        #         foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_$(k)_dPSTH$ext")),figfmt)
        #     end
        # end
        # jldsave(joinpath(resultdir,"psth$(ii).jld2");cpsth,cdpsth,x,y,baseindex,siteid,exenv)

        # continue
        # Prepare LF
        d = "lf$ii"
        meta = dataset[d]["meta"]
        binfile = meta["fileName"]
        nsavedch = meta["nSavedChans"]
        nsample = meta["nFileSamp"]
        nch = nsavedch-1 # exclude sync channel
        hx,hy,hz = meta["probespacing"]
        pversion = meta["probeversion"]
        t0 = 0 # not concat drift correct binary yet
        exchmask = exchmasknp(dataset,imecindex=ii,datatype=d)
        pnrow,pncol = size(exchmask)
        depths = hy*(0:(pnrow-1))
        mmbinfile = mmap(binfile,Matrix{Int16},(nsavedch,nsample),0)
        synccondon = ref2sync(condon,dataset,ii)
        exenv["hy"] = hy
        exenv["t0"] = t0
        exenv["synccondon"] = synccondon
        exenv["exchmask"] = exchmask
        exenv["nch"] = nch
        exenv["nsample"] = nsample
        exenv["fs"] = meta["fs"]
        exenv["binfile"] = binfile

        # epoch LF
        epoch = [-40 150]
        epochs = synccondon.+t0.+epoch
        # All LF epochs, gain corrected(voltage), line noise(60,120,180Hz) removed, bandpass filtered,
        # with all channels mapped in the shape of probe where excluded channels are replaced with local average
        ys = fill2mask(epochsamplenp(mmbinfile,meta["fs"],epochs,1:nch;meta,bandpass=[1,100]),exchmask)
        # 2.5 fold downsampling(1kHz)
        ys = resample(ys,1/2.5,dims=3)
        fs = meta["fs"]/2.5
        baseindex = epoch2sampleindex([0 50],fs) # [-40, 10] ms

        # @views ys = cat((ys[:,i,:,:] for i in 1:pncol)...,dims=3)
        # @views clfp = Dict(condstring(r)=>
        #     dropdims(mean(ys[:,:,vcat((r.i .+ c*nrow(ctc) for c in 0:pncol-1)...)],dims=3),dims=3)
        #     for r in eachrow(cond))

        # if plot
        #     for c in 1:pncol
        #         cmcys = Dict(condstring(r)=>dropdims(mean(ys[:,c,:,r.i],dims=3),dims=3) for r in eachrow(cond))
        #         cmccsd = Dict(condstring(r)=>dropdims(mean(csd(ys[:,c,:,r.i],h=hy),dims=3),dims=3) for r in eachrow(cond))
        #         for k in keys(cmcys)
        #             plotanalog(cmcys[k];fs,cunit=:uv,hy)
        #             foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_Column_$(c)_$(k)_MeanLFP$ext")),figfmt)

        #             plotanalog(imfilter(cmccsd[k],Kernel.gaussian((1,1)));fs,hy)
        #             foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_Column_$(c)_$(k)_MeanCSD$ext")),figfmt)
        #         end
        #     end
        # end

        # LFP of same depth
        @views pys = ys[.!exchmask,:,:] # exclude and flat channels
        chpos = vcat(chpositionnp(pversion)[.!exchmask]...) # exclude and flat channel positions
        chgi = [findall(chpos[:,2].==up) for up in unique(chpos[:,2])] # group channels with same depth
        lfps = Array{Float64}(undef,length(chgi),size(pys)[2:end]...)
        @views foreach(i->lfps[i,:,:] = mean(pys[chgi[i],:,:],dims=1),eachindex(chgi)) # average of same depth
        @views clfp = Dict(condstring(r)=>
            dropdims(mean(lfps[:,:,r.i],dims=3),dims=3) # trial average
        for r in eachrow(cond))

        # CSD of same depth
        ccsd = Dict(k=>stfilter(csd(v,h=hy),temporaltype=:sub,ti=baseindex) for (k,v) in clfp)
        times = range(start=epoch[1],step=1/fs/SecondPerUnit,length=size(pys,2))

        if plot
            lfplim = mapreduce(x->maximum(abs.(x)),max,values(clfp))
            csdlim = mapreduce(x->maximum(abs.(x)),max,values(ccsd))
            for k in keys(clfp)
                plotanalog(clfp[k];x=times,hy,clims=(-lfplim,lfplim))
                foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_$(k)_LFP$ext")),figfmt)
                plotanalog(ccsd[k];x=times,color=:RdBu,hy,clims=(-csdlim,csdlim))
                foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_$(k)_dCSD$ext")),figfmt)
            end
        end
        jldsave(joinpath(resultdir,"lf$(ii).jld2");clfp,ccsd,times,depths,fs,baseindex,siteid,exenv)


        # # power spectrum of same depth
        # nw = 2
        # @views ppss,psfreqs = powerspectrum(pys[:,baseindex[end]+1:end,:],fs;freqrange=[1,100],nw)
        # pss = Array{Float64}(undef,length(chgi),size(ppss)[2:end]...)
        # @views foreach(i->pss[i,:,:] = mean(ppss[chgi[i],:,:],dims=1),eachindex(chgi)) # average of same depth
        # tps = dropdims(mean(pss,dims=2),dims=2) # freq average
        # @views cfps = Dict(condstring(r)=>
        #     dropdims(mean(pss[:,:,r.i],dims=3),dims=3) # trial average
        # for r in eachrow(cond))

        # # local coherence
        # @views lcs,lcfreqs = localcoherence(pys[:,baseindex[end]+1:end,:],chpos,fs;freqrange=[1,100],lr=55,sigma=25,chgroupdim=2,nw)
        # tlc = dropdims(mean(lcs,dims=2),dims=2) # freq average
        # @views cflc = Dict(condstring(r)=>
        #     dropdims(mean(lcs[:,:,r.i],dims=3),dims=3) # trial average
        # for r in eachrow(cond))

        # if plot
        #     plotanalog(tps;hy,color=:heat,n=mean(tps,dims=2),xlabel="Trial",xunit="",cunit=:fr,layer)
        #     foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_tPower.lf$ext")),figfmt)
        #     plotanalog(tlc;hy,xlabel="Trial",xunit="",color=:heat,cunit=:fr,n=mean(tlc,dims=2),layer)
        #     foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_tCoherence.lf$ext")),figfmt)
        #     fpslims = extrema(mapreduce(extrema,union,values(cfps)))
        #     flclims = extrema(mapreduce(extrema,union,values(cflc)))
        #     for k in keys(cfps)
        #         plotanalog(cfps[k];hy,x=psfreqs,xlabel="Frequency",xunit="Hz",clims=fpslims,color=:heat,cunit=:fr,n=mean(cfps[k],dims=2),layer)
        #         foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_$(k)_fPower.lf$ext")),figfmt)
        #         plotanalog(cflc[k];hy,x=lcfreqs,xlabel="Frequency",xunit="Hz",clims=flclims,color=:heat,cunit=:fr,n=mean(cflc[k],dims=2),layer)
        #         foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_$(k)_fCoherence.lf$ext")),figfmt)
        #     end
        # end
        # jldsave(joinpath(resultdir,"lf$(ii)+.jld2");tps,cfps,psfreqs,tlc,cflc,lcfreqs,depths,fs,baseindex,siteid,exenv)
    end

    # Unit Spike Trian of Epochs
    if plot
        epochext = max(preicidur,suficidur)
        epochs = [condon.-epochext condoff.+epochext]
        for i in eachindex(unitspike)
            ys,ns,ws,is = epochspiketrain(unitspike[i],ref2sync(epochs,dataset,unitsync[i]),isminzero=true,shift=epochext)
            title = "IMEC$(unitsync[i])_$(unitgood[i] ? "S" : "M")U$(unitid[i])_SpikeTrian"
            plotspiketrain(ys,timeline=[0,conddur],title=title)
            foreach(ext->savefig(joinpath(resultdir,"$title$ext")),figfmt)
        end
    end
end
