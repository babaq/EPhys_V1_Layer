# EPhys_V1_Layer

1. Each `Flash2Color` test is processed through batching in **"BatchTests_SpikeGLX.jl"**, during which AP and LF metrics result is saved in `ap0.jld2` and `lf0.jld2`, and spike info is saved in `spike.jld2`.
1. Spike metrics is processed in **"SpikeInfo.jl"**, during which result is saved in `unitdepthfeature.jld2`.
1. Manual layer assignment, layer normalization and statistics in **"LayerEstimation.jl"**.
