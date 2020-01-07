# Circadian Simulation

|Time Series|Phase Plot|
|---|---|
|<img src="img/anim_timeSeries.gif">|<img src="img/anim_phasePlot.gif">|

Author: Shao-Ting Chiu (stevengogogo4321@gmail.com), 丁啟祐

---
## Abstract

The circadian rhythm is a biological oscillator that behaves limit cycle and is synchronized with solar time. In this presentation, we simulates the circadian model(Goldbeter,1996) to describe the modulation of circadian clock by changing the degradation rate of PER protein. 

## Documentation

- <a href="doc/CircadianClock.pdf">PDF slides</a>

## Outline




<center>
<img src="img/MODEL_circadian_goldbeter1996.png">
</center>

*Source: Ingalls, B. P. (2013).*

<img src="img/math_gold1996.png">

## Simulation

### Circadian clock doesn't always oscillate under any degradation rate of PER

|Time Series|Phase Plot|
|---|---|
|<img src="img/timeSeries_three_0.2.png">|<img src="img/phasePlot_criticalpoint_0.2.png">|
|<img src="img/timeSeries_three_0.6.png">|<img src="img/phasePlot_criticalpoint_0.6.png">|
|<img src="img/timeSeries_three_0.95.png">|<img src="img/phasePlot_criticalpoint_0.95.png">|
|<img src="img/timeSeries_three_1.9.png">|<img src="img/phasePlot_criticalpoint_1.9.png">|
|<img src="img/timeSeries_three_3.2.png">|<img src="img/phasePlot_criticalpoint_3.2.png">|


### Negative Feedback between Nucleus PER and mRNA

<img src="img/animTime_movingCircleD.gif">

### Delay and Cascading Phosphorylation

<img src="img/animTime_movingCircle_p012.gif">

## Scripts

- Mathematical model of PER model
    - <a href="MODEL_PerGoldbeter1996.jl">MODEL_PerGoldbeter1996.jl</a>
- Plotting
    - <a href="plot_PERGoldbeter1996.jl">plot_PERGoldbeter1996.jl</a>
- Interactive Plotting
    - <a href="/Interactive_Plotting.jl">Interactive_Plotting.jl</a>

## Installation

The model is implemented in [the Julia language v1.3.0](https://julialang.org/).. Please install the following packages by [Pkg module](https://docs.julialang.org/en/v1/stdlib/Pkg/index.html):

- ModelingTookkit
- DifferentialEquations
- Plots

## References
1. Ingalls, B. P. (2013). Mathematical modeling in systems biology: an introduction. MIT press. ([link](https://books.google.com.tw/books?hl=zh-TW&lr=&id=OYr6AQAAQBAJ&oi=fnd&pg=PR5&dq=Ingalls,+B.+Mathematical+Modeling+in+Systems+Biology+-+an+Introduction.&ots=ucgsG0-NAA&sig=gXJPRtpiAQDzyLYYcGOT5CQBFVc&redir_esc=y#v=onepage&q=Ingalls%2C%20B.%20Mathematical%20Modeling%20in%20Systems%20Biology%20-%20an%20Introduction.&f=false))
1. Goldbeter, A. (1997). Biochemical oscillations and cellular rhythms: the molecular bases of periodic and chaotic behaviour. Cambridge university press. ([link](https://books.google.com.tw/books?hl=en&lr=&id=dKk0I-KMDJIC&oi=fnd&pg=PP1&ots=WVtd4X2-1N&sig=9pWRpEfrLnXo7kqvaTOfQBCpMUU&redir_esc=y#v=onepage&q&f=false))
