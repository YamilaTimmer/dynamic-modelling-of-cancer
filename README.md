# Dynamic modelling of cancer using different growth models

## Authors

---

- Yamila Timmer
- Jarno Duiker

## Relevance

---

## Description

---

In this package contains different growth models that all have their flaws and plusses. The 
package is designed for dynamic modeling of cancer or tumors. 


## Ordinary Differential Equations (ODE)    
This tumor growth simulation contains various types of ordinary differential equations (or ODE) that can be used to describe tumor growth over time. The models vary in complexity and realism and each have their own use-cases, below we elaborate further on each model.  

### Linear Growth model
Linear models suggest that there will be a stable daily growth, this makes it a simple model to use.
This model has a great flaw, as tumor cells start growing more when the larger they become, they only are limited by nutrion,
oxygen and bloodvessels. Therefore this model is not to be seen as accurate.

The mathematical description of this model is as follows:  $V(t) = C \cdot t + V_{0}$ where C = per day growth rate. 
The graph that will come out of this function is a straight line as seen in figure 1.


*Figure 1: Linear growth model*

### Lineair limited growth model
This model treat the tumour like 2 fases, when the tumour is grows very slowly and the bigger it gets the more constant  the growth speed gets (turning linear)
It will look like an exponential growth after which it turns into a more linear growth as to be seen in figure x.


the formula is: $$ \frac{dV}{dt} = c \cdot \frac{V}{V + d} $$

![linearlim.png](Img%2Flinearlim.png)
---

### Exponential model  
Exponential models are similar to linear models, with the difference here being that the volume of the tumor is multiplied with the growth rate, instead of adding them together. This makes for an exponential increase, which can be described with: $$V_{\text{t}} = V + c \cdot V$$ , The resulting plot can be seen in figure 2.
  
This model can be used to describe size increase when a tumor is still in the early stages. However, the longer a tumor persists, the more the doubling time increases. This has various underlying factors such as an increase in the average cell-cycle time or loss of cells due to apoptosis [(Gerlee, 2013)](https://doi.org/10.1158/0008-5472.can-12-4355). This model is thus not suited to describe tumor growth over long tem periods.  

![exponentialgrowth.png](Img%2Fexponentialgrowth.png)

*Figure 2: Exponential model*

  ---
### Exponential flattening model  
Realistically, a tumor can not keep growing indefinitely as there are physical and physiological limitations. Some of the models described on this page do not keep this limitation in mind, and allow for the tumor to "grow" indefinitely. In an organism there will always be a maximum volume for the tumor, which depends on factors such as the tumor's access to resources and to "free" space to grow in [(Murphy et al., 2016)](https://doi.org/10.1186/s12885-016-2164-x). The exponential flattening model does keep this limitation in mind and the tumor growth per time unit is described as:  $$V_{\text{t}} = V + c \cdot V \cdot (1.0 - V/v_{max})$$ Here the increase in volume is determined by how close the tumor is to reaching its maximum volume, the closer it gets to this threshold, the more growth will slow down.

![exponentailflat.png](Img%2Fexponentailflat.png)

*Figure 3: Exponential flattening model*

---
### Mendelsohn Growth 
In contrast to the Exponential flattening model, the Mendelsohn growth model believes that growth is umilited (Bindhammer, n.d.) . therefore implying that a tumor cell could replicate and grow bigger forever.
In reality however, tumor growth does slow down due to the host not being able to provide some of it's crucial growth needs. These being oxygen, nutrients and bloodvessels (providing oxygen and nutrients).

Mendelsohn proposed the formula: $\frac{dV}{dt} = c \cdot V^{d}$ (Gerlee, 2013) c is the tumour growth rate, this gets multiplied by the volume to the power of d (The allometric factor). This gives a plot like in figure 4.

![mendelsohn.png](Img%2Fmendelsohn.png)

*Figure 4: Mendelsohn Growth model*

---
### Von Bertalanffy Model  
This ODE takes into account growth in relation to surface area (surface rule model), which states cell growth should be proportional to its surface area [(Chan et al., 2023)](https://doi.org/10.47611/jsrhs.v12i4.5202). The model assumes that the net growth rate does not only consist of tumor cell proliferation but also of tumor cell death [(Botmann & Dobrovolny, 2025)](https://doi.org/10.3389/fams.2025.1542617). In order to apply this, in the equation, volume is powered to 2/3. The bigger the surface area, the higher the amount of nutrients/energy for the cells to absorb and the faster the tumor will grow.   
  
This model has been succesfully implemented to predict tumor growth in literature, e.g. by [Heesterman et al. (2018)](https://doi.org/10.1055/s-0038-1667148  ) and can be described using:  
$$Vt = c \cdot V^{2/3} - V_{max} \cdot V$$

![bertalanfyy.png](Img%2Fbertalanfyy.png)

---  
  
### Gompertz model  
The Gompertz model was originally made to predict human mortality curves, but turned out to be a very suitable model to predict cancer growth, as it seems to provide the best predictions for e.g. breast and lung cancer growth [(Murphy et al., 2016)](https://doi.org/10.1186/s12885-016-2164-x). This model takes into account growth velocity, or the change of weight/height over time which is useful for monitoring growth [(Zanotti & Faria, 2025)](https://doi.org/10.56238/edimpacto2025.041-005).   
  
Per time unit, this model is described with: $$V_t = c \cdot V \cdot \ln\Big(\frac{V_{max}}{V}\Big)$$  

![gompert.png](Img%2Fgompert.png)

---
###  Logistic Growth 
The Exponentail Growth model has some limitations in predicting the long term growth rate of cancer cell proliferation (Tabassum et al., 2019).
This is the reason the logistic model was introduced. It explains the behaviour of a later stage cancer cell better. This model was first introduced
in 1883 by Francois verhulst.

The version we will implement is the following:
$\frac{dV}{dt} = C \cdot V \cdot (V_{max} - V)$

The max tumour vollume occurs due to the lack of bloodvessels and fight over nutrients (Chan et al., 2023). 

![loggrowth.png](Img%2Floggrowth.png)

### Montroll growth model  
A model that says tumors can grow quickly when there is nutrients, space and bloodvessels but the growth will slow when these resources become scares (Rodrigues, 2024)(Goel et al., 1971).
It will look like a log function, because the growth will quickly rise and then level out.


The formula for this model is: $ \frac{dV}{dt} = c \cdot V \cdot \left(\frac{V}{d_{\text{max}} - V}\right) $

The graph will look like this:
![montrol.png](Img%2Fmontrol.png)

---

---
### Allee Effect
Models in tumor growth often assume exponential growth kinetics at low cell population. In recent pre & clinical observations of tumor initiation
or recurrence indicate the presence of tumor growth kinetics in which growth rates scale positvely with cell numbers.
Recent observations however suggest a cooperative growth pattern also known as the allee effect. Here growth rates increase with cell numbers at low densities (Johnson et al., 2019).

The formula for this model is: $\frac{dV}{dt} = C \cdot (V - V_{min})$

This wil result in a graph like in Figure x
![alleegrowth.png](Img%2Falleegrowth.png)

---


### Surface limited growth model
This model states that when essential nutrient of any kind gets limited, the growth rate of an individual cell will be proportional to its surface area.
rather than to its volume. The decrease in dimensionality from volume to surface is expected to favor the smaller cells.

The model deals with cell growth under unusual nutritional conditions, and the predicitons on how cell replication cycle is assiumed to behave when there is unusual nutritional conditions (Grover, 1988). 

$$\frac{dV}{dt} = c \cdot \frac{V}{(V + d)^{1/3}}$$

![surfacelim.png](Img%2Fsurfacelim.png)
---

# Sources  
- Murphy, H., Jaafari, H., & Dobrovolny, H. M. (2016). Differences in predictions of ODE models of tumor growth: a cautionary example. _BMC Cancer_, _16_(1), 163. https://doi.org/10.1186/s12885-016-2164-x  
- Heesterman, B., Bokhorst, J., De Pont, L., Verbist, B., Bayley, J., Van Der Mey, A., Corssmit, E., Hes, F., Van Benthem, P., & Jansen, J. (2018). Mathematical Models for Tumor Growth and the Reduction of Overtreatment. _Journal Of Neurological Surgery Part B Skull Base_, _80_(01), 072–078. https://doi.org/10.1055/s-0038-1667148  
- Chan, K., Kao, C., Gordinier, J., & Ganden, K. (2023). Treatment Optimization for Tumor Growth by Ordinary Differential Equations. _Journal Of Student Research_, _12_(4). https://doi.org/10.47611/jsrhs.v12i4.5202  
- Botmann, N. K. G., & Dobrovolny, H. M. (2025). Assessing the role of model choice in parameter identifiability of cancer treatment efficacy. _Frontiers in Applied Mathematics And Statistics_, _11_. https://doi.org/10.3389/fams.2025.1542617  
- Zanotti, Y. P., & Faria, H. A. M. (2025). CLASSICAL MATHEMATICAL MODELS OF POPULATION GROWTH FOR PREDICTING CELL CULTURE IN BIOREACTORS. In _Seven Editora eBooks_. https://doi.org/10.56238/edimpacto2025.041-005  
- Gerlee, P. (2013). The Model Muddle: In Search of Tumor Growth Laws. _Cancer Research_, _73_(8), 2407–2411. https://doi.org/10.1158/0008-5472.can-12-4355
- Chan, K., Kao, C.-Y., Gordinier, J., & Ganden, K. (2023). Treatment Optimization for Tumor Growth by Ordinary Differential Equations. Journal of Student Research, 12(4). https://doi.org/10.47611/jsrhs.v12i4.5202
- Bindhammer, M. (n.d.). From the Mendelsohn model to the Gompertz and logistic growth law.
- Talkington, A., & Durrett, R. (2015). Estimating Tumor Growth Rates In Vivo. Bulletin of Mathematical Biology, 77(10), 1934–1954. https://doi.org/10.1007/s11538-015-0110-8
- Tabassum, S., Rosli, N. B., & Binti Mazalan, M. S. A. (2019). Mathematical Modeling of Cancer Growth Process: A Review. Journal of Physics: Conference Series, 1366(1), 012018. https://doi.org/10.1088/1742-6596/1366/1/012018
- Johnson, K. E., Howard, G., Mo, W., Strasser, M. K., Lima, E. A. B. F., Huang, S., & Brock, A. (2019). Cancer cell population growth kinetics at low densities deviate from the exponential growth model and suggest an Allee effect. PLoS Biology, 17(8), e3000399. https://doi.org/10.1371/journal.pbio.3000399
- Rodrigues, J. A. (2024). Using Physics-Informed Neural Networks (PINNs) for Tumor Cell Growth Modeling. Mathematics, 12(8), 1195. https://doi.org/10.3390/math12081195
- Goel, N. S., Maitra, S. C., & Montroll, E. W. (1971). On the Volterra and Other Nonlinear Models of Interacting Populations. Reviews Of Modern Physics, 43(2), 231–276. https://doi.org/10.1103/revmodphys.43.231
- Grover, N. (1988). Surface-limited growth: A model for the synchronization of a growing bacterial culture through periodic starvation. Journal Of Theoretical Biology, 134(1), 77–87. https://doi.org/10.1016/s0022-5193(88)80303-5