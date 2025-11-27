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

![lingrow.png](Img%2Flingrow.png)

*Figure 1: Linear growth model*

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

---  
  
### Gompertz model  
The Gompertz model was originally made to predict human mortality curves, but turned out to be a very suitable model to predict cancer growth, as it seems to provide the best predictions for e.g. breast and lung cancer growth [(Murphy et al., 2016)](https://doi.org/10.1186/s12885-016-2164-x). This model takes into account growth velocity, or the change of weight/height over time which is useful for monitoring growth [(Zanotti & Faria, 2025)](https://doi.org/10.56238/edimpacto2025.041-005).   
  
Per time unit, this model is described with: $$V_t = c \cdot V \cdot \ln\Big(\frac{V_{max}}{V}\Big)$$  

---
###  Logistic Growth 
The Exponentail Growth model has some limitations in predicting the long term growth rate of cancer cell proliferation (Tabassum et al., 2019).
This is the reason the logistic model was introduced. It explains the behaviour of a later stage cancer cell better. This model was first introduced
in 1883 by Francois verhulst.

The version we will implement is the following:
$\frac{dV}{dt} = C \cdot V \cdot (V_{max} - V)$

The max tumour vollume occurs due to the lack of bloodvessels and fight over nutrients (Chan et al., 2023). 


---
### Allee Effect
Models in tumor growth often assume exponential growth kinetics at low cell population. In recent pre & clinical observations of tumor initiation
or recurrence indicate the presence of tumor growth kinetics in which growth rates scale positvely with cell numbers.
Recent observations however suggest a cooperative growth pattern also known as the allee effect. Here growth rates increase with cell numbers at low densities (Johnson et al., 2019).

The formula for this model is: $\frac{dV}{dt} = C \cdot (V - V_{min})$

---
### Montroll model  

---
### Lineair limited growth model

$$ \frac{dV}{dt} = c \cdot \frac{V}{V + d} $$

---
### Surface limited growth model

$$\frac{dV}{dt} = c \cdot \frac{V}{(V + d)^{1/3}}$$

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
