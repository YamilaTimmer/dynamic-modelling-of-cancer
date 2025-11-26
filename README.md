## Ordinary Differential Equations (ODE)    
This tumor growth simulation contains various types of ordinary differential equations (or ODE) that can be used to describe tumor growth over time. The models vary in complexity and realism and each have their own use-cases, below we elaborate further on each model.  
### Exponential model  
Exponential models are similar to linear models, with the difference here being that the volume of the tumor is multiplied with the growth rate, instead of adding them together. This makes for an exponential increase, which can be described with: $$V_{\text{t}} = V + c \cdot V$$  
  
This model can be used to describe size increase when a tumor is still in the early stages. However, the longer a tumor persists, the more the doubling time increases. This has various underlying factors such as an increase in the average cell-cycle time or loss of cells due to apoptosis [(Gerlee, 2013)](https://doi.org/10.1158/0008-5472.can-12-4355). This model is thus not suited to describe tumor growth over long tem periods.  
  
  ---
### Exponential flattening model  
Realistically, a tumor can not keep growing indefinitely as there are physical and physiological limitations. Some of the models described on this page do not keep this limitation in mind, and allow for the tumor to "grow" indefinitely. In an organism there will always be a maximum volume for the tumor, which depends on factors such as the tumor's access to resources and to "free" space to grow in [(Murphy et al., 2016)](https://doi.org/10.1186/s12885-016-2164-x). The exponential flattening model does keep this limitation in mind and the tumor growth per time unit is described as:  $$V_{\text{t}} = V + c \cdot V \cdot (1.0 - V/v_{max})$$ Here the increase in volume is determined by how close the tumor is to reaching its maximum volume, the closer it gets to this threshold, the more growth will slow down.  
  
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
### Montroll model  

---
### Lineair limited growth model

# Sources  
- Murphy, H., Jaafari, H., & Dobrovolny, H. M. (2016). Differences in predictions of ODE models of tumor growth: a cautionary example. _BMC Cancer_, _16_(1), 163. https://doi.org/10.1186/s12885-016-2164-x  
- Heesterman, B., Bokhorst, J., De Pont, L., Verbist, B., Bayley, J., Van Der Mey, A., Corssmit, E., Hes, F., Van Benthem, P., & Jansen, J. (2018). Mathematical Models for Tumor Growth and the Reduction of Overtreatment. _Journal Of Neurological Surgery Part B Skull Base_, _80_(01), 072–078. https://doi.org/10.1055/s-0038-1667148  
- Chan, K., Kao, C., Gordinier, J., & Ganden, K. (2023). Treatment Optimization for Tumor Growth by Ordinary Differential Equations. _Journal Of Student Research_, _12_(4). https://doi.org/10.47611/jsrhs.v12i4.5202  
- Botmann, N. K. G., & Dobrovolny, H. M. (2025). Assessing the role of model choice in parameter identifiability of cancer treatment efficacy. _Frontiers in Applied Mathematics And Statistics_, _11_. https://doi.org/10.3389/fams.2025.1542617  
- Zanotti, Y. P., & Faria, H. A. M. (2025). CLASSICAL MATHEMATICAL MODELS OF POPULATION GROWTH FOR PREDICTING CELL CULTURE IN BIOREACTORS. In _Seven Editora eBooks_. https://doi.org/10.56238/edimpacto2025.041-005  
- Gerlee, P. (2013). The Model Muddle: In Search of Tumor Growth Laws. _Cancer Research_, _73_(8), 2407–2411. https://doi.org/10.1158/0008-5472.can-12-4355