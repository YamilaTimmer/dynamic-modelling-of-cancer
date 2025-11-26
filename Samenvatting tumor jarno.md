
# Tumor Growth Models Summary

## 1. Linear Growth
* **Formula:** $V(t) = C \cdot t + V_{0}$
* **Rate of Change:** $\frac{dV}{dt} = C$
* **Concept:** This suggests that the same amount of growth occurs each day. This makes it a simple formula to understand and it's implementation will be simple too.
* **Critique:** In reality, tumor cells grow faster when there are more cells. Therefore, this model is not realistic.  This model does not take that into account, also it does not look for nutrient availability nor bloodvessels.

---

## 2. Mendelsohn Growth (Power Law)
* **Formula:** $\frac{dV}{dt} = c \cdot V^{d}$
* **Concept:** Describes "unlimited growth".   This model implies growth continues forever.
* **Critique:** In reality, growth slows at a certain point due to nutrient and vascular (blood vessel) limitations.

---

## 3. Logistic Growth (1838 Model)
* **Formula:** $\frac{dV}{dt} = C \cdot V \cdot (V_{max} - V)$
* **Concept:**  Created in 1838 to simulate realistic conditions.Takes into account that growth will stagnate at some point.
* **Mechanism:** Stagnation occurs due to a lack of nutrients. 
* **Critique:** May not be realistic for every tumour type 

---

## 4. Allee Effect
* **Formula:** $\frac{dV}{dt} = C \cdot (V - V_{min})$
* **Concept:** States that the growth rate of a tumor can be dependent on its cell population. This means that if there is little cells the tumour might not grow but even decline and therefor the tumour could die.

---

## 5. Surface / 2-Dim Growth
* **Context:** "Oppervlakte" (Surface Area).
* **Mechanism:** Occurs very early in development, without blood vessels. Growth is dependent on $O_2$ (oxygen) and nutrients diffusing from the flesh around it. Because the nutrients get diffused through said cell surfrace, the growth will be limited by the surface area and not the volume.

---

## 6. General Limitations: Noise
* **Issue:** Growth models generally do not take "noise" into account.
* **Consequence:** If combining multiple scans, this lack of noise accounting causes problems.


## Code

```python
# possible upset
class LinearGrowth
# growth rate

def __init__(self, growthRate, startvolume, time):
	self.C = growthRate
	self.V0 = startvolume
	self.time = time
	
	

def LinearGrowth(self):
	ResultLinGrow = []
	for i in self.time:
		V = self.C * self.time + self.V0
		ResultLinGrow.append(V)
		
	return ResultLinGrow
	
def MendelsohnGrowth(self):
	ResultMendelGrowth = []
	# nog niet de wiskunde manier gevonden
	for i in self.time:
		V = 

def Logistic Growth(self)
	# not correct yet
	ResultLogicGrowth = []
	V = 0
	for i in self.time:
		V = self.C * V * (Vmax - V)
		ResultLogicGrowth.append(V)

```